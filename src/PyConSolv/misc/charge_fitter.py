"""Constrained least-squares charge fitter for FFTK-style water-interaction
fitting.

Given a batch of QM-derived target interaction energies (one per
ligand+water probe) and the ligand's CGenFF Lennard-Jones parameters,
finds the per-atom partial charges that best reproduce the QM values
under:

  - an equality constraint on total charge
  - optional symmetry groups (atoms inside a group share one charge)
  - per-atom magnitude bounds
  - a small L2 regulariser pulling each atom's charge toward its
    initial (CGenFF) value, to suppress wild solutions when the data
    is information-poor

The MM interaction energy used as the objective's prediction is the
ligand <-> water cross-term: pairwise Coulomb plus a 12-6 Lennard-Jones
in CHARMM convention (Lorentz-Berthelot combining; epsilon stored as a
negative number in PRM files but used as |epsilon| in the kernel).
Intra-molecular ligand and intra-water terms cancel out exactly when
subtracting E(ligand) + E(water) from E(complex), so they are not
included in the evaluator.
"""
from dataclasses import dataclass
from typing import List, Optional

import numpy as np
from scipy.optimize import minimize


COULOMB_K = 332.0716   # kcal A / (mol e^2), CHARMM/AMBER convention

TIP3P_CHARGES = np.array([-0.834, 0.417, 0.417])
TIP3P_LJ = np.array([
    [-0.1521, 1.7682],   # O
    [-0.0460, 0.2245],   # H
    [-0.0460, 0.2245],   # H
])


@dataclass
class FitResult:
    charges: np.ndarray
    initial_charges: np.ndarray
    qm_targets: np.ndarray
    mm_predictions: np.ndarray
    residuals: np.ndarray
    rmsd: float
    objective: float
    n_iter: int
    success: bool
    message: str


def mmCrossEnergy(coords_A: np.ndarray, charges_A: np.ndarray,
                   lj_A: np.ndarray,
                   coords_B: np.ndarray, charges_B: np.ndarray,
                   lj_B: np.ndarray) -> float:
    """Pairwise Coulomb + 12-6 LJ between two atom sets.

    coords_*  : (N, 3) Angstroms
    charges_* : (N,) in elementary charge units
    lj_*      : (N, 2) per-atom [epsilon, Rmin/2] in CHARMM convention.
                Epsilon is taken as |epsilon|; PRM files write it as a
                negative number by convention but the kernel uses the
                magnitude.

    Returns interaction energy in kcal/mol.
    """
    diffs = coords_A[:, None, :] - coords_B[None, :, :]
    r = np.linalg.norm(diffs, axis=-1)

    e_coul = COULOMB_K * np.sum(np.outer(charges_A, charges_B) / r)

    eps_ij = np.sqrt(np.outer(np.abs(lj_A[:, 0]), np.abs(lj_B[:, 0])))
    rmin_ij = lj_A[:, 1, None] + lj_B[None, :, 1]
    ratio6 = (rmin_ij / r) ** 6
    e_lj = np.sum(eps_ij * (ratio6 * ratio6 - 2.0 * ratio6))

    return float(e_coul + e_lj)


def _buildExpansion(n_atoms: int, symmetry_groups: list) -> tuple:
    """Return (independent_indices, group_sizes, expand_fn).

    expand_fn maps a reduced parameter vector
       x = [group_0_charge, ..., group_{G-1}_charge,
            indep_0_charge, ..., indep_{I-1}_charge]
    to the full per-atom charge array.
    """
    symmetry_groups = symmetry_groups or []
    grouped = set()
    for g in symmetry_groups:
        if any(i in grouped for i in g):
            raise ValueError('overlapping symmetry groups')
        grouped.update(g)
    independent = [i for i in range(n_atoms) if i not in grouped]
    group_sizes = np.array([len(g) for g in symmetry_groups], dtype=float)
    n_groups = len(symmetry_groups)

    def expand(x: np.ndarray) -> np.ndarray:
        q = np.zeros(n_atoms)
        for k, group in enumerate(symmetry_groups):
            for i in group:
                q[i] = x[k]
        for k, i in enumerate(independent):
            q[i] = x[n_groups + k]
        return q

    return independent, group_sizes, expand


def _initialReduced(initial_charges: np.ndarray,
                     symmetry_groups: list,
                     independent: list) -> np.ndarray:
    symmetry_groups = symmetry_groups or []
    n_groups = len(symmetry_groups)
    x0 = np.zeros(n_groups + len(independent))
    for k, g in enumerate(symmetry_groups):
        x0[k] = float(np.mean(initial_charges[list(g)]))
    for k, i in enumerate(independent):
        x0[n_groups + k] = float(initial_charges[i])
    return x0


def fitCharges(
    qm_results: list,
    ligand_coords: np.ndarray,
    ligand_elements: list,
    initial_charges: np.ndarray,
    ligand_lj: np.ndarray,
    water_charges: np.ndarray = None,
    water_lj: np.ndarray = None,
    total_charge: float = 0.0,
    symmetry_groups: list = None,
    regularizer: float = 1e-1,
    weights: list = None,
    charge_bounds: tuple = (-2.0, 2.0),
    method: str = 'SLSQP',
    max_iter: int = 200,
    ftol: float = 1e-8,
) -> FitResult:
    """Fit ligand charges to FFTK-scaled QM interaction energies.

    qm_results       : list of InteractionResult (or duck-typed objects
                       with .geometry, .n_ligand_atoms, .e_int_scaled).
    ligand_coords    : (N, 3) ligand-only Cartesian, Angstroms.
    initial_charges  : (N,) starting CGenFF charges.
    ligand_lj        : (N, 2) per-atom [epsilon, Rmin/2] in CHARMM units.
    water_charges    : (3,) water charges; defaults to TIP3P.
    water_lj         : (3, 2) water LJ; defaults to TIP3P.
    total_charge     : enforced sum(q).
    symmetry_groups  : list of index lists; atoms in a group share one q.
    regularizer      : L2 weight pulling each q toward initial_charges.
                       Default 0.1 is tuned for real probes covering one or
                       a few polar sites; under-determined charges (atoms
                       far from any probe) stay near CGenFF. Drop to ~1e-3
                       only when probes cover the whole molecule.
    weights          : per-probe weights for the data term (default 1.0).
    charge_bounds    : (lo, hi) per-atom magnitude bounds.
    method           : scipy optimizer; SLSQP supports eq constraints.
    """
    initial_charges = np.asarray(initial_charges, dtype=float)
    ligand_lj = np.asarray(ligand_lj, dtype=float)
    water_charges = (np.asarray(water_charges, dtype=float)
                     if water_charges is not None else TIP3P_CHARGES)
    water_lj = (np.asarray(water_lj, dtype=float)
                if water_lj is not None else TIP3P_LJ)

    n_atoms = len(initial_charges)
    independent, group_sizes, expand = _buildExpansion(n_atoms,
                                                          symmetry_groups)
    x0 = _initialReduced(initial_charges, symmetry_groups, independent)

    weights = (np.asarray(weights, dtype=float) if weights is not None
               else np.ones(len(qm_results)))
    targets = np.array([r.e_int_scaled for r in qm_results])

    def predict(q: np.ndarray) -> np.ndarray:
        out = np.empty(len(qm_results))
        for k, r in enumerate(qm_results):
            n_lig = r.n_ligand_atoms
            water_xyz = r.geometry[n_lig:n_lig + 3]
            out[k] = mmCrossEnergy(ligand_coords, q, ligand_lj,
                                   water_xyz, water_charges, water_lj)
        return out

    def objective(x: np.ndarray) -> float:
        q = expand(x)
        e_mm = predict(q)
        data = float(np.sum(weights * (targets - e_mm) ** 2))
        reg = float(regularizer * np.sum((q - initial_charges) ** 2))
        return data + reg

    n_groups = len(symmetry_groups) if symmetry_groups else 0

    def total_charge_constraint(x: np.ndarray) -> float:
        s = 0.0
        if n_groups:
            s += float(np.sum(x[:n_groups] * group_sizes))
        s += float(np.sum(x[n_groups:]))
        return s - total_charge

    bounds = [(float(charge_bounds[0]), float(charge_bounds[1]))
              for _ in range(len(x0))]
    constraints = [{'type': 'eq', 'fun': total_charge_constraint}]

    result = minimize(
        objective, x0, method=method, bounds=bounds,
        constraints=constraints,
        options={'maxiter': max_iter, 'ftol': ftol},
    )

    q_final = expand(result.x)
    predictions = predict(q_final)
    residuals = targets - predictions
    rmsd = float(np.sqrt(np.mean(residuals ** 2))) if len(residuals) else 0.0

    return FitResult(
        charges=q_final,
        initial_charges=initial_charges.copy(),
        qm_targets=targets,
        mm_predictions=predictions,
        residuals=residuals,
        rmsd=rmsd,
        objective=float(result.fun),
        n_iter=int(getattr(result, 'nit', 0)),
        success=bool(result.success),
        message=str(getattr(result, 'message', '')),
    )
