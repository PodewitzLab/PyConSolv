import os
import matplotlib.pyplot as plt
import numpy as np
from tkinter import *
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

from .inputparser import XYZ


class GUI:
    def __init__(self, window, path: str, files: list):
        """
        Class to present a GUI to the user where charges for each fragment can be entered

        Parameters:
            :param window: TKinter window
            :param string path: path to pdb files of fragments
            :param list files: names of files for which charges need to be assigned

        Class variables:
            - self.path = path to pdb files of fragments
            - self.files = names of files for which charges need to be assigned
            - self.window = TKinter window
            - self.box = text entry box
            - self.button = ok button to proceed to next structure
            - self.structureid = stores the index of the current structure from the list of files
            - self.charges = stores value of charges for all fragments
        """
        self.path = path
        self.files = files
        self.window = window
        self.box = Entry(window)
        self.box.insert(0, "-1")
        self.box.pack()
        self.button = Button(window, text="Assign Charge", command=self.getValue)
        self.button.pack()
        self.plot(self.path + '/' + self.files[0])
        self.structureid = 0
        self.charges = np.zeros(len(self.files))

    def plot(self, file: str):
        """
        Create and display an image in the TKinter window

        Parameters:
            :param string file: pdb file to be displayed

        Class variables:
            - self.canvas = tkinter canvas that displays the image
        """
        fig = plt.figure(figsize=(6, 6))
        ax = plt.gca()
        ax.imshow(self.getPDBasImage(file))
        ax.axis('off')
        self.canvas = FigureCanvasTkAgg(fig, master=self.window)
        self.canvas.get_tk_widget().pack()
        self.canvas.draw()

    def redraw(self, file: str):
        """
        Redraw image on canvas

        Parameters:
            :param string file: pdb file to be displayed

        Class variables:
        """
        self.box.delete(0, END)
        self.box.insert(0, "-1")
        ax = plt.gca()
        ax.imshow(self.getPDBasImage(file))
        ax.axis('off')
        self.canvas.draw()

    def getValue(self):
        """
        Get value entered in the entry box. If the last structure has a charge assigned, it closes the window
        and writes the charge map file.

        Parameters:

        Class variables:
        """
        self.charges[self.structureid] = self.box.get()
        if self.structureid < (len(self.files) - 1):
            self.structureid += 1
            self.redraw(self.path + '/' + self.files[self.structureid])
        else:
            self.writeChargeMap()
            self.quit()

    def writeChargeMap(self):
        """
        Write chargeMap file, containing the charges for each fragment, to the chargeMap.dat file
        Format is: rows of 'filename' 'charge'

        Parameters:

        Class variables:
        """
        print('Charges have been mapped to fragments as follows:\n')
        f = open(self.path + '/chargeMap.dat', 'w')
        for i in range(len(self.files)):
            line = self.files[i].split('.pdb')[0], self.charges[i]  # Name charge
            print('{} -> {}'.format(*line))
            f.write('{} {}\n'.format(*line))
        f.close()
        print('\n')
        print('Map written to {}'.format(self.path + '/chargeMap.dat\n'))

    def getPDBasImage(self, file: str, type: str = 'pdb'):
        """
        Creates a 2D image of a pdb or xyz file using RDkit and returns a numpy array containing an image

        Parameters:
            :param string file: pdb or xyz file containing the fragment
            :param string type: file type, either 'pdb' or 'xyz'

        Class variables:

        Returns:
            - mol_img = image as a 3D numpy array (RGB)
        """
        if type == 'pdb':
            mol = Chem.MolFromPDBFile(file, removeHs=False)
        else:
            mol = Chem.MolFromXYZFile(file)
        AllChem.Compute2DCoords(mol)
        for atom in mol.GetAtoms():
            if atom.GetSymbol() != 'C':
                atom.SetProp("atomLabel", atom.GetSymbol())
        mol_img = Chem.Draw.MolToImage(mol, size=(600, 600))
        mol_img = np.asarray(mol_img)
        return mol_img

    def quit(self):
        """
        Close window and quit tkinter

        Parameters:

        Class variables:
        """
        self.window.quit()


class FragmentReviewGUI:
    def __init__(self, window, filepath: str, radius: float):
        """
        Class to present a GUI to the user for reviewing an extracted fragment before parametrization

        Parameters:
            :param window: TKinter window
            :param string filepath: path to the XYZ file of the fragment
            :param float radius: radius used for extraction (for display purposes)

        Class variables:
            - self.filepath = path to the XYZ fragment file
            - self.radius = radius used for extraction
            - self.window = TKinter window
            - self.confirmed = whether the user confirmed the fragment
        """
        self.filepath = filepath
        self.radius = radius
        self.window = window
        self.confirmed = False

        # Title label
        self.title_label = Label(window, text=f"Fragment extracted with radius {radius} Å",
                                  font=('Helvetica', 12, 'bold'))
        self.title_label.pack(pady=10)

        # Info label
        self.info_label = Label(window, text="Review the extracted substructure below.\n"
                                              "Capping hydrogens have been added at cut points.",
                                 font=('Helvetica', 10))
        self.info_label.pack(pady=5)

        # Plot the structure
        self.plot(filepath)

        # Button frame
        self.button_frame = Frame(window)
        self.button_frame.pack(pady=10)

        # Confirm button
        self.confirm_button = Button(self.button_frame, text="Confirm & Continue",
                                      command=self.confirm, bg='green', fg='white',
                                      font=('Helvetica', 11, 'bold'), padx=20, pady=5)
        self.confirm_button.pack(side=LEFT, padx=10)

        # Cancel button
        self.cancel_button = Button(self.button_frame, text="Cancel",
                                     command=self.cancel, bg='red', fg='white',
                                     font=('Helvetica', 11, 'bold'), padx=20, pady=5)
        self.cancel_button.pack(side=LEFT, padx=10)

    def plot(self, file: str):
        """
        Create and display an image of the XYZ file in the TKinter window

        Parameters:
            :param string file: xyz file to be displayed

        Class variables:
            - self.canvas = tkinter canvas that displays the image
        """
        fig = plt.figure(figsize=(8, 8))
        ax = plt.gca()
        img = self.getXYZasImage(file)
        if img is not None:
            ax.imshow(img)
        else:
            ax.text(0.5, 0.5, 'Could not render structure',
                    ha='center', va='center', fontsize=14)
        ax.axis('off')
        ax.set_title(f'Substructure Preview', fontsize=12)
        self.canvas = FigureCanvasTkAgg(fig, master=self.window)
        self.canvas.get_tk_widget().pack()
        self.canvas.draw()

    def getXYZasImage(self, file: str):
        """
        Creates a 2D image of an xyz file using RDkit and returns a numpy array containing an image.
        Bonds are determined using the XYZ class with the database txt files.

        Parameters:
            :param string file: xyz file containing the fragment

        Returns:
            - mol_img = image as a 3D numpy array (RGB), or None if failed
        """
        try:
            # Use the XYZ class to detect bonds (same as rest of codebase)
            db_path = os.path.split(os.path.split(__file__)[0])[0]
            xyz = XYZ(db_file=db_path + '/db/atom-radius.txt',
                      db_metal_file=db_path + '/db/metal-radius.txt')
            xyz.readXYZ(file)
            xyz.calculateDistanceMatrix()
            xyz.generateAdjacencyMatrix()
            xyz.generateLinkList()

            # Load molecule with RDKit
            mol = Chem.MolFromXYZFile(file)
            if mol is None:
                return None

            # Add bonds using the connectivity from XYZ class
            rwmol = Chem.RWMol(mol)
            bonded = []

            # Add organic bonds from linkList
            for i in range(len(xyz.linkList)):
                if len(xyz.linkList[i]) > 0:
                    for j in xyz.linkList[i]:
                        if [j, i] in bonded or [i, j] in bonded:
                            continue
                        try:
                            rwmol.AddBond(i, j, Chem.BondType.SINGLE)
                            bonded.append([i, j])
                        except Exception as e:
                            pass  # Bond might already exist

            # Add metal bonds (format: "metal_idx @ElementLigand_idx ligand_idx")
            for bond_info in xyz.metalBonds:
                parts = bond_info.split()
                if len(parts) >= 3:  # Fixed: was >= 4, but format has 3 parts
                    metal_idx = int(parts[0])
                    ligand_idx = int(parts[-1])
                    if [metal_idx, ligand_idx] not in bonded and [ligand_idx, metal_idx] not in bonded:
                        try:
                            rwmol.AddBond(metal_idx, ligand_idx, Chem.BondType.SINGLE)
                            bonded.append([metal_idx, ligand_idx])
                        except Exception as e:
                            pass

            mol = rwmol.GetMol()
            AllChem.Compute2DCoords(mol)

            for atom in mol.GetAtoms():
                if atom.GetSymbol() != 'C':
                    atom.SetProp("atomLabel", atom.GetSymbol())

            mol_img = Chem.Draw.MolToImage(mol, size=(700, 700))
            mol_img = np.asarray(mol_img)
            return mol_img
        except Exception as e:
            print(f'Error rendering structure: {e}')
            return None

    def confirm(self):
        """
        User confirmed the fragment, set flag and close window
        """
        self.confirmed = True
        self.window.quit()

    def cancel(self):
        """
        User cancelled, set flag and close window
        """
        self.confirmed = False
        self.window.quit()

    def is_confirmed(self) -> bool:
        """
        Return whether the user confirmed the fragment

        Returns:
            - bool: True if confirmed, False if cancelled
        """
        return self.confirmed
