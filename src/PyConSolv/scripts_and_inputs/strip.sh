#!/usr/bin/env bash

echo "Preparing to analyse trajectory"
echo "please enter the id of the atoms that should be used for alignment in the format: 1,2,4-10"

read atoms

echo "Aligning simulation on $atoms"

sed -i "s/replaceme/$atoms/g" align.in
sed -i "s/replaceme/$atoms/g" align_dry.in

cpptraj < align.in
cpptraj < dry_sim.in
cpptraj < align_dry.in
cpptraj < cluster_kmeans.in

echo "done"