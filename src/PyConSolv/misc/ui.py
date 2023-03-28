import matplotlib.pyplot as plt
import numpy as np
from tkinter import *
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw


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

    def getPDBasImage(self, file: str):
        """
        Creates a 2D image of a pdb file using RDkit and returns a numpy array containing an image

        Parameters:
            :param string file: pdb file containing the fragment

        Class variables:

        Returns:
            - mol_img = image as a 3D numpy array (RGB)
        """
        pdb = Chem.MolFromPDBFile(file, removeHs=False)
        AllChem.Compute2DCoords(pdb)
        for atom in pdb.GetAtoms():
            if atom.GetSymbol() != 'C':
                atom.SetProp("atomLabel", atom.GetSymbol())
        mol_img = Chem.Draw.MolToImage(pdb, size=(600, 600))
        mol_img = np.asarray(mol_img)
        return mol_img

    def quit(self):
        """
        Close window and quit tkinter

        Parameters:

        Class variables:
        """
        self.window.quit()
