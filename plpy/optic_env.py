import pandas as pd
import numpy as np
from scipy.interpolate import interp1d

class Multilayer(object):
    """
    MultilayerArch has a set of multilayer information
    for transfer matrix calculation in a given condition.

    """

    def __init__(self, layer_names=[], thicknesses=[], n=[[]], vac_wavelenth=[]):
        """
        The object Multilayer receives only necessory components for
        construct optical multilayer configuration to calculate
        """

        if np.shape(self.layer_names) == np.shape(self.thicknesses):
            if np.shape(self.thicknesses) == (np.shape(self.n.y)[0],):
                self.layers = np.shape(self.layer_names)[0]
            else:
                raise ValueError("thicknesses and refractive indecies array must be equal "
                                 "in length along 0 axis %s %s" %
                                 (np.shape(self.thicknesses), np.shape(self.n.y)))
        else:
            raise ValueError("layer names and thicknesses array must be equal in length %s %s" %
                             (np.shape(self.layer_names), np.shape(self.thicknesses)))

        self.layer_names = np.array(layer_names)
        self.thicknesses = np.array(thicknesses)
        self.n = interp1d(vac_wavelenth, n)



    def __repr__(self):
        string = ''
        for idx in range(self.layers):
            string += str(idx+1) + 'th layer thickness of ' + self.layer_names[idx] + \
                ': ' + str(self.thicknesses[idx]) + ' with complex n: ' + \
                str(isinstance(self.n(self.n.x[0])[idx], complex)) + \
                ('\n' * (idx != self.layers - 1)) # add return
        return string

    def __getitem__(self, index=0):
        if index < self.layers:
            indexed_layer = Multilayer([self.layer_names[index]], [self.thicknesses[index]],
                                       [self.n.y[index]], self.n.x)
            return indexed_layer
        else:
            raise IndexError

    def __add__(self, other):
        if isinstance(other, Multilayer):
            layer_names = np.append(self.layer_names, other.layer_names)
            thicknesses = np.append(self.thicknesses, other.thicknesses)
            n = np.append(self.n.y, other.n.y, axis=0)
            vac_wavelength = self.n.x  # it is needed to correct x
            self = Multilayer(layer_names, thicknesses, n, vac_wavelength)
            return self
        else:
            raise TypeError('can only concatenate Multilayer (not "%s") to Multilayer' % (str(type(other).__name__)))

    def pop(self, index):
        pop_item = self[index]
        self.layers -= 1
        self.layer_names = np.delete(self.layer_names, index)
        self.thicknesses = np.delete(self.thicknesses, index)
        self.n = interp1d(self.n.x, np.delete(self.n.y, index, axis=0))
        return pop_item

a = Multilayer(['ITO', 'PEDOT'], [200, 50],
               [[1.5+0.01j, 1.5+0.001j, 1.6+0.001j], [1.3+0.05j, 1.4+0.05j, 1.5+0.05j]], [300, 500, 700])