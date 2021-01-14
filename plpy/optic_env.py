import pandas as pd
import numpy as np
from scipy.interpolate import interp1d

class Multilayer(object):
    """
    MultilayerArch has a set of multilayer information
    for transfer matrix calculation in a given condition.

    """

    def __init__(self, layer_names=[], thicknesses=[], n=[[]], vac_wavelenth=[],
                 **kwargs):
        """
        The object Multilayer receives only necessory components for
        construct optical multilayer configuration to calculate
        """
        #TODO: parameter check, and allow to receive interp1d n values and/or n and k seperately

        ## parameters check
        # for refractive index
        if len(np.shape(n)) == 1:
            n = np.reshape(n,(1,len(n)))

        # layer name data type
        if np.array(layer_names).dtype.kind == "U":
            pass
        else:
            raise ValueError("layer names should be string, not %s" % (np.array(layer_names).dtype))

        # numer of layers
        if np.shape(layer_names) == np.shape(thicknesses):
            pass
        else:
            raise ValueError("layer names and thicknesses array must be equal in length %s %s" %
                             (np.shape(layer_names), np.shape(thicknesses)))
        if np.shape(thicknesses)[0] == np.shape(n)[0]:
            self.layers = np.shape(layer_names)[0]
        else:
            raise ValueError("thicknesses and refractive indecies array must be equal "
                             "in length along 0 axis %s %s" %
                             (np.shape(thicknesses), np.shape(n)))

        # kwargs
        if 'active_layer' in kwargs:
            self.active_layer = kwargs['active_layer']

        # declare Multilayer parameters
        self.layer_names = np.array(layer_names)
        self.thicknesses = np.array(thicknesses)
        self.n = interp1d(vac_wavelenth, n)



    def __repr__(self):
        string = ''
        for idx in range(self.layers):
            string += "%dth: %s \t(%d) \t%s%s" \
                      % (idx + 1, self.layer_names[idx], self.thicknesses[idx],
                         type(self.n(self.n.x[0])[idx]).__name__,
                         ('\n' * (idx != self.layers - 1))) # return condition
        return string

    def __str__(self):
        # =
        #print(string)
        return self.__repr__()

    def __getitem__(self, index=0):
        if index < self.layers:
            indexed_layer = Multilayer([self.layer_names[index]], [self.thicknesses[index]],
                                       [self.n.y[index]], self.n.x)
            return indexed_layer
        else:
            raise IndexError

    def __add__(self, other):
        if isinstance(other, Multilayer):
            print("add %s and %s" % (self.__repr__(), other.__repr__()))
            layer_names = np.append(self.layer_names, other.layer_names)
            thicknesses = np.append(self.thicknesses, other.thicknesses)
            #TODO: interpolation is need to the addition of n and wavelength
            n = np.append(self.n.y, other.n.y, axis=0)
            vac_wavelength = self.n.x  # it is needed to correct x
            self = Multilayer(layer_names, thicknesses, n, vac_wavelength)
            return self
        else:
            raise TypeError('can only concatenate Multilayer (not "%s") to Multilayer'
                            % (type(other).__name__))

    def pop(self, index):
        print("pop %dth layer" % (index+1))
        pop_item = self[index]
        self.layers -= 1
        self.layer_names = np.delete(self.layer_names, index)
        self.thicknesses = np.delete(self.thicknesses, index)
        self.n = interp1d(self.n.x, np.delete(self.n.y, index, axis=0))
        return pop_item


# test instance

a = Multilayer(['ITO', 'PEDOT'], [200, 50],
               [[1.5+0.01j, 1.5+0.001j, 1.6+0.001j], [1.3+0.05j, 1.4+0.05j, 1.5+0.05j]],
               [300, 500, 700], active_layer = 2 )

b = np.complex128(1j)