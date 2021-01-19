import pandas as pd
import numpy as np
import copy
from scipy.interpolate import interp1d

class Multilayer(object):
    """
    MultilayerArch has a set of multilayer information
    for transfer matrix calculation in a given condition.
    """

    def __init__(self, layer_names=list(), thicknesses=list(), refractive_index=interp1d,
                 **kwargs):
        """
        The object Multilayer receives only necessory components for
        construct optical multilayer configuration to calculate
        :param layer_names:
        :param thicknesses:
        :param refractive_index:
        :param kwargs: vac_wavelength, active_layer
        """
        # TODO: parameter check, and allow to receive n and k seperately

        # declare Multilayer parameters
        # kwargs
        if 'active_layer' in kwargs:
            self.active_layer = kwargs['active_layer']
        if 'vac_wavelength' in kwargs:
            vac_wavelength = kwargs['vac_wavelength']

        # method parameters
        self.layer_names = np.array(layer_names,ndmin = 1)
        self.thicknesses = np.array(thicknesses,ndmin = 1)
        if isinstance(refractive_index, interp1d):
            self.n = refractive_index
        elif isinstance(refractive_index, (list, np.ndarray)) and not vac_wavelength is None:
            self.n = interp1d(vac_wavelength, np.array(refractive_index, ndmin=2))
        else:
            raise TypeError("")
        self.layers = len(self.layer_names)

        ## parameters check
        if self.layer_names.ndim != 1:
            raise ValueError("the dinemsion of layer_names must be 1")
        if self.thicknesses.ndim != 1:
            raise ValueError("the dinemsion of thicknesses must be 1")
        if self.n.y.ndim != 2:
            raise ValueError("the dinemsion of refractive indices must be less than 3")
        # length: layers = len(layer_names == thicknesses == n.y(# of row))
        # length: n.x == n.y(# of column) (maybe, included in interp1d)
        if self.thicknesses.shape[0] != self.layers:
            raise ValueError("the dimension of thicknesses (%s) must be equal to %s" %
                             (self.thicknesses.shape[0], self.layers))
        if self.n.y.shape[0] != self.layers:
            raise ValueError("the rows of refractive index (%s) must be equal to %s" %
                             (self.n.y.shape[0], self.layers))
        # data types
        """if self.layer_names.dtype.kind != "U":
            raise ValueError("layer name must be string, not %s" % (np.array(layer_names).dtype))"""

    def __repr__(self):
        """
        Privides representative structure of the object
        :return: string
        """
        out_list = ['%s_%d' % (self.layer_names[i], self.thicknesses[i])
                    for i in range(self.layers)]
        return '/ '.join(out_list)

    def __str__(self):
        """
        Provides detail structure with representative refractive index at i = 0
        :return: string in detail
        """
        out_list = ['%dth \t %s \t(%d nm) \t%.2f%s%.2f at %dnm (%s)' %
                    (i + 1, self.layer_names[i], self.thicknesses[i], self.n.y[i, 0].real,
                     '+-'[int(self.n.y[i, 0].imag < 0)], abs(self.n.y[i, 0].imag), self.n.x[0] ,
                     ['transparent', 'opaque'][int(self.n.y[i, 0].imag > 0.5)])
                    for i in range(self.layers)]
        return '\n'.join(out_list)

    def __getitem__(self, val):
        """
        Indexing object, could be sliced
        :param val:
        :return: index/sliced object
        """
        if not (isinstance(val, int) or isinstance(val, slice)):
            raise TypeError('parameter must be int or slice not %s' % (type(val)))
        if isinstance(val, int):
            if val >= self.layers:
                raise IndexError('index %d is out of bounds for the object %s layers' % (val, self.layers))

        indexed = Multilayer(self.layer_names[val], self.thicknesses[val], self.n.y[val],
                             vac_wavelength = self.n.x)
        return indexed

    def __add__(self, other):
        """
        Merge layer objects
        :param other:
        :return: merged layers
        """
        if not isinstance(other, Multilayer):
            raise TypeError('can only concatenate Multilayer (not "%s") to Multilayer'
                            % (type(other).__name__))

        # print('add "%s" and "%s"' % (self.__repr__(), other.__repr__()))
        layer_names = np.append(self.layer_names, other.layer_names)
        thicknesses = np.append(self.thicknesses, other.thicknesses)
        #TODO: interpolation is need to the addition of n and wavelength
        n = np.append(self.n.y, other.n.y, axis=0)
        vac_wavelength = self.n.x  # it is needed to correct x
        result = Multilayer(layer_names, thicknesses, interp1d(vac_wavelength, n))
        return result

    def __mul__(self, val):
        """

        :param val:
        :return:
        """
        if not isinstance(val, int):
            raise TypeError("Invalid multiplization type (%s)" % (type(val).__name__))

        # print('multiplies "%s" by %d' % (self.__repr__(), other))
        result = copy.copy(self)
        for i in range(val - 1):
            result = result + self
        return result


    def __sub__(self, other):
        """

        :param other:
        :return:
        """
        if not isinstance(other, Multilayer):
            raise TypeError('only Multilayer could be substratable')

        rep_self, rep_other = repr(self), repr(other)
        # print('substract "%s" from "%s"' % (rep_other, rep_self))
        index = rep_self.find(rep_other)
        if index == -1:
            raise LookupError('The structure "%s" does not include "%s"' % (rep_self, rep_other))
        result = copy.copy(self)
        len = other.layers
        for i in range(len): # pop layers from last index
            result.pop(index + len - 1 - i)

        return result

    def __floordiv__(self, other):
        """
        divide layers if the layers are iterative
        :param other: integer
        :return:
        """
        if not isinstance(other, int):
            raise TypeError("Invalid floor division type (%s)" % (type(other).__name__))
        mod = self.layers % other
        if mod != 0:
            raise ValueError("The Multilayer could not be divide by %d" % (other))

        div = self.layers // other
        div_result = []
        for i in range(other):
            div_result.append(self[i * div : (i+1) * div])
        same = True
        for i in range(other - 1):
            same = same and (div_result[i] == div_result[i+1])
        if same:
            return div_result[0]
        else:
            raise ValueError("The Multilayer is not iterative.")

    def __reversed__(self):
        """
        Reverse layers
        :return: reversed(self)
        """
        self.layer_names = self.layer_names[::-1]
        self.thicknesses = self.thicknesses[::-1]
        self.n = interp1d(self.n.x, self.n.y[::-1])
        return self

    def __eq__(self, other):
        if not isinstance(other, Multilayer):
            raise TypeError('could compare with %s' % (type(other)))

        return repr(self) == repr(other)

    def pop(self, index):
        """

        :param index:
        :return:
        """
        # error handling could be at indexing (__getitem__)
        # print("pop %dth layer" % (index+1))
        pop_item = self[index]
        self.layers -= 1
        self.layer_names = np.delete(self.layer_names, index)
        self.thicknesses = np.delete(self.thicknesses, index)
        self.n = interp1d(self.n.x, np.delete(self.n.y, index, axis=0))
        return pop_item

# test instance
a = Multilayer(['ITO', 'PEDOT'], [200, 50],
                  [[1.5+0.001j, 1.5+0.0001j, 1.6+0.0001j], [1.3+0.05j, 1.4+0.05j, 1.5+0.05j]],
                  vac_wavelength=[300, 500, 700], active_layer = 2)

b = a * 3