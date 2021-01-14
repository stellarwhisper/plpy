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
        out_list = list(self.layer_names)
        for i in range(self.layers):
            out_list[i] += '_%d' % (self.thicknesses[i])
        return '/ '.join(out_list)

    def __str__(self):
        #string = ''
        out_list = []
        for i in range(self.layers):
            if self.n.y[i, 0].imag > 0.5:
                type_refractive_index = 'opaque'
            else:
                type_refractive_index = 'limpid'
            n_idx = self.n.y[i, 0]

            #string += "%dth: \t%s \t(%d nm) \tn = %.2f%s%.2fi, ... (%s)%s" % \
            #          (i + 1, self.layer_names[i], self.thicknesses[i],
            #           n_idx.real, '+-'[int(n_idx.imag < 0)], n_idx.imag,
            #           type_refractive_index,  #type(self.n(self.n.x[0])[idx]).__name__,
            #           '\n' * (i != self.layers - 1))     # return condition
            out_list.append(''.join([str(i+1),'th \t', self.layer_names[i], ' \t(',
                                  str(self.thicknesses[i]), ' nm) \tn = ',
                                  '%.2f %s %.2f' % (n_idx.real, '+-'[int(n_idx.imag < 0)], n_idx.imag),
                                  ' \t', type_refractive_index]))
        return '\n'.join(out_list)

    def __getitem__(self, val):
        #print(val)
        if isinstance(val, slice):
            self.layer_names = self.layer_names[val]
            self.thicknesses = self.thicknesses[val]
            self.n = interp1d(self.n.x, self.n.y[val])
            self.layers = np.shape(self.layer_names)[0]
            return self

        elif isinstance(val, int):
            if val < self.layers:
                indexed_layer = Multilayer([self.layer_names[val]],
                                           [self.thicknesses[val]],
                                           # if the object could receive interp1d
                                           # scipy object, replace here
                                           [self.n.y[val]], self.n.x)
                return indexed_layer
            else:
                raise IndexError("index out of range")
        else:
            raise TypeError("Invalid index type (%s)" % (type(val).__name__))

    def __add__(self, other):
        if isinstance(other, Multilayer):
            #print('add "%s" and "%s"' % (self.__repr__(), other.__repr__()))
            layer_names = np.append(self.layer_names, other.layer_names)
            thicknesses = np.append(self.thicknesses, other.thicknesses)
            #TODO: interpolation is need to the addition of n and wavelength
            n = np.append(self.n.y, other.n.y, axis=0)
            vac_wavelength = self.n.x  # it is needed to correct x
            result = Multilayer(layer_names, thicknesses, n, vac_wavelength)
            return result
        else:
            raise TypeError('can only concatenate Multilayer (not "%s") to Multilayer'
                            % (type(other).__name__))

    def __mul__(self, other):
        if isinstance(other, int):
            #print('multiplies "%s" by %d' % (self.__repr__(), other))
            result = self[::1]
            for i in range(other-1):
                result = result + self
            return result
        else:
            raise TypeError("Invalid multiplization type (%s)" % (type(val).__name__))

    def __sub__(self, other):
        if isinstance(other, Multilayer):
            rep_self, rep_other = self.__repr__(), other.__repr__()
            print('substract "%s" from "%s"' % (rep_other, rep_self))
            index = rep_self.find(rep_other)
            if index == -1:
                raise LookupError('The structure "%s" does not include "%s"' % (rep_self, rep_other))
            else:

                layer_names = np.append(self.layer_names[:index], self.layer_names[index+other.layers:])
                thicknesses = np.append(self.thicknesses[:index], self.thicknesses[index+other.layers:])
                # TODO: interpolation is need to the addition of n and wavelength
                n = np.append(self.n.y[:index], self.n.y[index+other.layers:], axis=0)
                vac_wavelength = self.n.x  # it is needed to correct x
                result = Multilayer(layer_names, thicknesses, n, vac_wavelength)
                return result

    def __floordiv__(self, other):
        """
        TODO: It is needed to be check
        :param other:
        :return:
        """
        if isinstance(other, int):
            mod = self.layers % other
            if mod == 0:
                div = self.layers // other
                div_result = []
                for i in range(div):
                    div_result.append(self[:i*mod])
                decision = True
                for i in range(div-1):
                    decision = decision and (div_result[i] == div_result[i+1])
                if decision:
                    return div_result[0]
                else:
                    raise ValueError("The Multilayer is not iterative.")
            else:
                raise  ValueError("The Multilayer could not be divide by %d" % (other))
        else:
            raise TypeError("Invalid floor division type (%s)" % (type(val).__name__))


    def __reversed__(self):
        self.layer_names = self.layer_names[::-1]
        self.thicknesses = self.thicknesses[::-1]
        self.n = interp1d(self.n.x, self.n.y[::-1])
        return self

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
               [[1.5+0.001j, 1.5+0.0001j, 1.6+0.0001j], [1.3+0.05j, 1.4+0.05j, 1.5+0.05j]],
               [300, 500, 700], active_layer = 2 )

b = np.complex128(1j)

a*3