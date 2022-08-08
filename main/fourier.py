"""
fourier.py: This module includes the Fourier transform which is the main function in this project.
"""
 
__name__ = "fourier"
__author__ = "Radoslaw Kluczewski"
__license__ = "MIT"
__version__ = "1.1"
__status__ = "accomplished" 

from astropy.timeseries import LombScargle
import numpy as np


class FourierTransform:
    """
    Class includes the Fourier transformation function.
    """

    def __init__(self, time, amplitude, normalization_point):
        self.time = time
        self.amplitude = amplitude
        self.normalization_point = normalization_point


    def fourier_transform(self):
        """
        This function is a Fourier transform for unequally spaced data.

        :param time: float of JDHel values
        :param amplitude: float of ppt (parts per thousand), which is the equivalent of the relative flux intensity of a star
        :param normalization_point: float number to normalize to 1 data, the base value is 6.1
        :returns: floating point frequency and amplitude list
        """

        self.frequency, self.power = LombScargle(self.time, self.amplitude, normalization='standard').autopower(
                                                 minimum_frequency=0, maximum_frequency=50)
        self.amplitude = np.sqrt(abs(self.power)) * self.normalization_point
        return self.frequency, self.amplitude
    

def fourier_transformation_without_fortran(time, amplitude, normalization_point=6.1):
    """
    The function of calculating the Fourier transform from basic data.
    """
    
    output = FourierTransform(time, amplitude, normalization_point=6.1)
    return output.fourier_transform()