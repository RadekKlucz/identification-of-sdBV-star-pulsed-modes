from astropy.timeseries import LombScargle
import numpy as np


class FourierTransform:
    """
    """

    def __init__(self, time, amplitude, normilize_number):
        self.time = time
        self.amplitude = amplitude
        self.normilize_number = normilize_number


    def fourier_transform(self):
        """This founction 
        """

        self.frequency, self.power = LombScargle(self.time, self.amplitude, normalization='standard').autopower(
                                                 minimum_frequency=0, maximum_frequency=50)
        self.amplitude = np.sqrt(abs(self.power)) * self.normilize_number
        self.norm = np.linalg.norm(self.amplitude)
        self.normal_array = self.amplitude / self.norm
        return self.frequency, self.amplitude
    

def fourier_transformation_without_fortran(time, amplitude, normalize_number=6.1):
    output = FourierTransform(time, amplitude, normilize_number=6.1)
    return output.fourier_transform()

