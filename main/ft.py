from astropy.timeseries import LombScargle
from load import LoadAndFind
import numpy as np


class FourierTransform(LoadAndFind):
    ''''''

    def __init__(self):
        super().__init__()

    def fourier_transform(self):
        ''''''
        self.frequency, self.power = LombScargle(self.frequencies, self.ppt, normalization='standard').autopower(
            minimum_frequency=0, maximum_frequency=50)
        # power to amplitude
        self.amplitude = np.sqrt(abs(self.power)) / max(self.power)
        self.norm = np.linalg.norm(self.amplitude)
        self.normal_array = self.amplitude / self.norm
        return self.frequency, self.normal_array, self.amplitude
    

def fourier_transformation_without_fortran():
    output = FourierTransform()
    return output.fourier_transform()
