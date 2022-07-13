from scipy.signal import find_peaks
# from comands import Calculate_the_noise_level
import pandas as pd
import numpy as np
import function as fc
import math
import matplotlib.pyplot as plt
import ft 


class LoadAndFind:
    """
    """

    def __init__(self):
        # self.ft_data = pd.read_csv(
         #   './output/ft50.trf', sep='\s+', header=None)  # data ft50.trf

        # self.frequencies = self.ft_data[0]  # list axis x in data
        # self.ppt = list(self.ft_data[1])  # list axis y in data
        self.ft_data = pd.read_csv('./data/sdB93.dane', header=None, sep='\s+')
        self.data_column_1 = self.ft_data[0]
        self.data_column_2 = self.ft_data[1]

        self.frequencies = ft.fourier_transformation_without_fortran(self.data_column_1, self.data_column_2)[0]
        self.ppt = ft.fourier_transformation_without_fortran(self.data_column_1, self.data_column_2)[1] 

        self.noise = np.average(self.ppt) * 4
        self.peaks = find_peaks(
            self.ppt, height=self.noise)  # finding peaks or put Calculate_the_noise_level() in height variable
        # list of the heights of the peaks
        self.height = self.peaks[1]['peak_heights']
        # list of the peaks position
        self.peak_position = self.frequencies[self.peaks[0]]

        # cycles per second to Hertz conversion
        self.periods = list(map(lambda x: 86400 / x, list(self.frequencies)))
        # peaks to periods conversion
        self.peaks_periods = list(
            map(lambda x: 86400 / x, list(self.peak_position)))

        self.export_data_periods = pd.DataFrame(
            {'periods': self.periods, 'ppt': self.ppt})  # export data in periods

    def export(self):
        return self.export_data_periods.to_csv("./output/ft50_periods.trf", sep='\t', header=None)


class SplitData(LoadAndFind):
    """
    """

    def __init__(self):
        super().__init__()
        self.parts_frequencies = np.array_split(
            self.frequencies, 25)  # split data in frequencies
        self.parts_ppt = np.array_split(self.ppt, 25)  # split data in ppt
        # averages signals in small intervals
        self.sigmas = fc.Function.average_list(self.parts_ppt)
        self.thresholds = np.array(self.sigmas) * 4  # sigma 4 to intervals

        self.lines = list()
        for i in range(25):
            self.lines.append([(i * 2, self.thresholds[i]),
                              (i * 2 + 2, self.thresholds[i])])


class DistanceAndHistograms(LoadAndFind):
    """
    """

    def __init__(self):
        super().__init__()


    def calculate_the_distance_between_peaks(self):
        self.distance_all = list()
        for i in self.peaks_periods:
            for j in self.peaks_periods[1:]:
                self.distance_all.append(math.dist([i], [j]))
        return self.distance_all


    def limiting(self, value_1, value_2, distance=None):
        """
        """

        if distance is None:
            self.distance = self.calculate_the_distance_between_peaks()
        else:
            self.distance = distance

        self.limiting_full = list()
        for i in self.distance:
            if i >= value_1 and i <= value_2:
                self.limiting_full.append(i)
        return self.limiting_full


    def fitting_gauss(self, data, start_interval, end_interval, start_point_fit, end_point_fit, bins=50):
        """
        """

        self.y, self.x, _ = plt.hist(data, bins, density=1,
                                     alpha=0.5, color='green')  # trimming the axis of the histogram
        self.start_trim_1 = math.floor(len(self.y) * start_interval)
        self.stop_trim_1 = math.floor(len(self.y) * end_interval)
        self.y = self.y[self.start_trim_1:self.stop_trim_1]
        self.x = self.x[self.start_trim_1:self.stop_trim_1]

        self.H, self.A, self.x0, self.sigma = fc.Function.gauss_fit(
            self.x, self.y)  # gauss parameters
        self.x0 = self.x0 + (self.x[1] - self.x[0]) / 2

        self.xfit = np.arange(start_point_fit, end_point_fit)  # fitting gauss
        self.yfit = list()
        for i in self.xfit:
            self.yfit.append(fc.Function.gauss(i, self.H, self.A, self.x0, self.sigma))
        return [self.xfit, self.yfit, self.x0, self.sigma]


class MachingMods(LoadAndFind):
    """
    """

    def __init__(self, mean, sigma, mean_1, sigma_1, mean_2, sigma_2, list):
        super().__init__()
        self.mean = mean
        self.sigma = sigma
        self.mean_1 = mean_1
        self.sigma_1 = sigma_1
        self.mean_2 = mean_2
        self.sigma_2 = sigma_2
        self.list = list

    def assing_mods(self):
        """
        """

        self.L_List = list()
        for i in self.list:
            List_l = str()
            if fc.Function.is_mode(i, self.mean, self.sigma, self.list):
                List_l += '1 '
            if fc.Function.is_mode(i, self.mean_1, self.sigma_1, self.list):
                List_l += '2 '
            if fc.Function.is_mode(i, self.mean_2, self.sigma_2, self.list):
                List_l += '3 '
            self.L_List.append(List_l)
        self.List_filtred = pd.DataFrame((list(self.peak_position), self.peaks_periods,
                                    self.height, self.L_List)).T  # filttring empty rows
        self.List_filtred[3].replace('', np.nan, inplace=True)
        self.List_filtred.dropna(subset=[3], inplace=True)
        self.List_reset = self.List_filtred.reset_index(drop=True)
        self.List_filtred_latex = self.List_reset.to_latex(index=True)
        return self.List_reset, self.List_filtred_latex


class Confirm(LoadAndFind):
    """
    """
    def __init__(self):
        super().__init__()


    def run_fourier_transformation(self):
        """
        """

        self.exported = self.export_data_periods.loc[self.export_data_periods['periods'] <= 20000]
        self.export_1 = self.exported.sort_values('periods')
        self.data_confirm = ft.fourier_transformation_without_fortran(list(self.export_1['periods']), 
                                                                      self.ppt[:len(self.export_1)], 
                                                                      normalize_number=1)
        self.export_1.to_csv('./output/data_sorted.trf', sep='\t',
                             header=None, index=None)
        return self.data_confirm


def Mods(mean, sigma, mean_1, sigma_1, mean_2, sigma_2, list):
    """
    """

    variables = MachingMods(mean, sigma, mean_1,
                            sigma_1, mean_2, sigma_2, list)
    return variables.assing_mods()
