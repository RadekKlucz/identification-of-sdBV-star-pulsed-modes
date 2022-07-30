"""
load.py: The module with permanent important data and functions used in the main notebook.
"""

__name__ = "load"
__author__ = "Radoslaw Kluczewski"
__license__ = "MIT"
__version__ = "1.8"
__status__ = "production"

from scipy.signal import find_peaks
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
import function as fc
import fourier as ft
from comands import calculate_the_noise_level


class LoadAndFind:
    """
    The class contains variables and functions for loading data and finding peaks.
    """

    def __init__(self, data, method):
        self.data = data
        self.method = method.lower()

    def load_function(self):
        """
        Data loading function and conversion to the required data in the notebook and exporting data in periods to the output folder.        

        :input param data: string input your data file
        :input param method: string choose method (p >>python<< or f >>fortran<<)
        :returns: list of frequencies (axis x in data), list of ppt (periods per time, axis y in data), 
                  float value of noise, list of peaks, list of peaks heights, 
                  list of peaks position, cycles per second converted to hertz list, 
                  list of peak conversions into periods, data file for exported periods 
        """

        if self.method == 'p':
            self.ft_data = pd.read_csv(self.data, header=None, sep='\s+')
            self.data_column_1 = self.ft_data[0]
            self.data_column_2 = self.ft_data[1]
            self.frequencies_with_zero = ft.fourier_transformation_without_fortran(
                self.data_column_1, self.data_column_2)[0]
            self.frequencies = np.array(
                list(filter(lambda x: x != 0, self.frequencies_with_zero)))
            self.ppt = ft.fourier_transformation_without_fortran(
                self.data_column_1, self.data_column_2)[1][1:]
            self.noise = np.average(self.ppt) * 4
        else:
            self.ft_data = pd.read_csv(
                './output/ft50.trf', sep='\s+', header=None)
            self.frequencies = list(self.ft_data[0])
            self.ppt = list(self.ft_data[1])
            self.noise = calculate_the_noise_level() * 4

        self.peaks = find_peaks(self.ppt, height=self.noise)
        self.height = self.peaks[1]['peak_heights']
        self.frequencies_array = np.array(self.frequencies)
        self.peak_position = self.frequencies_array[self.peaks[0]]
        self.periods = list(map(lambda x: 86400 / x, list(self.frequencies)))
        self.peaks_periods = list(
            map(lambda x: 86400 / x, list(self.peak_position)))
        self.export_data_periods = pd.DataFrame(
            {'periods': self.periods, 'ppt': self.ppt})
        self.export_data_periods.to_csv(
            "./output/data/fourier_data_periods.trf", sep='\t', header=None)

        return self.frequencies, self.ppt, self.noise, self.peaks, self.height, self.peak_position, self.periods, self.peaks_periods, \
            self.export_data_periods


class SplitData(LoadAndFind):
    """
    The class contains a function to split data into small coubicks.
    This class has the potential to develop in the future.
    """

    def __init__(self, data, method):
        super().__init__(data, method)
        """
        The method contains valiables used to split data. 

        :param part_frequencies: split data in frequencies
        :param parts_ppt: split data in ppt
        :param sigmas: averages signals in small intervals
        :param thresholds: sigma 4 to intervals
        :param lines: list of lines to graph
        """
        self.parts_frequencies = np.array_split(load_data(data, method)[0], 25)
        self.parts_ppt = np.array_split(load_data(data, method)[1], 25)
        self.sigmas = fc.Function.average_list(self.parts_ppt)
        self.thresholds = np.array(self.sigmas) * 4
        self.lines = list()
        for i in range(25):
            self.lines.append([(i * 2, self.thresholds[i]),
                              (i * 2 + 2, self.thresholds[i])])


class DistanceAndHistograms(LoadAndFind):
    """
    The class contains functions for drawing histograms.
    """

    def __init__(self, data, method):
        super().__init__(data, method)

    def calculate_the_distance_between_peaks(self):
        """
        The function with the calculated distance between the peaks in the graph. 

        :returns: list of distances between the peaks
        """
        self.distance_all = list()
        self.peaks_periods_1 = load_data(self.data, self.method)[7]
        for i in self.peaks_periods_1:
            for j in self.peaks_periods_1[1:]:
                self.distance_all.append(math.dist([i], [j]))
        return self.distance_all

    def limiting(self, value_1, value_2, distance=None):
        """
        The function to limiting data in the histograms.

        :param value_1: float lower limit of the data range
        :param value_2: float upper bound of the data range
        :param distance: float value to select the preferred distance between peaks (e.g. 150)
        returns: list of values from a given range
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
        Function for fitting Gaussian functions to histograms.

        :param data: input list data
        :param start_interval: float beginning of the data interval
        :param end_interval: float end of data range
        :param start_point_fit: float start point of function matching
        :param end_point_fit: float end point of function matching
        :param bins: float number of bins on the chart
        :results: list of floating point parameters to match the histogram
        """

        self.y, self.x, _ = plt.hist(
            data, bins, density=1, alpha=0.5, color='green')
        self.start_trim_1 = math.floor(len(self.y) * start_interval)
        self.stop_trim_1 = math.floor(len(self.y) * end_interval)
        self.y = self.y[self.start_trim_1:self.stop_trim_1]
        self.x = self.x[self.start_trim_1:self.stop_trim_1]
        self.H, self.A, self.x0, self.sigma = fc.Function.gauss_fit(
            self.x, self.y)
        self.x0 = self.x0 + (self.x[1] - self.x[0]) / 2
        self.xfit = np.arange(start_point_fit, end_point_fit)
        self.yfit = list()
        for i in self.xfit:
            self.yfit.append(fc.Function.gauss(
                i, self.H, self.A, self.x0, self.sigma))
        return [self.xfit, self.yfit, self.x0, self.sigma]


class MachingMods(LoadAndFind):
    """
    The class contains a function for mod associating.
    """

    def __init__(self, mean, sigma, mean_1, sigma_1, mean_2, sigma_2, list, data, method):
        super().__init__(data, method)
        self.mean = mean
        self.sigma = sigma
        self.mean_1 = mean_1
        self.sigma_1 = sigma_1
        self.mean_2 = mean_2
        self.sigma_2 = sigma_2
        self.list = list

    def assing_mods(self):
        """
        The function of assigning modes and numbers of detected peaks. 
        For the function to work properly, you must use 2 parameters from the three histograms.

        :param mean: floating point mean gauss matching parameter
        :param sigma: flaoting variancy gauss matching parameter
        :param mean_1: floating point mean gauss matching parameter
        :param sigma_1: flaoting variancy gauss matching parameter
        :param mean_2: floating point mean gauss matching parameter
        :param sigma_2: flaoting variancy gauss matching parameter
        :param list: list of data to assigning modes
        :input param data: string input your data file
        :input param method: string choose method (p >>python<< or f >>fortran<<)
        :results: l number list and l number list in latex format
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
        self.List_filtred = pd.DataFrame((list(load_data(self.data, self.method))[5], load_data(self.data, self.method)[7],
                                          load_data(self.data, self.method)[4], self.L_List)).T
        self.List_filtred[3].replace('', np.nan, inplace=True)
        self.List_filtred.dropna(subset=[3], inplace=True)
        self.List_reset = self.List_filtred.reset_index(drop=True)
        # self.List_filtred_latex = self.List_reset.to_latex(index=True)  # uncomment if you want to export to latex
        return self.List_reset


class Confirm(LoadAndFind):
    """
    The class has a function for confirming mods in the chart.
    """

    def __init__(self, data, method):
        super().__init__(data, method)

    def run_fourier_transformation_to_confirm(self):
        """
        The function that performs a Fourier transform to confirm the result with the input data in periods.

        :results: list of confirming floating data to be plotted on the graph
        """
        self.export_data_periods_1 = load_data(self.data, self.method)[8]
        self.ppt_1 = load_data(self.data, self.method)[1]
        self.exported = self.export_data_periods_1.loc[self.export_data_periods_1['periods'] <= 20000]
        self.export_1 = self.exported.sort_values('periods')
        self.data_confirm = ft.fourier_transformation_without_fortran(
            list(self.export_1['periods']), self.ppt_1[:len(self.export_1)], 1)
        self.export_1.to_csv('./output/data/data_confirm.trf',
                             sep='\t', header=None, index=None)
        return self.data_confirm


def mods(mean, sigma, mean_1, sigma_1, mean_2, sigma_2, list, data, method):
    """
    The function of assigning modes and numbers of detected peaks. 
    For the function to work properly, you must use 2 parameters from the three histograms.

    :param mean: floating point mean gauss matching parameter
    :param sigma: flaoting variancy gauss matching parameter
    :param mean_1: floating point mean gauss matching parameter
    :param sigma_1: flaoting variancy gauss matching parameter
    :param mean_2: floating point mean gauss matching parameter
    :param sigma_2: flaoting variancy gauss matching parameter
    :param list: list of data to assigning modes
    :input param data: string input your data file
    :input param method: string choose method (p >>python<< or f >>fortran<<)
    :results: l number list and l number list in latex format
    """

    variables = MachingMods(mean, sigma, mean_1, sigma_1,
                            mean_2, sigma_2, list, data, method)
    return variables.assing_mods()


def load_data(data, method):
    """
    The function of loading all the necessary data into the chart. Also returns the file with data in periods. 
    There are frequencies at index 0, index 1 is ppt, index 2 is noise, index 3 is peaks, index 4 is height, 
    index 5 is peak position, index 6 is periods, index 7 is peaks periods

    :input param data: string input your data file
    :input param method: string choose method (p >>python<< or f >>fortran<<)
    :returns: list of frequencies, ppt, noise, peaks, height, peak_position, periods, peaks_periods, file with data in periods
    """
    variables = LoadAndFind(data, method)
    return variables.load_function()
