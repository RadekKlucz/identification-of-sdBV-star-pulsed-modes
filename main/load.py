"""
load.py: The module with permanent important data and functions used in the main notebook.
"""
 
__name__ = "load"
__author__ = "Radoslaw Kluczewski"
__license__ = "MIT"
__version__ = "1.7"
__status__ = "production" 

from scipy.signal import find_peaks
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
import sys
import function as fc
import fourier as ft
# from comands import Calculate_the_noise_level  # uncomment this line if you want to run fortran files


class LoadAndFind:
    """
    The class contains variables and functions for loading data and finding peaks.
    """

    def __init__(self):
        """
        Data loading function and conversion to the required data in the notebook.        
        
        :param ft_data: list of  imported data
        :param frequencies: list of frequencies (axis x in data)
        :param ppt: list of ppt (periods per time, axis y in data)
        :param noise: float value of noise
        :param peaks: list of peaks 
        :param height: list of peaks heights 
        :param peak_position: list of peaks position
        :param periods: cycles per second converted to hertz list
        :param peaks_periods: list of peak conversions into periods
        :param export_data_periods: export the list to a file in the output folder
        """
        
        self.ft_data = pd.read_csv('./data/sdB93.dane', header=None, sep='\s+')
        self.data_column_1 = self.ft_data[0]
        self.data_column_2 = self.ft_data[1]
        self.frequencies_w = ft.fourier_transformation_without_fortran(self.data_column_1, self.data_column_2)[0]
        self.frequencies = np.array(list(filter(lambda x: x != 0, self.frequencies_w)))
        self.ppt = ft.fourier_transformation_without_fortran(self.data_column_1, self.data_column_2)[1][1:]
        self.noise = np.average(self.ppt) * 4
        #self.ft_data = pd.read_csv('./output/data/ft50.trf', sep='\s+', header=None)
        #self.frequencies = list(self.ft_data[0]) 
        #self.ppt = list(self.ft_data[1])  
        #self.noise = Calculate_the_noise_level() * 4  # uncoment this lines if you want evaluate fortran files
        self.peaks = find_peaks(self.ppt, height=self.noise)  
        self.height = self.peaks[1]['peak_heights']
        self.peak_position = self.frequencies[self.peaks[0]]
        self.periods = list(map(lambda x: 86400 / x, list(self.frequencies)))
        self.peaks_periods = list(map(lambda x: 86400 / x, list(self.peak_position)))
        self.export_data_periods = pd.DataFrame({'periods': self.periods, 'ppt': self.ppt})


    def export(self):
        """
        The function of exporting data in periods to the output folder.

        :returns: data file for exported periods 
        """

        return self.export_data_periods.to_csv("./output/fourier_data_periods.trf", sep='\t', header=None)


class SplitData(LoadAndFind):
    """
    The class contains a function to split data into small coubicks.
    This class has the potential to develop in the future.
    """

    def __init__(self):
        super().__init__()
        """
        The function contains valiables used to split data. 

        :param part_frequencies: split data in frequencies
        :param parts_ppt: split data in ppt
        :param sigmas: averages signals in small intervals
        :param thresholds: sigma 4 to intervals
        :param lines: list of lines to graph
        """
        self.parts_frequencies = np.array_split(self.frequencies, 25)
        self.parts_ppt = np.array_split(self.ppt, 25) 
        self.sigmas = fc.Function.average_list(self.parts_ppt)
        self.thresholds = np.array(self.sigmas) * 4
        self.lines = list()
        for i in range(25):
            self.lines.append([(i * 2, self.thresholds[i]), (i * 2 + 2, self.thresholds[i])])


class DistanceAndHistograms(LoadAndFind):
    """
    The class contains functions for drawing histograms.
    """

    def __init__(self):
        super().__init__()


    def calculate_the_distance_between_peaks(self):
        """
        The function with the calculated distance between the peaks in the graph. 

        :returns: list of distances between the peaks
        """
        self.distance_all = list()
        for i in self.peaks_periods:
            for j in self.peaks_periods[1:]:
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

        self.y, self.x, _ = plt.hist(data, bins, density=1, alpha=0.5, color='green')  
        self.start_trim_1 = math.floor(len(self.y) * start_interval)
        self.stop_trim_1 = math.floor(len(self.y) * end_interval)
        self.y = self.y[self.start_trim_1:self.stop_trim_1]
        self.x = self.x[self.start_trim_1:self.stop_trim_1]
        self.H, self.A, self.x0, self.sigma = fc.Function.gauss_fit(self.x, self.y)  
        self.x0 = self.x0 + (self.x[1] - self.x[0]) / 2
        self.xfit = np.arange(start_point_fit, end_point_fit)
        self.yfit = list()
        for i in self.xfit:
            self.yfit.append(fc.Function.gauss(i, self.H, self.A, self.x0, self.sigma))
        return [self.xfit, self.yfit, self.x0, self.sigma]


class MachingMods(LoadAndFind):
    """
    The class contains a function for mod associating.
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
        The function of assigning modes and numbers of detected peaks. 
        For the function to work properly, you must use 2 parameters from the three histograms.

        :param mean: floating point mean gauss matching parameter
        :param sigma: flaoting variancy gauss matching parameter
        :param mean_1: floating point mean gauss matching parameter
        :param sigma_1: flaoting variancy gauss matching parameter
        :param mean_2: floating point mean gauss matching parameter
        :param sigma_2: flaoting variancy gauss matching parameter
        :param list: list of data to assigning modes
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
        self.List_filtred = pd.DataFrame((list(self.peak_position), self.peaks_periods,
                                    self.height, self.L_List)).T  
        self.List_filtred[3].replace('', np.nan, inplace=True)
        self.List_filtred.dropna(subset=[3], inplace=True)
        self.List_reset = self.List_filtred.reset_index(drop=True)
        # self.List_filtred_latex = self.List_reset.to_latex(index=True)  # uncomment if you want to export to latex
        return self.List_reset


class Confirm(LoadAndFind):
    """
    The class has a function for confirming mods in the chart.
    """
    def __init__(self):
        super().__init__()


    def run_fourier_transformation_to_confirm(self):
        """
        The function that performs a Fourier transform to confirm the result with the input data in periods.

        :results: list of confirming floating data to be plotted on the graph
        """

        self.exported = self.export_data_periods.loc[self.export_data_periods['periods'] <= 20000]
        self.export_1 = self.exported.sort_values('periods')
        self.data_confirm = ft.fourier_transformation_without_fortran(list(self.export_1['periods']), self.ppt[:len(self.export_1)], 1)
        self.export_1.to_csv('./output/data_confirm.trf', sep='\t', header=None, index=None)
        return self.data_confirm


def mods(mean, sigma, mean_1, sigma_1, mean_2, sigma_2, list):
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
    :results: l number list and l number list in latex format
    """

    variables = MachingMods(mean, sigma, mean_1, sigma_1, mean_2, sigma_2, list)
    return variables.assing_mods()