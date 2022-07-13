"""
function.py: This module contains all the necessary functions used in the project.
"""
 
__name__ = "function"
__author__ = "Radoslaw Kluczewski"
__license__ = "MIT"
__version__ = "1.3"
__status__ = "production" 

import numpy as np
from scipy.optimize import curve_fit


class Function:
    '''
    This class contains the basic math functions and functions used in both the load module and the master script.
    More on the Gaussian function can be found here: <https://en.wikipedia.org/wiki/Gaussian_function>    
    '''

    def gauss(x, H, A, mu, sigma):
        '''
        This function represent the mathematical form of the Gaussian function.

        :param x: float point of values
        :param H: float parameter
        :param A: float parameter
        :param mu: float parameter
        :param sigma: float parameter
        :results: Mathematically express a Gaussian function to fit on a graph.
        '''
        return H + A * np.exp(- (x - mu) ** 2 / (2 * sigma ** 2))


    def gauss_fit(x_1, y_1):
        '''
        This function fit the Gaussian function on a graph and return best fit values for parameters fo the gicen model.
        
        :param x_1: float value of a point on the axis x
        :param y_1: float value of a point on the axis y
        :param mean_1: float expected value
        :param sigma_1: float variance
        :returns: popt list of the best fit values for parameters
        '''
        mean_1 = sum(x_1 * y_1) / sum(y_1)
        sigma_1 = np.sqrt(sum(y_1 * (x_1 - mean_1) ** 2)) / sum(y_1)
        popt, _ = curve_fit(Function.gauss, x_1, y_1, 
                               p0=[min(y_1), max(y_1), mean_1, sigma_1])
        return popt


    def is_mode(j, mean_2, sigma_2, list):
        '''
        This function identifies l numbers such as 1, 2, 3.

        :param j: float value from the given list
        :param mean_2: float parameter calculated form the gauss fit 
        :param sigma_2: float parameter calculated from the gauss fit
        :param list: list of values 
        :results: List with identified numbers l
        '''
        for i in list:
            if True:
                sigma_2 = abs(sigma_2) * 1
                if j - mean_2 - sigma_2 < i < j - mean_2 + sigma_2 or j + mean_2 - sigma_2 < i < j + mean_2 + sigma_2:
                    return True
        return False


    def split_data(point_1, point_2):
        '''
        Function to split data on the small list.
        
        :param point_1: float first point
        :param point_2: float last point
        :returns: tuple of split data  
        '''
        value_1, value_2 = divmod(len(point_1), point_2)
        return (point_1[i * value_1 + min(i, value_2): (i + 1) * value_1 + min(i + 1, value_2)] for i in range(point_2))


    def average_list(list):
        '''
        Function with calculate average level form data list.

        :param list: list of float values
        :returns: averages list of values
        '''
        averages = []
        for i in list:
            averages.append(np.average(i))
        return averages