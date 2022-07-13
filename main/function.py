import numpy as np
from scipy.optimize import curve_fit


class Function:
    '''Clasa funkcji opis co zawiera. Dodac odnosniki do teorii, skad wzialem to'''

    def __init__(self, x, H, A, mu, sigma):  # to nie działa !!!
        self.x = x
        self.H = H
        self.A = A
        self.mu = mu
        self.sigma = sigma

    def gauss(x, H, A, mu, sigma):
        '''funkcja gaussa'''
        return H + A * np.exp(- (x - mu) ** 2 / (2 * sigma ** 2))

    def gauss_fit(x_1, y_1):
        '''funkcja dopasowujaca gaussa'''
        mean_1 = sum(x_1 * y_1) / sum(y_1)
        sigma_1 = np.sqrt(sum(y_1 * (x_1 - mean_1) ** 2)) / sum(y_1)
        popt, pcov = curve_fit(Function.gauss, x_1, y_1, p0=[
                               min(y_1), max(y_1), mean_1, sigma_1])
        return popt

    def is_mode(j, mean_2, sigma_2, list):
        '''funkcja do identyfikacji modow, liczb l = 1, 2, 3'''
        for i in list:
            if True:
                sigma_2 = abs(sigma_2) * 1
                if j - mean_2 - sigma_2 < i < j - mean_2 + sigma_2 or j + mean_2 - sigma_2 < i < j + mean_2 + sigma_2:
                    return True
        return False

    def split_data(point_1, point_2):
        '''do dzielenia danych'''
        value_1, value_2 = divmod(len(point_1), point_2)
        return (point_1[i * value_1 + min(i, value_2): (i + 1) * value_1 + min(i + 1, value_2)] for i in range(point_2))

    def average_list(list):
        '''Function with calculate average level form data list'''
        averages = []
        for i in list:
            averages.append(np.average(i))
        return averages


# def Gauss(x, H, A, mu, sigma):
 #   variables = Function(x, H, A, mu, sigma)
  #  return variables.gauss(x, H, A, mu, sigma)
