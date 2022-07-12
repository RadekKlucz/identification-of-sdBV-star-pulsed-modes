#from concurrent.futures.process import _check_system_limits
from sqlite3 import ProgrammingError
from subprocess import run
from os import popen
from playsound import playsound


class CheckSystem:
    def __init__(self, system):
        self.system = system

    def check_system(self):
        while True:
            if self.system == 'linux':
                print('Good son!')
                playsound('./mp3-files/u-got-that.mp3')
                break
            if self.system != 'linux':
                print('Again!')
                playsound('./mp3-files/error.mp3')
                self.system = input()
            if self.system == 'windows':
                print('Seriously? Do not use this script!')
                playsound('./mp3-files/ultra.mp3')
                break


class CommandFt:
    """This class defines functions that are then 
    executed in bash on an object that has the 
    commands needed to run the fourier transform program"""

    def __init__(self, permissions, compile, program_ft_with_the_data):
        self.permissions = permissions
        self.compile = compile
        self.program_ft_with_the_data = program_ft_with_the_data

    def compile_the_files(self):
        """A function that takes arguments to grant 
        permissions and compile fortran files"""
        return run(self.permissions, shell=True, capture_output=True), run(self.compile, shell=True, capture_output=True)

    def run_ft(self):
        """A function that takes the command argument 
        of the fortran fourier transform program with the data"""
        return run(self.program_ft_with_the_data, shell=True, capture_output=False,)


class CommandNoise:
    """This class defines a noise count function which 
    is then performed in bash for the ft50 data resulting 
    from the calculation of the Fourier transform"""

    def __init__(self, program_noise_with_ft50_data):
        self.program_noise_with_ft50_data = program_noise_with_ft50_data

    def run_noise(self):
        """A function that calculates the noise level 
        using the noise program, which then assigns 
        a sigma 4 value to the variable"""
        output = popen(self.program_noise_with_ft50_data).read()
        # if you want 3 sigma muptiply by 3/4
        noise_string = output[79:90].replace('\n', '')
        noise_4_sigma = float(noise_string)
        return noise_4_sigma


def Fourier_transformation(data, programe='./fortran-files/jkft50 '):
    """A function that executes a class CommandFt for a command object, 
    the result of which is the ft50.trf file with data after 
    the Fourier transform has been computed."""
    #check_system = input('Enter your system (Windows or Linux):')
    # if check_system == 'Linux':
    command = CommandFt('chmod u+x ./compall', 'cd fortran-files && ./compall',
                        programe + data)  # change data if you want
    return command.compile_the_files(), command.run_ft()
    # else:
    #   return print('This version is not compatible with Windows. Use the fortran script and load data with the loader')


def Calculate_the_noise_level():
    """A function that executes a class on a command object 
    that evaluates to the noise level for the ft50.trf data. 
    The range of calculations is from 0 to 50 cycles due to 
    the limitation of the jkft50 program"""
    #check_system = input('Enter your system (Windows or Linux):')
    # if check_system == 'Linux':
    comand_1 = CommandNoise('./fortran-files/ftnoise ./output/ft50.trf 0 50')
    return comand_1.run_noise()
    # else:
    #    return print('This version is not compatible with Windows. Use the fortran script and load data with the loader')


def Check_system():
    comand_2 = CheckSystem(
        input('What system are you using? (Linux or Windows):').lower())
    return comand_2.check_system()
