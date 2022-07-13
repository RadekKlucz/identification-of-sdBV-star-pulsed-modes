from subprocess import run
from os import popen


class CommandFt:
    """This class defines functions that are then 
    executed in bash on an object that has the 
    commands needed to run the fourier transform program.
    """

    def __init__(self, permissions, compile, program_ft_with_the_data):
        self.permissions = permissions
        self.compile = compile
        self.program_ft_with_the_data = program_ft_with_the_data

    def compile_the_files(self):
        """A function that takes arguments to grant 
        permissions and compile fortran files.
        """

        return run(self.permissions, shell=True, capture_output=True), \
               run(self.compile, shell=True, capture_output=True)

    def run_ft(self):
        """A function that takes the command argument 
        of the fortran fourier transform program with the data.
        """

        return run(self.program_ft_with_the_data, shell=True, capture_output=False)


class CommandNoise:
    """This class defines a noise count function which 
    is then performed in bash for the ft50 data resulting 
    from the calculation of the Fourier transform.
    """

    def __init__(self, program_noise_with_ft50_data):
        self.program_noise_with_ft50_data = program_noise_with_ft50_data

    def run_noise(self):
        """A function that calculates the noise level 
        using the noise program, which then assigns 
        a sigma 4 value to the variable.
        """

        self.output = popen(self.program_noise_with_ft50_data).read()
        self.noise_string = self.output[79:90].replace('\n', '')
        self.noise_4_sigma = float(self.noise_string)  # If you want 3 sigma muptiply by 3/4
        return self.noise_4_sigma


def Fourier_transformation(data, programe='./fortran-files/jkft50 '):
    """A function that executes a class CommandFt for a command object, 
    the result of which is the ft50.trf file with data after 
    the Fourier transform has been computed.
    """

    command = CommandFt('chmod u+x ./compall', 'cd fortran-files && ./compall',
                        programe + data)
    return command.compile_the_files(), command.run_ft()


def Calculate_the_noise_level():
    """A function that executes a class on a command object 
    that evaluates to the noise level for the ft50.trf data. 
    The range of calculations is from 0 to 50 cycles due to 
    the limitation of the jkft50 program.
    """

    comand_1 = CommandNoise('./fortran-files/ftnoise ./output/ft50.trf 0 50')
    return comand_1.run_noise()
