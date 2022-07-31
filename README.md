# Identification-of-sdBV-star-pulsed-modes

This is a data project that has been used in astronomical research that has helped identify the modes of stellar pulses. This Python project automates things that humans have done manually in the past. A brief theoretical introduction is available in the report catalog in English and Polish.

In the basic form of the project, the fortran language was used to calculate the Fourier transform to identify the modules. The new version has the ability to calculate the Fourier transform in a python without a fortran. This method is recommended, but a first version is also available.

## 🧑‍💻 Technology stack

* Docker
* Jupiter notebook
* Python
  * pandas
  * numpy
  * matplotlib
  * scipy
  * astropy
* Fortran
* LaTeX

## 🔥 How to run with Docker

1. Clone the repository to your computer,
2. Open PowerShell/bash or another terminal at the folder with project and run the following command:

```cmd
docker compose up --build
```

3. Copy the last link from the terminal to your favorite browser to run jupiter notebook,
4. Enjoy the data project.

### 🐋 Container

Base Jupyter Notebook Stack -- python environment to run all files and modules

## 🎆 How to run without Docker

This method is recommended for people who prefer to use fortran files to calculate Fourier transform and noise using fortran files. There is also a method with no fortran files, requiring only Python modules. The steps are bellow:

1. Install all required modules in Python,
2. Install gfortran on your computer if you want to use fortran (optional),
3. Clone the repository to your computer,
4. Open main.ipynb in your favorite text editor and run it with the specified p (python) or f (fortran) method.

## 🌠 Features

🌟 **Docker containerisation and orchestration**

🌟 **Natural Language Processing**

🌟 **Data Manipulation with pandas and numpy**

🌟 **Visualizing data with a matplotlib**

## 📁 Directory Structure

    ├───main
    │   ├───data
    │   ├───fortran-files
    │   └───output
    └───raports
        └───latex

## 📧 Contact

[![LinkedIn](https://i.stack.imgur.com/gVE0j.png) Radosław Kluczewski](https:///www.linkedin.com/in/radoslaw-kluczewski) 
&nbsp;
[![GitHub](https://i.stack.imgur.com/tskMh.png) RadekKlucz](https://github.com/RadekKlucz)

## License

[![Licence](https://img.shields.io/github/license/Ileriayo/markdown-badges?style=for-the-badge)](./LICENSE)
