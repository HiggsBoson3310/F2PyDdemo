# F2PyDdemo
This is a simple example of how to use Fortran code in Python by using the numpy package and the F2Py module.

The code computes the Rayleigh and Raman scattering cross sections (in cm^2) for Hydrogenic states. It does it by the Dalgarno-Lewis method. The calculation follows Sadeghpour et a [1]. And the details of the calculation can be found in the notes in this repo.

The following instructions allow you to run the code in a LINUX (UNIX) machine with Python 3 installed and with the current version of numpy (1.21.1) and scipy (1.7.0).
For the Halstead machine in Purdue this is loaded by using the commands 

```
module load use.own
module load anaconda/2020.11-py38
```
Make sure you have the gnu compilers activated, since there are some compiling issues with the default intel compilers that we are still trying to solve. 

To do so in the Purdue machine run

```
module load gcc
```

After downloading this repo and decompressing it go to the main directory, with the FORTRAN files (COUL90.FOR, COULN.FOR and MAABgensub.f) and the auxliary signature files (coul90.pyf, couln.pyf and gensub.pyf)open a terminal and run the following:

```
f2py -c coul90.pyf COUL90.FOR --quiet
f2py -c couln.pyf COULN.FOR --quiet
f2py -c gensub.pyf MAABgensub.f --quiet
```

The quiet command prevents the compiler from issuing end of line warnings. This generates three .so files that contain the packages (or function wrappers) to be imported by python.

To verify the installation of the packages run:

```
python
import gensub
import coul92
import couln
```

It should run smoothly. The .pyf signature files contain the definition of the subroutines and functions in each one of the fortran source files (which can then be called as python functions) and the variables that are used in with identification as an input or output variable. The syntax of the pyf files is similar to Fortran 90. 

To explore the defined functions and the requiered variables run the following right after the import lines while still in python:

```
print(gensub.seaton.__seaton__)
```

This will print the necessary input for the function and the output. There is no description, but the notation is the same from the source fortran codes. 

To exit python back into the linux terminal type 

```
exit()
```

To generate data, go into the HydrogenPolarizability.py file and write the name of the test to run at the end. There is a section of the code that describes all the programmed tests.

Just type the name of the test to run. For instance to generate the plots to compare with MacNamara et al [2] we want to run test4. In the HydrogenPolarizability.py file go to the end and type 

```
test4()
```

After saving, go to the terminal in the folder of the repo and type 

```
python HydrogenPolarizability.py
```
this takes a while to run, but generates the cross sections in a dense mesh and stores them in the "CrossSections" directory. This generates files with the different cross sections. The first column is always the photon energy (x axis) and the other columns are explained in the code. The files preloaded in the directory correspond to the data from the MacNamara paper. The name of the files indicates the involved states. 

Othr tests involve sanity checks for the Dalgarno Lewis radial function (test0), Table 1 from the Sadeghpour et al paper cited above (test1), a test on the continuity of the Greens function in energy for fixed values of r and r prime (testgf).

References:
[1] Sadeghpour, H. R., & Dalgarno, A. (1992). Rayleigh and Raman Scattering by Hydrogen and Caesium. Journal of Physics B: Atomic, Molecular and Optical Physics, 25(22), 4801–4809
[2] McNamara et al. (2018). Efficient calculation of Rayleigh and Raman scattering. Physical Review A, 98(4) 
