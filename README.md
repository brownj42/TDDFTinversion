# TDDFT-inversion program

This program is the software that performs the numerical TDDFT potential
inversion technique outline in [Solver for the electronic V-representation
problem of time-dependent density functional
theory](https://arxiv.org/abs/1904.10958). The main code is written in Fortran
but can be called from Python. Example main programs are available for Fortran
in the [fortran](/fortran/) folder and for Python 3 in the [python](/python/)
folder.

# Dependancies

The program needs
* [ARPACK](https://www.caam.rice.edu/software/ARPACK/)
* BLAS (OpenBLAS with OpenMP parallization or MKL are good choices) 
* LAPACK
* Fortran compiler (gfortran or ifort)

If you wish to use the Python wrapper
* C compiler (gcc or icc)
* C++ compiler (g++ or ipp)
* Python 3.6+
* Recent version of numpy which includes f2py
* [f90wrap](https://github.com/jameskermode/f90wrap)

# Compilation

 Inside the [python](/python/) and [fortran](/fortran/) folders are makefiles
 that should compile the respective programs. You may have to edit the
 Makefiles to point the appropriate ARPACK BLAS and LAPACK libary object
 files. The makefiles were tested with the standard libary locations on an
 OpenSuse desktop.

# Documentation 

There is some information on the use of the software in the
[doc](/doc/) folder
