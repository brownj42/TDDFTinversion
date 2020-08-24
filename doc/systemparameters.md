# systemparameters derived type

This describes each of the values located in the systemparameters derived
type.

## Initialization
This derived type (labelled sysparams here) is initialized using np1 (the number of grid points)

Fortran
```f90
type(systemparameters) :: sysparams
integer :: np1
call init_systemparameters(sysparams,np1)
```
Python
```python
np1=numpy.int32( )
sysparams=TDDFTinversion.derivedtypes.init_systemparameters(np1)
```

The systemparameters derived type will now have the np1 value defined

### np1

#### Description
The number of grid points in each dimension

#### Data type
Fortran
```f90
integer 
```
Python
```python
numpy.int32
```

## Mandatory values defined by user in program that calls library

After initialization, user must add the following parameters to the derived
type using commands such as

Fortran
```f90
sysparams%value=...
```
Python
```python
sysparams.value=...
```



### nd

#### Description
The number of dimensions that each particle lives in. For first quantization,
only 1,2, or 3 dimensions is allowed. For second quantization, only 1
dimension is allowed.

#### Data type
Fortran
```f90
integer 
```
Python
```python
numpy.int32
```






### npart

#### Description
The number of particles in the total system


#### Data type
Fortran
```f90
integer
```
Python
```python
numpy.int32
```

### xmin

#### Description
The min position of the grid points that define the lattice with
equal spacing of length n

#### Data type
Fortran
```f90
real(8) 
```
Python
```python
numpy.float64
```

### xmax

#### Description
The min position of the grid points that define the lattice with
equal spacing of length n

#### Data type
Fortran
```f90
real(8) 
```
Python
```python
numpy.float64
```

### dth

#### Description
The goal time step for the program. If an advance of the KS system fails, the
actual time step will be less than this for the next attempt. However, the
program will attempt to increase the time step back to this user defined
timestep on future calls of the program.

#### Data type
Fortran
```f90
real(8) 
```
Python
```python
numpy.float64
```

### T

#### Description
The one dimensional kinetic energy operator. For the sinc basis one can use
Fortran
```f90
call buildkeo(sysparams)
```
Python
```python
TDDFTinversion.keomod.buildkeo(sysparams)
```

after defining xmax and xmin

#### Data type
Fortran
```f90
real(8), dimension(np1,np1)
```
```python
numpy.empty([np1,np1],dtype=numpy.float64)
```

## Values that may be initially defined by the user

These are values that may be initially defined by the user but should not be
changed after the first call to propagate the KS system.

### ct

#### Description
The current time in the propagation. Initialized to 0 by program.

#### Data type
Fortran
```f90
real(8) 
```
Python
```python
numpy.float64
```

### dt

#### Description
The current time-step size. Initially set to dth by program. It is often
beneficial to start the propagation with a smaller time step so as to obtain
the time derivative of the potential more accurately.

#### Data type
Fortran
```f90
real(8) 
```
Python
```python
numpy.float64
```

### quantization

#### Description
If quantization = 1 then first quantization is used. If quantization =2 then
second quantization is used. First quantization can be used for #nd=1,2,3
while second quantization can only be used for nd=1. Default is 1.

#### Data type
Fortran
```f90
integer
```
Python
```python
numpy.int32
```

## Values that may be altered at any time.

### energy

#### Description
The energy of KS system is rescaled to this value by shifting the calculated
potential. If there is no energy imparted on the full system by a
time-dependant potential than this value should be chosen to be some initial
value and then not changed. Otherwise, a good choice may be to measure the
full system's energy at the beginning and end of each time step and use that
value. Default is 0.

#### Data type
Fortran
```f90
real(8)
```
Python
```python
numpy.float64
```

### energynew

#### Description
The energy of KS system after time step dt is rescaled to this value by shifting the calculated
potential. If there is no energy imparted on the full system by a
time-dependant potential than this value should be chosen to be some initial
value and then not changed. Otherwise, a good choice may be to measure the
full system's energy at the beginning and end of each time step and use that
value. Default is 0. After call to fill_systemparameters, value is equal to energy.

#### Data type
Fortran
```f90
real(8)
```
Python
```python
numpy.float64
```

### pinv0minresqlp1

#### Description
The parameter that determines which method is used to invert the
potential. The recommended and initial setting is 1 which uses MINRES-QLP. If
you are in 1D, setting pinv0minresqlp1 will utilize the Morse-Penrose pseudo
inverse.

#### Data type
Fortran
```f90
integer 
```
Python
```python
numpy.int32
```

### singlet

#### Description
The parameter defines the symmetry for a two particle first quantized system with singlet symmetry. 0 means singlet symmetry is not enforced. 1 means singlet symmetry is enforced. Default is 0.

#### Data type
Fortran
```f90
integer 
```
Python
```python
numpy.int32
```

### triplet

#### Description
The parameter defines the symmetry for a two particle first quantized system with triplet symmetry. 0 means triplet symmetry is not enforced. 1 means triplet symmetry is enforced. Default is 0.

#### Data type
Fortran
```f90
integer 
```
Python
```python
numpy.int32
```


### dvksmax

#### Description
The maximum time derivative of the inverted potential at a given point. The
initial value is 1000.

#### Data type
Fortran
```f90
real(8)
```
```python
numpy.float64
```

### newfile

#### Description
The flag that determines whether the output of data will be in a newfile or
appended to the previous files. If set to 1, the old files will be
overwritten. Otherwise, the new data will be appended to the old files.The description of the output files is in the
main [README.md](README.md) newfile is initialized to 1 and then after every
time data is written reverts to 0. 

#### Data type
Fortran
```f90
integer 
```
Python
```python
numpy.int32
```

### outputdata

#### Description
The flag that determines whether the current Kohn-Sham values will be
outputted to data files. If set to 1, the data will be written to
file. Otherwise, no data will be written.

#### Data type
Fortran
```f90
integer 
```
Python
```python
numpy.int32
```


## Values calculated by program

The following are systemparameter values that are accessible to the user but
should not be altered. They are determined by the program by
Fortran
```f90
call fill_systemparameters(sysparams)
```
Python
```python
TDDFTinversion.derivedtypes.fill_systemparameters(sysparams)
```

### ntot1

#### Description

The number of grid points that defines each particle. np1<sup>nd</sup>



#### Data type
Fortran
```f90
integer 
```
Python
```python
numpy.int32
```



### ntot2

#### Description
The number of grid points that defines two particles. np1<sup>2 nd</sup>


#### Data type
Fortran
```f90
integer
```
Python
```python
numpy.int32
```

### ntot

#### Description
The number of grid points that defines the full system wavefunction.   
If in first quantization: ntot=ntot1<sup>npart</sup>  
If in second quantization: ntot=ntot1 choose npart

#### Data type
Fortran
```f90
integer
```
Python
```python
numpy.int32
```

### dx

####
The grid spacing. (xmax-xmin)/(np1-1)

#### Data type
Fortran
```f90
real(8)
```
Python
```python
numpy.float64
```

### xlattice

#### Description
The equally spaced grid points for one dimension defined by xmin,xmax and np1.

#### Data type
Fortran
```f90
real(8), dimension(sysparams%np1)
```
Python
```python
numpy.empty([sysparams%np1,1],dtype=numpy.float64)
```

