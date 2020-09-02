# ksvalues derived type

The derived type that stores all vectors related to the KS system

## Initialization

If you would like to initialize the orbitals using your own method, the arrays
can be allocated using

Fortran
```f90
type(ksvalues) :: ksvals
call init_ksvalues(sysparams,ksvals)
```
Python
```python
ksvals=TDDFTinversion.derivedtypes.init_ksvalues(sysparams)
```

If you would like to allow the program to initialize the KS orbitals
Fortran
```f90
type(ksvalues) :: ksvals
call initializekssystem(sysparams,sharevals,exactdensity,fullvals,ksvals)
```
Python
```python
ksvals=TDDFTinversion.initial_states.initializeksystem(sysparams,sharedvals,exactdensity,fullvals)
```

where exactdensity is a real(8) or numpy.float64 vector of length
sysparams%ntot1 or sysparams.ntot1 which has the exact density and fullvals is
the fullvalues derived type.

## Mandatory variables defined by user

### phi

#### Description
The ntot1 by npart array that stores the orbitals of the KS system. Before
beginning the propagation of the KS system. These orbitals should be defined
by the user such that the density and derivative of the density match the
target system. If orbitals are defined such that the density matches but the
derivative of the density does not. The program will try to assign the phases
to correct this but this will not work well unless the initial derivative of
the density is close to the target derivative of the density.

#### Data type
Fortran
```f90
complex(8), dimension(ntot1,npart)
```
Python
```python
numpy.empty([ntot1,npart],dtype=numpy.complex128)
```


## Variables calculated by program

### dp

#### Description
The density of the KS system at time systemparameters%ct

#### Data type
Fortran
```f90
real(8), dimension(ntot1)
```
Python
```python
numpy.empty([ntot1,npart],dtype=numpy.float64)
```

### vks

#### Description
The full 1-body KS potential at time systemparameters%ct

#### Data type
Fortran
```f90
real(8), dimension(ntot1)
```
Python
```python
numpy.empty([ntot1,npart],dtype=numpy.float64)
```

### vksh

#### Description
The full 1-body KS potential with sharedvalues%v1 subtracted off for time from
the previous successful timestep

#### Data type
Fortran
```f90
real(8), dimension(ntot1)
```
Python
```python
numpy.empty([ntot1],dtype=numpy.float64)
```

### vkshh

#### Description
The full 1-body KS potential with sharedvalues%v1 subtracted off for time from
the previous to previous successful timestep

#### Data type
Fortran
```f90
real(8), dimension(ntot1)
```
Python
```python
numpy.empty([ntot1],dtype=numpy.float64)
```

### vhar

#### Description
The hartree potential calculated using dp and sharedvalues%vinteract at the
current time

### vksh

#### Description
The full 1-body KS potential with sharedvalues%v1 subtracted off for time from
the previous successful timestep

#### Data type
Fortran
```f90
real(8), dimension(ntot1)
```
Python
```python
numpy.empty([ntot1],dtype=numpy.float64)
```

### dvks

#### Description
The current estimate of the derivative of the vks potential with respect to
time

#### Data type
Fortran
```f90
real(8), dimension(ntot1)
```
Python
```python
numpy.empty([ntot1],dtype=numpy.float64)
```

### dvksh

#### Description
The previous to current time derivative of the vks potential.

#### Data type
Fortran
```f90
real(8), dimension(ntot1)
```
Python
```python
numpy.empty([ntot1],dtype=numpy.float64)
```

