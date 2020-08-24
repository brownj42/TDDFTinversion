# sharedvalues derived type

This describes each of the values located in the sharedvalues derived type

## Initialization
The arrays located in sharedvalues are allocated using the systemparameters
derived type that has been properly initialized

Fortran
```f90
use derivedtypes
type(sharedvalues) :: sharedvals
call init_sharedvalues(sysparams,sharedvals)
```
```python
sharedvals=TDDFTinversion.derivedtypes.init_sharedvalues(sysparams)
```

## Mandatory values defined by user

### v1

#### Description

The array for the initial external 1-body potential. 

#### Data type
```f90
real(8), dimension(sysparams%ntot1)
```
```python
numpy.empty([sysparams%ntot1],dtype=numpy.float64)
```

### vin

#### Description

The time independant one-body potential

#### Data type
```f90
real(8), dimension(sysparams%ntot1)
```
```python
numpy.empty([sysparams%ntot1],dtype=numpy.float64)
```

### vinteract

#### Description

The 2-body interaction potential

#### Data type

```f90
real(8), dimension(sysparams%ntot2)
```
```python
numpy.empty([sysparams%ntot2],dtype=numpy.float64)
```



