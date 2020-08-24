# derivedtypes  

This is the module that contains the definitions of the derived types and the subroutines that initialize and destroy the allocatable arrays therein.  

## Derived types

There are four derived types that are used in the program.

* [systemparameters](systemparameters.md): This contains all the information that defines the system
being studied (mandatory)
* [ksvalues](ksvalues.md): This contains all the information about the Kohn-Sham system (mandatory)
* [sharedvalues](sharedvalues.md): This contains the one-body potential, time-independant 1-body (mandatory)
potential, the 2-body interaction potential. (mandatory)
* [fullvalues](fullvalues.md): This contains the full n-body wavefunction and the n-body (optional)
potential

## Subroutines

### Initialization  

This subroutine allocates the arrays for the kinetic energy operator and lattice grid values in the systemparameters derived type. The number of 1d grid points is the necessary input. At this time, all dimensions have the same number of grid points.

```f90
subroutine init_systemparameters(sysparams,number_of_1d_points)
    type(systemparameters), intent(out) :: sysparams
    integer, intent(in) :: number_of_1d_points

```  

This subroutine initializes the arrays in the ksvalues (Kohn-Sham values) derived type.

```f90
subroutine init_ksvalues(sysparams,ksvals)
    type(systemparameters), intent(in) :: sysparams
    type(ksvalues), intent(out) :: ksvals
```  

This subroutine initializes the arrays in the sharedvalues derived type.

```f90
subroutine init_sharedvalues(sysparams,sharedvals)
    type(systemparameters), intent(in) :: sysparams
    type(sharedvalues), intent(out) :: sharedvals
```  

This subroutine initializes the non-user defined values in the systemparameters derived type

```f90
subroutine fill_systemparameters(sysparams)
    use nchoosekmod
    type(systemparameters), intent(inout) :: sysparams
```

This subroutine initializes the arrays in the fullvalues derived type

```f90
subroutine init_fullvalues(sysparams,fullvals)
        type(systemparameters), intent(in) :: sysparams
        type(fullvalues), intent(out) :: fullvals
```  

### Destroy subroutines (releases allocated memory)

```f90
!Destroys the arrays in systemparameters derived types
subroutine destroy_systemparameters(sysparams)
    type(systemparameters), INTENT(inout) :: sysparams
```

```f90
!Destroys the arrays in the sharedvalues derived types
subroutine destroy_sharedvalues(sharedvals)
    type(sharedvalues), intent(inout) :: sharedvals
```

```f90
!Destroys the arrays in the fullvalues derived type
subroutine destroy_fullvalues(fullvals)
    type(fullvalues), intent(inout) :: fullvals
```

```f90
!Destroys the arrays in the ksvalues derived types
subroutine destroy_ksvalues(ksvals)
    type(ksvalues), intent(inout) :: ksvals
```
