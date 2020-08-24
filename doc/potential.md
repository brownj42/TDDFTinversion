# potential

This is the module that contains the subroutines that generate the example potentials. It is unnecessary to call any of these functions to run the program.

## Subroutines

```f90
!Generates the 1body potential given different input parameters.
subroutine generate_1bodypot(sysparams,sharedvals)
    use derivedtypes
    type(systemparameters), intent(in) :: sysparams
    type(sharedvalues), intent(inout) :: sharedvals
```

```f90
!Generates the 2body potential for different input parameters
subroutine generate_2bodypot(sysparams,sharedvals)
    use derivedtypes
    type(systemparameters), intent(in) :: sysparams
    type(sharedvalues), intent(inout) :: sharedvals
```

```f90
!Generates the full n_body wavefunction, really only used for the 1D 2-electron example
subroutine generate_nbodypot(sysparams,sharedvals,fullvals)
    use derivedtypes
    type(systemparameters), intent(in) :: sysparams
    type(sharedvalues), intent(in) :: sharedvals
    type(fullvalues), intent(inout) :: fullvals
```

```f90
!adds the driving potential for the 1D 2-electron example
subroutine add_driving_potential(sysparams,sharedvals,fullvals)
    use derivedtypes
    type(systemparameters), intent(in) :: sysparams
    type(sharedvalues), intent(inout) :: sharedvals
    type(fullvalues), intent(inout) :: fullvals
```
