# density

The module that calculates the density and its first and second time-derivative.

## subroutines

```f90
subroutine fullwf_density(sysparams,yh,density)
    !Calculates the fullwavefunction density given systemparameters derived type and wavefunction
    !yh with size sysparams%ntot
    type(systemparameters), intent(in) :: sysparams
    complex(8), intent(in) :: yh(:)
    real(8), intent(out) :: density(sysparams%ntot1)
```

```f90
subroutine calcddnx(sysparams,sharedvals,ntot1,psi,v,ddnx) 
    !Calculates the second time derivative of the density ddnx of length sysparams%ntot1
    !from psi of length sysparams%ntot
    !with potential v 
    type(systemparameters), intent(in) :: sysparams
    type(sharedvalues), intent(in) :: sharedvals 
    integer, intent(in) :: ntot1
    real(8), intent(in) :: v(:)
    complex(8), intent(in) :: psi(:)
    real(8), intent(out) :: ddnx(sysparams%ntot1)
```

```f90
subroutine calcdnx(sysparams,sharedvals,ntot1,psi,v,dnx)
    !Calculates the second time derivative of the density dnx of length sysparams%ntot1
    !from psi of length sysparams%ntot
    !with potential v
    type(systemparameters), intent(in) :: sysparams
    type(sharedvalues), intent(in) :: sharedvals 
    integer, intent(in) :: ntot1
    real(8), intent(in) :: v(:)
    complex(8), intent(in) :: psi(:)
```