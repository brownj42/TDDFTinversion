# propagate

This is the module that is used to advance individual quantum orbitals. There is both a time-independant and time dependant propagator.

## Subroutines

```f90
!the subroutine used to propagate the wavefunction by dt given in sysparams%dt
subroutine advancewf(sysparams,sharedvals,niter,vks,yin,yout)
    type(systemparameters), intent(in) :: sysparams
    type(sharedvalues), intent(in) :: sharedvals
    integer, intent(in)  :: niter !number of matrix-vector products for Lanczos
    real(8), intent(in) :: Vks(:) !Potential vector
    complex(8), intent(in) :: yin(:) !input wavefunction at time t
    complex(8), intent(out) :: yout(:) !output wavefunction at time t+dt
```

```f90
!subroutine used to advance the wavefunction yin by sysparams%dt to yout
!The time derivative of the potential at t and t+dt are optional, if they are present, a higher order
!method is used to advance the system in time.
subroutine advancewftd(sysparams,sharedvals,niter,vin,vks,vksnew,dvks,dvksnew,yin,yout)
   use derivedtypes
   use hamiltonian_mod
   type(systemparameters), intent(in) :: sysparams
   type(sharedvalues), intent(in) :: sharedvals
   integer, intent(in)  :: niter !number of iterations used for Short iterative Lanczos
   real(8), intent(in) :: vin(:) !time independant part of potential
   real(8), intent(in) :: Vks(:) !time-dependant part of potential at time t
   real(8), intent(in) :: vksnew(:) !time-dependant part of potential at time t+dt
   real(8), intent(in), optional :: dvks(:) !time derivative of potential at time t
   real(8), intent(in), optional :: dvksnew(:) !time derivative of potential at time t+dt
```

```f90
subroutine advanceKSsystem(dpe,dpenew,dnx,ddnx,ddnxnew,sysparams,KSvals,sharedvals,info)
!The subroutine that implements that main method to advance the KSorbitals and calculate the 
!corresponding TDDFT potential. The five vectors required to advance from ct to ct+dt are:
!dnx: The time derivative of the density at time t
!dpe: The density at time t
!dpenew: The density at time t+dt
!ddnx: The second time-derivative of the density at time t
!ddnxnew: The second time-derivative of the density at time t+dt
!if after attempted time step, info=1 the orbitals have advanced and vks at time t+dt is stored in
!KSvals along with the new orbitals
!If after attempted time step, info=0 The time step failed and new dpenew and ddnxnew are requested at
!time t+dt where dt is now halved from the previous step.
   real(8), intent(in) :: dnx(:),dpe(:),dpenew(:),ddnx(:),ddnxnew(:)
   type(systemparameters), intent(inout) :: sysparams
   type(KSvalues), intent(inout) :: KSvals
   type(sharedvalues), intent(in) :: sharedvals
   integer, intent(out) :: info
```
