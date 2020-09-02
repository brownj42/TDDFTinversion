# initial_states

This is the module that implements the initialization of the KS orbitals for the example systems given in the paper.

## Subroutines

```f90
!Generates the full system wavefunction for the examples in the !paper.
subroutine initializefullsystem(sysparams,fullvals)
    type(systemparameters), intent(in) :: sysparams
    type(fullvalues), intent(inout) :: fullvals
```

```f90
! returns the one electron ground state in the first KS orbital
! KSvals%phi(:,1)= the one electron ground state given the
! sharedvals%v1(:) potential and sysparams%T(:,:) KEO
subroutine oneelectrongroundstate(sysparams,sharedvals,KSvals)
    use derivedtypes
    use hamiltonian_mod
    type(systemparameters), intent(in) :: sysparams
    type(sharedvalues), intent(in) :: sharedvals
    type(KSvalues), intent(out) :: KSvals
```

```f90
!Initializes the KSvals%phi orbitals given the parameters in
!sysparams and sharedvals to match the full electron density dpe
!fullvals, optional used for case when only one electron by
!setting KSvals%phi(:,1)=fullvals%psi(:)
subroutine initializeKSsystem(sysparams,sharedvals,dpe,fullvals,KSvals)
    use derivedtypes
    use density_mod
    type(systemparameters), intent(in) :: sysparams
    real(8), intent(in) :: dpe(:)
    type(sharedvalues), intent(in) :: sharedvals
    type(KSvalues), intent(out) :: KSvals
    type(fullvalues), intent(in), optional :: fullvals
```

```f90
!calculates all eigenstates of a 1D Hamiltonian with KEO matrix T and
!diagonal potential v of length np1
subroutine calcalleigenstates(np1,T,v,eigenstates)
     integer, intent(in) :: np1
     real(8), intent(in) :: T(:,:)
     real(8), intent(in) :: v(:)
     real(8), intent(out) :: eigenstates(np1,np1)
```
