module derivedtypes 
  implicit none
  type systemparameters
     integer :: nd  !number of dimensions each particle lives in
     integer :: np1 !number of grid points in each dimension
     integer :: npart !number of particles
     integer :: ntot1 !number of grid points for each particle
     integer :: ntot2 !number of grid points for two particles
     integer :: ntot  !number of grid points for full system
     real(8) :: dx !grid spacing
     real(8) :: ct !current time
     real(8) :: dt !current time step
     real(8) :: dth !goal time step	
     real(8) :: xmax,xmin !max grid value, min grid value
     real(8) :: energy !energy of KS system at ct
     real(8) :: energynew !energy of KSsystem at time ct+dt
     integer :: dvks_sw !switch for propagator, don't change
     integer :: dvksh_sw !switch for propagator, don't change
     real(8) :: dvksmax !max time-derivative of potential
     integer :: pinv0minresqlp1 !which inversion algorithm
     integer :: quantization !first =1 or second =2 quantization
     integer :: newfile !generate new output file
     integer :: outputdata !output data after successful time step
     integer :: loop !don't change, internal counter
     integer :: occupy_case  !used for the initialization of the KS orbitals 
     integer :: singlet,triplet !singlet,triplet=1 turns on singlet symmetry
     real(8), allocatable :: xlattice(:) !1D grid points
     real(8), allocatable :: T(:,:) !1D kinetic energy operator
  end type systemparameters

  type ksvalues
     complex(8), allocatable :: phi(:,:) !KS orbtial
     real(8), allocatable :: vks(:) !current Vks potential
     real(8), allocatable :: vksh(:) !previous Vks potential - external 1-body potential
     real(8), allocatable :: vkshh(:) !second previous Vks potential - 1-body potential
     real(8), allocatable :: vhar(:) !hartree potential
     real(8), allocatable :: dvks(:) !time-derivative of potential
     real(8), allocatable :: dvksh(:) !previous time-derivative of potential
     real(8), allocatable :: dp(:) !current KS density
  end type ksvalues

  type fullvalues
     complex(8), allocatable :: psi(:) !full n-body wavefunction size sysparams%ntot
     real(8), allocatable :: v(:) !2-body or full n-body potential depending on quantization
  end type fullvalues
  
  type sharedvalues
     real(8), allocatable :: v1(:) !1-body external potential
     real(8), allocatable :: vin(:) !time-independant 1-body external potential
     real(8), allocatable :: vinteract(:) !interaction potential
  end type sharedvalues

contains
  
    subroutine init_systemparameters(sysparams, np1)
       type(systemparameters), intent(out) :: sysparams
       INTEGER(4), INTENT(in) :: np1
       sysparams%np1=np1
       sysparams%dvksmax=1000.d0
       sysparams%ct=0.d0
       sysparams%xmin=0.d0
       sysparams%xmax=dble(np1)
       sysparams%energy=0.d0
       sysparams%energynew=0.d0
       sysparams%pinv0minresqlp1=1
       sysparams%quantization=1
       sysparams%newfile=1
       sysparams%outputdata=1
       sysparams%loop=1
       sysparams%dvks_sw=1
       sysparams%dvksh_sw=1
       sysparams%occupy_case=0
       sysparams%singlet=0
       sysparams%triplet=0 
       ALLOCATE(sysparams%T(np1,np1),sysparams%xlattice(np1))
    end subroutine init_systemparameters

    subroutine destroy_systemparameters(sysparams)
       type(systemparameters), INTENT(inout) :: sysparams
       if (allocated(sysparams%T)) deallocate(sysparams%T)
       if (allocated(sysparams%xlattice)) deallocate(sysparams%xlattice)
    end subroutine destroy_systemparameters

    subroutine fill_systemparameters(sysparams)
       use nchoosekmod
       type(systemparameters), intent(inout) :: sysparams
       integer :: jin,prev
       sysparams%ntot1=sysparams%np1**sysparams%nd
       sysparams%ntot2=sysparams%ntot1**min(sysparams%npart,2)
       if (sysparams%quantization==1) then
          sysparams%ntot=sysparams%ntot1**sysparams%npart
       else
          jin=0
          prev=0 
          sysparams%ntot=0
          call countnchoosek(sysparams%ntot1,sysparams%npart,jin,prev,sysparams%ntot)
       end if
       sysparams%energynew=sysparams%energy
       sysparams%dx=(sysparams%xmax-sysparams%xmin)/(sysparams%np1-1)
       do jin=1,sysparams%np1
          sysparams%xlattice(jin)=(jin-1)*sysparams%dx+sysparams%xmin
       end do
       sysparams%dt=sysparams%dth
    end subroutine fill_systemparameters

    subroutine init_sharedvalues(sysparams,sharedvals)
        type(systemparameters), intent(in) :: sysparams
        type(sharedvalues), intent(out) :: sharedvals
        integer :: ntot1,ntot2
        ntot1=sysparams%ntot1
        ntot2=sysparams%ntot2
        allocate(sharedvals%v1(ntot1))
        allocate(sharedvals%vin(ntot1))
        allocate(sharedvals%vinteract(ntot2))
     end subroutine init_sharedvalues

     subroutine destroy_sharedvalues(sharedvals)
        type(sharedvalues), intent(inout) :: sharedvals
        if (allocated(sharedvals%v1)) deallocate(sharedvals%v1)
        if (allocated(sharedvals%vin)) deallocate(sharedvals%vin)
        if (allocated(sharedvals%vinteract)) deallocate(sharedvals%vinteract)
     end subroutine destroy_sharedvalues

     subroutine init_fullvalues(sysparams,fullvals)
        type(systemparameters), intent(in) :: sysparams
        type(fullvalues), intent(out) :: fullvals
        integer :: ntot,ntot2
        ntot=sysparams%ntot
        if (sysparams%quantization==1) then
           allocate(fullvals%v(ntot))
           allocate(fullvals%psi(ntot))
        else
           ntot2=sysparams%ntot2
           allocate(fullvals%v(ntot2))
           allocate(fullvals%psi(ntot))
        end if
     end subroutine init_fullvalues

     subroutine destroy_fullvalues(fullvals)
        type(fullvalues), intent(inout) :: fullvals
        if (allocated(fullvals%v)) deallocate(fullvals%v)
        if (allocated(fullvals%psi)) deallocate(fullvals%psi)
     end subroutine destroy_fullvalues

     subroutine init_ksvalues(sysparams,ksvals)
        type(systemparameters), intent(in) :: sysparams
        type(ksvalues), intent(out) :: ksvals
        integer :: ntot1,npart
        ntot1=sysparams%ntot1
        npart=sysparams%npart
        allocate(KSvals%phi(ntot1,npart))
        allocate(KSvals%vks(ntot1),KSvals%vksh(ntot1),KSvals%vkshh(ntot1))
        allocate(KSvals%vhar(ntot1),KSvals%dvks(ntot1),KSvals%dvksh(ntot1))
        allocate(KSvals%dp(ntot1))
        KSvals%vks=0.d0
        KSvals%vksh=0.d0
        KSvals%vkshh=0.d0
        KSvals%vhar=0.d0
        KSvals%dvks=0.d0
        KSvals%dvksh=0.d0
     end subroutine init_ksvalues

     subroutine destroy_ksvalues(ksvals)
        type(ksvalues), intent(inout) :: ksvals
        if (allocated(ksvals%phi)) deallocate(ksvals%phi)
        if (allocated(ksvals%vks)) deallocate(ksvals%vks)
        if (allocated(ksvals%vksh)) deallocate(ksvals%vksh)
        if (allocated(ksvals%vkshh)) deallocate(ksvals%vkshh)
        if (allocated(ksvals%vhar)) deallocate(ksvals%vhar)
        if (allocated(ksvals%dvks)) deallocate(ksvals%dvks)
        if (allocated(ksvals%dvksh)) deallocate(ksvals%dvksh)
        if (allocated(ksvals%dp)) deallocate(ksvals%dp)
     end subroutine destroy_ksvalues

end module derivedtypes
