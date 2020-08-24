Program tddftinversion
  use derivedtypes 
  use readin_data_mod
  use propagate 
  use density
  use initial_states
  use potential 
  use KSadvance_mod
  implicit none
  Complex(8), allocatable :: psinew(:)
  
  real(8), allocatable :: dpe(:),dpenew(:)
  real(8), allocatable :: ddnx(:),dnx(:),ddnxnew(:)
  Integer :: npart,ntot1,nadvancesteps
  real(8) :: ct,dt
  integer :: loop
  integer :: info
  character(50) :: filename
  integer :: niter
  type(systemparameters) :: sysparams
  type(sharedvalues) :: sharedvals
  type(fullvalues) :: fullvals
  type(KSvalues) :: KSvals
  filename='parameters.in' 
  call readparameters(readin=.True.,filename=filename,sysparams=sysparams)
  
  nadvancesteps=1200
 
  
  ntot1=sysparams%ntot1  
  npart=sysparams%npart
 
  !allocate a bunch of vectors 
  allocate(dnx(sysparams%ntot1),ddnx(sysparams%ntot1),psinew(sysparams%ntot)) 
  allocate(ddnxnew(sysparams%ntot1))

  call init_sharedvalues(sysparams,sharedvals)
  call generate_1bodypot(sysparams,sharedvals)
  call generate_2bodypot(sysparams,sharedvals)

  call init_fullvalues(sysparams,fullvals)
  call generate_nbodypot(sysparams,sharedvals,fullvals)

  call initializefullsystem(sysparams,fullvals)
 
  allocate(dpe(sysparams%ntot1),dpenew(sysparams%ntot1)) 
  !Calculate target density for initial state
  call fullwf_density(sysparams,fullvals%psi,dpe)

  call initializeKSsystem(sysparams,sharedvals,dpe,fullvals,KSvals)
 
  if (npart==2) then
     !add driving potential
     call add_driving_potential(sysparams,sharedvals,fullvals)
  end if

  !initialize time to 0
  ct=0.d0
  sysparams%ct=0.d0

  !krylov subspace size for propagation
  niter=30
   

  !Calculate 2nd derivative of density for initial state 
  !Remove known time-independant potential
  doloop:do loop=1,nadvancesteps
     ct=sysparams%ct
     dt=sysparams%dt
     print*,' '
     print*,'Starting step for time',ct
     print*,'Time step size',dt
     

     !calculate 1st derivative of density
     call calcdnx(sysparams,ntot1,fullvals%psi,fullvals%v,dnx)
     ! calculate full wf density
     call fullwf_density(sysparams,fullvals%psi,dpe)
     !calculate 2nd derivative of density
     call calcddnx(sysparams,ntot1,fullvals%psi,fullvals%v,ddnx)
     !propagate full npart system and calculate new density
     call advancewf(sysparams,sharedvals,niter,fullvals%v,fullvals%psi,psinew)
     !calculate new fullwf density
     call fullwf_density(sysparams,psinew,dpenew)
     !calculate new second derivative of full system reduced density
     call calcddnx(sysparams,ntot1,psinew,fullvals%v,ddnxnew)
     call advanceKSsystem(dpe,dpenew,dnx,ddnx,ddnxnew,loop,sysparams,KSvals,sharedvals,info)
    
     if (info==1) then
        fullvals%psi=psinew
     end if 
  end do doloop


end program tddftinversion
