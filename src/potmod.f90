Module potential
  implicit none
contains
  !Generates the 1body potential given different input parameters.
  subroutine generate_1bodypot(sysparams,sharedvals)
    use derivedtypes 
    type(systemparameters), intent(in) :: sysparams
    type(sharedvalues), intent(inout) :: sharedvals
    integer :: i,nd,npart,np1,ntot1,ntot2,ntot
    nd=sysparams%nd
    np1=sysparams%np1
    ntot1=sysparams%ntot1
    ntot2=sysparams%ntot2
    ntot=sysparams%ntot
    npart=sysparams%npart
   
    if (sysparams%quantization==1) then 
       !generate external potential
       if (nd==1.and.npart==2) then
          sharedvals%v1=pot1d_2(np1,sysparams%xlattice)
       elseif (npart==1.and.nd==1) then
          sharedvals%v1=pot1d(np1,sysparams%xlattice)
       elseif (nd==3.and.npart==1) then
          sharedvals%v1=pot3d(np1,ntot1,sysparams%xlattice)
       end if
       !Initialize time-independant 1-body potential
       if (npart>1) then
          sharedvals%vin=sharedvals%v1
       else
          sharedvals%vin=0.d0
       end if
    else
       do i=1,ntot1
          sharedvals%v1(i)=dble(i-(ntot1+1)/2)**2/20.d0!1.d0
       end do
       sharedvals%vin=sharedvals%v1
    end if

  end subroutine generate_1bodypot

  !Generates the 2body potential for different input parameters
  subroutine generate_2bodypot(sysparams,sharedvals,sc)
    use derivedtypes
    type(systemparameters), intent(in) :: sysparams
    type(sharedvalues), intent(inout) :: sharedvals
    real(8), intent(in), optional :: sc
    integer :: nd,npart,np1,ntot1,ntot2,ntot
    integer :: i,j
    nd=sysparams%nd
    np1=sysparams%np1
    ntot1=sysparams%ntot1
    ntot2=sysparams%ntot2
    ntot=sysparams%ntot
    npart=sysparams%npart
       
    
    if (sysparams%quantization==1) then
       if (npart==2.and.nd==1) then
          !generate full two-body potential
          if (present(sc)) then
             sharedvals%vinteract=pot_interaction(np1,ntot2,sysparams%xlattice,sc)
          else
             sharedvals%vinteract=pot_interaction(np1,ntot2,sysparams%xlattice,0.1d0) 
          end if
       elseif(npart==1) then
          sharedvals%vinteract=0.d0
       end if
    else
       sharedvals%vinteract=0.d0
       do i=1,ntot1
          do j=1,ntot1
             if (i==j+1) then
                sharedvals%Vinteract((i-1)*np1+j)=1.d0
             elseif(i+1==j) then
                sharedvals%Vinteract((i-1)*np1+j)=1.d0
             elseif(i+np1-1==j) then
                sharedvals%Vinteract((i-1)*np1+j)=0.0d0
             elseif(i==j+np1-1) then
                sharedvals%Vinteract((i-1)*np1+j)=0.0d0
             end if
          end do
       end do
    end if
  end subroutine generate_2bodypot

  !Generates the full n_body wavefunction, really only used for the 1D 2-electron example
  subroutine generate_nbodypot(sysparams,sharedvals,fullvals)
    use derivedtypes
    type(systemparameters), intent(in) :: sysparams
    type(sharedvalues), intent(in) :: sharedvals
    type(fullvalues), intent(inout) :: fullvals
    integer :: nd,npart,np1,ntot1,ntot2,ntot
    nd=sysparams%nd
    np1=sysparams%np1
    ntot1=sysparams%ntot1
    ntot2=sysparams%ntot2
    ntot=sysparams%ntot
    npart=sysparams%npart
    
    if (sysparams%quantization==1) then 
       if (npart==2.and.nd==1) then
          !generate full two-body potential
          call pot2particle(ntot1,ntot2,sharedvals%v1,sharedvals%vinteract,fullvals%v)
       elseif(npart==1) then
          fullvals%v=sharedvals%v1
       end if
    else
       fullvals%v=sharedvals%vinteract
    end if
  end subroutine generate_nbodypot
    
  
  function pot3d(np1,ntot,xlattice) result(v)
    integer :: np1,ntot
    real(8) :: xlattice(np1)
    real(8) :: v(ntot)
    integer :: i,j,k,l
    l=0
    do i=1,np1
       do j=1,np1
          do k=1,np1
             l=l+1
             v(l)=0.5d0*(xlattice(i)**2+xlattice(j)**2+xlattice(k)**2)
             v(l)=v(l)+0.025d0*(xlattice(i)*xlattice(j)+&
                  2.d0*xlattice(k)*xlattice(i)+&
                  3.d0*xlattice(j)*xlattice(k))
          end do
       end do
    end do
  end function pot3d

  function pot1d(np1,xlattice) result(v)
    integer :: np1
    real(8) :: xlattice(np1)
    real(8) :: v(np1)
    integer :: i
    do i=1,np1
       v(i)=1.d0/2.d0*xlattice(i)**2
    end do
  end function pot1d

  function pot1d_2(np1,xlattice) result(v1)
    integer :: np1
    real(8) :: xlattice(np1)
    real(8) :: v1(np1)
    integer :: i
    real(8) :: alpha,beta
    alpha=5.d-11
    beta=1.3d-4
    do i=1,np1
       v1(i)=alpha*xlattice(i)**4*(xlattice(i)**6-beta/alpha)-0.8d0
    end do
  end function pot1d_2

  function pot_interaction(np1,ntot2,xlattice,sc) result(vinteract)
    integer :: np1,ntot2
    real(8) :: xlattice(np1)
    real(8) :: vinteract(ntot2)
    real(8) :: sc
    integer :: i,j,l
    l=0
    do i=1,np1
       do j=1,np1
          l=l+1
          vinteract(l)=1.d0/sqrt((xlattice(i)-xlattice(j))**2+sc)
       end do
    end do
  end function pot_interaction

  subroutine pot2particle(ntot1,ntot2,v1,vinteract,v)
    integer, intent(in) :: ntot1,ntot2
    real(8), intent(in) :: v1(ntot1),vinteract(ntot2)
    real(8), intent(out) :: v(ntot2)
    integer :: i,j,l
    l=0
    do i=1,ntot1
       do j=1,ntot1
          l=l+1
          v(l)=v1(i)+v1(j)+vinteract(l)
       end do
    end do
  end subroutine pot2particle
  
  !adds the driving potential for the 1D 2-electron example
  subroutine add_driving_potential(sysparams,sharedvals,fullvals,prefactor)
    use derivedtypes 
    type(systemparameters), intent(in) :: sysparams
    type(sharedvalues), intent(inout) :: sharedvals
    type(fullvalues), intent(inout) :: fullvals
    real(8), optional, intent(in) :: prefactor
    integer :: np1
    real(8), allocatable :: xlattice(:)
    real(8) :: prefac
    integer :: i,j,l
    l=0
    if (present(prefactor)) then
       prefac=prefactor
    else
       prefac=-0.1d0
    end if
    np1=sysparams%np1
    allocate(xlattice(np1))
    xlattice=sysparams%xlattice
    do i=1,np1
       do j=1,np1
          l=l+1
          fullvals%v(l)=sharedvals%vinteract(l)+&
                        sharedvals%vin(i)+sharedvals%vin(j)+&
                        prefac*(xlattice(i)+xlattice(j))
       end do
       sharedvals%v1(i)=sharedvals%vin(i)+prefac*xlattice(i)
       !sharedvals%vin(i)=sharedvals%vin(i)+prefac*xlattice(i)
    end do
  end subroutine add_driving_potential
  
  !adds the driving potential for the 1D 2-electron example
  subroutine calculate_driving_pot(sysparams,prefactor,potder)
    use derivedtypes 
    type(systemparameters), intent(in) :: sysparams
    real(8), intent(in) :: prefactor
    real(8), intent(inout) :: potder(:)
    integer :: np1
    real(8), allocatable :: xlattice(:)
    real(8) :: prefac
    integer :: i,j,l
    l=0
    prefac=prefactor
    np1=sysparams%np1
    allocate(xlattice(np1))
    xlattice=sysparams%xlattice
    do i=1,np1
       do j=1,np1
          l=l+1
          potder(l)=prefac*(xlattice(i)+xlattice(j))
       end do
    end do
  end subroutine calculate_driving_pot 
end Module potential
