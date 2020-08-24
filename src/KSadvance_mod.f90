module KSadvance_mod
  implicit none
contains
  subroutine advanceKSsystem(dpe,dpenew,dnx,ddnx,ddnxnew,sysparams,KSvals,sharedvals,info)
    use derivedtypes 
    use density_mod
    use hamiltonian_mod
    use propagate 
    use pot_invert_mod
    use assign_phases_mod
    use outputdata_mod
    real(8), intent(in) :: dnx(:),dpe(:),dpenew(:),ddnx(:),ddnxnew(:)
    type(systemparameters), intent(inout) :: sysparams
    type(KSvalues), intent(inout) :: KSvals
    type(sharedvalues), intent(in) :: sharedvals
    integer, intent(out) :: info
    type(systemparameters) :: sysparams1part
    integer :: np1,nd,ntot1,npart,ntot2
    integer :: niter,sw,swh
    integer :: loop
    complex(8), allocatable :: phi(:,:),phihold(:,:),phinew(:,:)
    complex(8), allocatable :: phinew2(:,:)
    complex(8), allocatable :: phiplus(:,:),phiminus(:)
    complex(8), allocatable :: f(:)
    real(8), allocatable :: vks(:),vksh(:),vkshh(:),dvks(:),dvksh(:)
    real(8), allocatable :: T(:,:),vin(:),dp(:),xlattice(:)
    real(8), allocatable :: v1(:),vhar(:)
    real(8) :: diff,density_error,dt,dth
    real(8) :: ct,enerh
    real(8) :: dvksmax
    integer :: i,j
    character(20) :: outform 
    complex(8), parameter :: j1=dcmplx(0.d0,1.d0)
    integer :: docalcvks,vksconverged
    integer :: pinv0minresqlp1
    niter=5
    !copy system parameters to local variables for ease of reading
    np1=sysparams%np1
    nd=sysparams%nd
    ntot1=sysparams%ntot1
    ntot2=sysparams%ntot2
    npart=sysparams%npart
    sw=sysparams%dvks_sw
    swh=sysparams%dvksh_sw
    ct=sysparams%ct
    dt=sysparams%dt
    dth=sysparams%dth
    loop=sysparams%loop
    dvksmax=sysparams%dvksmax
    pinv0minresqlp1=sysparams%pinv0minresqlp1
    if (pinv0minresqlp1==0.and.nd.ne.1) then
       print*, 'Can only use pseudoinverse in 1D'
       print*, 'stopping'
       stop
    end if
    allocate(T(np1,np1),xlattice(np1))
    allocate(vin(ntot1),v1(ntot1))
    call dcopy(ntot1,sharedvals%vin,1,vin,1)
    call dcopy(ntot1,sharedvals%v1,1,v1,1)
    xlattice=sysparams%xlattice
    T=sysparams%T
    allocate(sysparams1part%T(np1,np1),sysparams1part%xlattice(np1))
    sysparams1part=sysparams
    sysparams1part%npart=1
    sysparams1part%ntot=ntot1

    !copy Kohn-Sham variables to local variables for ease of reading
    allocate(phi(ntot1,npart),dp(ntot1),vhar(ntot1))
    do i=1,npart
        call zcopy(ntot1,KSvals%phi(:,i),1,phi(:,i),1)
    end do
    call dcopy(ntot1,KSvals%dp,1,dp,1)
    allocate(vks(ntot1),vkshh(ntot1),vksh(ntot1))
    !vks=KSvals%vks
    call dcopy(ntot1,KSvals%vks,1,vks,1)
    !vksh=KSvals%vksh
    call dcopy(ntot1,KSvals%vksh,1,vksh,1)
    !vkshh=KSvals%vkshh
    call dcopy(ntot1,KSvals%vkshh,1,vkshh,1)
    allocate(dvksh(ntot1),dvks(ntot1))
    !dvksh=KSvals%dvksh
    call dcopy(ntot1,KSvals%dvks,1,dvks,1)
    !dvks=KSvals%dvks
    call dcopy(ntot1,KSvals%dvksh,1,dvksh,1)
    vhar=KSvals%vhar
    !local variables for time-advance
    allocate(phihold(ntot1,npart))
    allocate(phinew(ntot1,npart),phinew2(ntot1,npart))
    allocate(phiminus(ntot1),phiplus(ntot1,npart))
    allocate(f(ntot1))


    !assign phases
    phihold=assign_phases(KSvals%phi,dnx,T,np1,nd,ntot1,npart,1.d-14)
    !shift current orbitals to phi 
    phi=phihold
    do i=1,npart
        call zcopy(ntot1,phihold(:,i),1,phi(:,i),1)
    end do
    !calculate vks potential
    if (loop.le.2) then
       call calcvks(niter,Np1,nd,ntot1,npart,T,pinv0minresqlp1,ddnx,phi,vks)
       !call calcdvks(np1,nd,ntot1,npart,swh+sw,enerh,dt,T,vin,phi,vkshh,vksh,vks,dvks)
       enerh=sysparams%energy
       call recenter_potential(np1,nd,ntot1,npart,enerh,T,phi,vks)
       dp=ks_density(ntot1,npart,phi)
       KSvals%vhar=calcvhar(ntot1,npart,ntot2,sharedvals%vinteract,dp)
       if (loop==2.and.sysparams%outputdata==1) then
           call printdata(np1,nd,ntot1,npart,xlattice,ct,vks,KSvals%vhar,v1,dp,dpe,phi,sysparams%newfile)
       end if
    else
       enerh=sysparams%energy
       !call calcvks(niter,Np1,nd,ntot1,npart,T,pinv0minresqlp1,ddnx,phi,vks)
       !call recenter_potential(np1,nd,ntot1,npart,enerh,T,phi,vks)
       !call calcdvks(np1,nd,ntot1,npart,swh+sw,enerh,dt,T,vin,phi,dvksmax,vkshh,vksh,vks,dvks)
    end if
    vkshh=vksh
    vksh=vks-vin
    
    !do first part of propagation for each particle
    vks=vks-vin
    do i=1,npart
       f=calcf1(np1,nd,ntot1,dt,T,vin,vks,dvks,phi(:,i))
       phiminus=phi(:,i)-dt*j1/2.d0*vks*phi(:,i)-f
       call advancewf(sysparams1part,sharedvals,niter,vin,phiminus,phiplus(:,i))
    end do

    !first approximation of next step
    do i=1,npart
       phinew(:,i)=(phiplus(:,i)+f)/(1.d0+dt*j1/2.d0*vks)
    end do

    !iterate and recalculate vks potential until converged
    dvksh=dvks
    diff=1.d0 
    docalcvks=1
    vksconverged=0
    swh=sw
    print*,'Starting advance of KS orbitals'
    enerh=sysparams%energynew
    iter:do j=1,100

       !take current approximation and calculate corresponding Vks potential
       if (docalcvks==1.and.diff.lt.1.d-32) then
          print*,'Orbitals converged, recalculate potential'
          call calcvks(niter,Np1,nd,ntot1,npart,T,pinv0minresqlp1,ddnxnew,phinew,vks)
          !recenter potential to keep energy constant
          call recenter_potential(np1,nd,ntot1,npart,enerh,T,phinew,vks)
          docalcvks=0
          
          !can try to approximate derivative of potential
          call calcdvks(np1,nd,ntot1,npart,sw,enerh,dt,T,vin,phinew,dvksmax,vkshh,vksh,vks,dvks)
          vks=vks-vin
       end if
       !dvks=0.d0

       !calculate new approximation to orbitals
       do i=1,npart
          f=calcf1(np1,nd,ntot1,dt,T,vin,vks,dvks,phinew(:,i))
          phinew2(:,i)=(phiplus(:,i)+f)/(1.d0+j1*dt/2.d0*vks)
       end do

       !calculate convergence error 
       diff=0.d0
       do i=1,npart
          diff=diff+dble(dot_product(phinew2(:,i)-phinew(:,i),phinew2(:,i)-phinew(:,i)))
       end do
       write(*,'(5X,A22,1x,es17.10)')'iteration difference =',diff/dble(ntot1)
       diff=diff/dble(ntot1) 
       !transfer new approximation to output KS orbitals
       phinew=phinew2
       if (j.gt.1.and.docalcvks==0.and.diff.lt.1.d-16) then
          if (vksconverged==0) then
             print*,'Potential Converged'
          end if
          docalcvks=0 !potential converged, run extra loops to converge wavefunction
          vksconverged=1
       else
          docalcvks=1 !potential not converged, recalculate potential 
       end if
       

       
       !check for max iterations taken
       if (j==100) then
          print*,'Potential did not converge, halving timestep'
          dt=dt/2.d0
          phinew=phi
          vks=vksh+vin
          dvks=dvksh
          swh=1
          sw=1
          vksconverged=0
          exit iter
       end if

       !if potential and orbitals converged exit
        if (vksconverged==1.and.diff.lt.1.d-32) then
           print*,'Orbitals and potential converged, Time step complete'
           vks=vks+vin
           if (loop.le.1) then
              print*,'restart with dvks saved'
              phinew=phi
              sw=1
              KSvals%dvks=dvks
              KSvals%vks=vks
              KSvals%vksh=vksh
              vksconverged=0
              exit iter
           end if
           
          
           
           if (dt.lt.dth.and.swh==0) then
              ct=ct+dt
              dt=dt*2.d0
              sw=1
              swh=1
           else
              ct=ct+dt
              sw=0
           end if
           dp=ks_density(ntot1,npart,phinew)
           KSvals%vhar=calcvhar(ntot1,npart,ntot2,sharedvals%vinteract,dp)
           if (sysparams%outputdata==1) then
              call printdata(np1,nd,ntot1,npart,xlattice,ct,vks,KSvals%vhar,v1,dp,dpe,phi,sysparams%newfile)
           end if

           
           !add back 1-body time-independant potential
           exit iter
        end if
     end do iter
     

     !calculate new KS density if time step is taken
     if (vksconverged==1) then
        dp=ks_density(ntot1,npart,phinew)
        density_error= dsqrt(dot_product(dp-dpenew,dp-dpenew)) 
        write(*,'(1x,A27,1x,i0,1x,A1,1x,es16.10)')'Density error for time step',loop,'=', density_error
     else
        density_error= dsqrt(dot_product(dp-dpe,dp-dpe)) 
        write(*,'(1x,A27,1x,i0,1x,A1,1x,es16.10)')'Density error for time step',loop,'=', density_error
     end if

     !Is particle number conserved?
     print*,'Check properties of orbitals'
     print*,'Particle number',sum(dp)
     if (npart==2) then
        print*,'Overlap Integral',abs(dot_product(phinew(:,1),phinew(:,2)))
     end if
     write(outform,'(a1,i0,a7)') '(',sysparams%npart,'es10.2)'
     print*,'Overlap matrix'
     do i=1,sysparams%npart
         write(*,outform) abs(matmul(dconjg(transpose(phinew)),phinew(:,i)))
     end do

     if (vksconverged==1) then
        !transfer output states to input for next loop
        KSvals%phi=phinew
        KSvals%vks=vks
        KSvals%vksh=vksh
        KSvals%vkshh=vkshh
        KSvals%vhar=calcvhar(ntot1,npart,ntot2,sharedvals%vinteract,dp)
        KSvals%dvks=dvks
        KSvals%dvksh=dvksh
        KSvals%dp=dp
        info=1
     else
        info=0
     end if
     sysparams%ct=ct
     sysparams%dt=dt
     sysparams%dvks_sw=sw
     sysparams%dvksh_sw=swh
     sysparams%loop=sysparams%loop+1
     

   end subroutine advanceKSsystem

 end module KSadvance_mod
