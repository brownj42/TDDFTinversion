Module initial_states
  implicit none
contains

  !Generates the full system wavefunction for the examples in the !paper.
  subroutine initializefullsystem(sysparams,fullvals)
    use derivedtypes
    use hamiltonian_mod
    type(systemparameters), intent(in) :: sysparams
    type(fullvalues), intent(inout) :: fullvals
    complex(8), allocatable :: eigf(:,:)
    real(8), allocatable :: eigval(:)
    integer :: np1,ntot,nd,npart
    real(8) :: dx
    npart=sysparams%npart
    nd=sysparams%nd
    np1=sysparams%np1
    ntot=sysparams%ntot
    dx=sysparams%dx

    if (sysparams%quantization==1) then
       print*,'Calculating full Wavefunction starting state'
       if (npart==2.and.nd==1) then
          !calculate eigenvector for full 2 particle system
          allocate(eigf(sysparams%ntot,1))
          allocate(eigval(1))
          call calceigenstates(sysparams,fullvals%v,30,1,eigf,eigval)
          print*,eigval
          fullvals%psi=eigf(:,1)
          deallocate(eigval,eigf)
          fullvals%psi=fullvals%psi/&
               dsqrt(dble(dot_product(fullvals%psi,fullvals%psi)))
       elseif (nd==1.and.npart==1) then
          !fullvals%psi=initial_state_1d(np1,ntot,sysparams%xlattice)
          call initial_state_1d(np1,ntot,sysparams%xlattice,fullvals%psi)
       elseif (nd==3.and.npart==1) then
          !fullvals%psi=initial_state_3d(np1,ntot,dx,sysparams%xlattice)
          call initial_state_3d(np1,ntot,dx,sysparams%xlattice,fullvals%psi)
          
       else
          print*,'no potential created for this choice of parameters'
          print*,'stopping'
          stop
       end if
    else
       fullvals%psi=dcmplx(0.d0,0.d0)
       fullvals%psi(1)=dcmplx(1.d0,0.d0)!/sqrt(2.d0)
    end if
  end subroutine initializefullsystem
  
  !returns one electron ground state in the first KS orbital
  subroutine oneelectrongroundstate(sysparams,sharedvals,KSvals)
    use derivedtypes
    use hamiltonian_mod
    type(systemparameters), intent(in) :: sysparams
    type(sharedvalues), intent(in) :: sharedvals
    type(KSvalues), intent(out) :: KSvals
    type(systemparameters) :: sysparams1

    sysparams1=sysparams
    sysparams1%npart=1
    sysparams1%ntot=sysparams%ntot1
    sysparams1%ntot2=sysparams%ntot1
    fas=0
    evens=.false.
    odds=.false.
    use2=.false.
    call init_ksvalues(sysparams,ksvals) 
    call calceigenstates(sysparams1,sharedvals%v1,25,1,KSvals%phi(:,1:1))
    KSvals%phi(:,2)=0.d0
  end subroutine oneelectrongroundstate 
 
  !Initializes the KSvals%phi orbitals given the parameters in
  !sysparams and sharedvals to match the full electron density dpe
  !fullvals, optional used for case when only one electron by
  !setting KSvals%phi(:,1)=fullvals%psi(:)
  subroutine initializeKSsystem(sysparams,sharedvals,dpe,fullvals,KSvals,occupy,cutoff)
    use derivedtypes 
    use density_mod
    type(systemparameters), intent(in) :: sysparams
    real(8), intent(in) :: dpe(:)
    type(sharedvalues), intent(in) :: sharedvals
    type(fullvalues), optional, intent(in) :: fullvals
    real(8), optional, intent(in) :: occupy(:,:)
    real(8), optional, intent(in) :: cutoff
    type(KSvalues), intent(out) :: KSvals
    type(systemparameters) :: sysparams1part
    real(8), allocatable :: fulldens(:,:)
    integer :: np1,npart,ntot1
    real(8), allocatable :: dp(:)
    complex(8), allocatable :: phi(:,:)
    real(8), allocatable :: U(:,:)
    integer :: i
    np1=sysparams%np1
    npart=sysparams%npart
    ntot1=sysparams%ntot1
    allocate(sysparams1part%T(np1,np1),sysparams1part%xlattice(np1))
    sysparams1part=sysparams
    sysparams1part%npart=1
    sysparams1part%ntot=ntot1
    sysparams1part%ntot2=ntot1
   
    call init_ksvalues(sysparams,ksvals) 
    allocate(phi(ntot1,npart),dp(ntot1))
    print*,'Calculating initial corresponding KS state'
    if (sysparams%quantization==1) then
       if (npart==2) then
          !Generate full one body reduced density matrix
          !allocate(fulldens(ntot1,ntot1))
          !fulldens=fullwf_densitymatrix(ntot1,npart,fullvals%psi)
          !call svd_initial(ntot1,npart,fulldens,phi)
          !deallocate(fulldens)
          if (present(occupy)) then
             call groundstate_initial(sysparams1part,dpe,sharedvals%v1,phi,KSvals%vks,KSvals%init,occupy=occupy)
          else
             if (present(cutoff)) then
                call groundstate_initial(sysparams1part,dpe,sharedvals%v1,phi,KSVals%vks,KSvals%init,cutoff=cutoff)
             else 
                call groundstate_initial(sysparams1part,dpe,sharedvals%v1,phi,KSVals%vks,KSvals%init)
             end if
          end if
       
       
          !output for debugging purposes 
          open(unit=232,file='wavebegin.dat')
          do i=1,ntot1
             write(232,'(3f20.10)')sysparams%xlattice(i),real(phi(i,1)),real(phi(i,2))
          end do
          close(232)
       
          !generate initial KS density
          dp=ks_density(ntot1,npart,phi)
       
          print*,'Properties of initial KS state'
          print*,'density error',sqrt(dot_product(dp-dpe,dp-dpe))
          print*,'particle one norm',dot_product(phi(:,1),phi(:,1))
          print*,'particle two norm',dot_product(phi(:,2),phi(:,2))
          print*,'overlap of one/two',dot_product(phi(:,1),phi(:,2))
          KSvals%phi=phi
       else
          if (present(fullvals)) then
             KSvals%phi(:,1)=fullvals%psi(:)
          else
             print*,'Warning, should use the same psi as the full system for one particle'
             KSvals%phi(:,1)=dpe
          end if
       end if
    else
      do i=1,sysparams%npart
         KSvals%phi(:,i)=dcmplx(0.d0,0.d0)
         KSvals%phi(i,i)=dcmplx(0.d0,-1.d0)**(i+1)
      end do
      allocate(U(sysparams%npart,sysparams%npart))
      call generateunitary(sysparams%npart,U)
      KSvals%phi(1:sysparams%npart,:)=matmul(U,KSvals%phi(1:sysparams%npart,:))
      deallocate(U)
    end if
    dp=ks_density(ntot1,npart,KSvals%phi)
    KSvals%dp=dp

  end subroutine initializeKSsystem

  subroutine generateunitary(npart,U)
     integer, intent(in) :: npart
     real(8), intent(out) :: U(npart,npart)
     complex(8) :: Uc(npart,npart),mid(npart,npart)
     real(8) :: w(npart)
     integer :: i
     integer :: lwork,lrwork,liwork
     real(8) :: qrwork(1)
     integer :: qiwork(1),info
     complex(8) :: qwork(1)
     integer, allocatable :: iwork(:)
     complex(8), allocatable  :: work(:)
     real(8), allocatable :: rwork(:)
     Uc=0.d0
     call random_number(w)
     do i =2,npart
        Uc(2:npart,1)=dcmplx(0.d0,-w(i))
        Uc(1,2:npart)=-Uc(2:npart,1)
     end do
     
     
     lwork=-1
     lrwork=-1
     liwork=-1
     call zheevd('V','U',npart,Uc,npart,W,qwork,lwork,qrwork,lrwork,qiwork,liwork,info)
     lwork=nint(dble(qwork(1)))
     allocate(work(lwork))
     lrwork=nint(qrwork(1))
     allocate(rwork(lrwork))
     liwork=qiwork(1)
     allocate(iwork(liwork))
     call zheevd('V','U',npart,Uc,npart,W,work,lwork,rwork,lrwork,iwork,liwork,info)
     deallocate(work,iwork,rwork)
     mid=0.d0
     do i=1,npart
        mid(i,i)=zexp(dcmplx(0.d0,-w(i)))
     end do
     
     U=dble(matmul(Uc,matmul(mid,dconjg(transpose(Uc)))))
  end subroutine generateunitary
  
  subroutine initial_state_1d(np1,ntot1,xlattice,yh)
    integer, intent(in) :: np1,ntot1
    real(8), intent(in) :: xlattice(np1)
    complex(8), intent(out) :: yh(ntot1)
    real(8) :: pi
    integer :: i
    pi=4.d0*datan(1.d0)
    do i=1,ntot1
       yh(i)=(exp(-(xlattice(i)-0.d0)**2/2.d0))*cos(dcmplx(1.5d0,0.d0)*1.d0*sqrt(pi)*(xlattice(i)-0.0d0))!&
    end do
    yh=yh/dsqrt(dble(dot_product(yh,yh)))
  end subroutine initial_state_1d
  
  !generate initial state given the full density matrix
  subroutine svd_initial(Ntot1,npart,fulldens,ys)
    integer, intent(in) :: Ntot1,npart
    real(8), intent(in) :: fulldens(ntot1,ntot1)
    complex(8), intent(out) :: ys(ntot1,npart)
    integer :: i,j
    integer :: lda,ldu,ldvt
    integer :: lwork,info
    character(1) :: Jobu,Jobvt
    real(8) :: qwork(1)
    real(8), allocatable :: s(:),A(:,:),u(:,:),vt(:,:),work(:)
    real(8), allocatable :: dpe(:)

    !Calculate SVD of full one body reduced density matrix 
    !The two vectors with the large coefficient are a good basis for for the starting vectors
    lda=ntot1
    ldu=ntot1
    ldvt=ntot1
    allocate(s(min(ntot1,ntot1)))
    allocate(A(lda,ntot1))
    allocate(u(ldu,ntot1))
    allocate(vt(ldvt,ntot1))
    Jobu='A'
    Jobvt='A'
    
    A=fulldens
    lwork=-1
    call dgesvd(Jobu,jobvt,ntot1,ntot1,A,lda,s,u,ldu,vt,ldvt,qwork,lwork,info)
    lwork=nint(qwork(1))
    allocate(work(lwork))
    call dgesvd(Jobu,jobvt,ntot1,ntot1,A,lda,s,u,ldu,vt,ldvt,work,lwork,info)
    deallocate(work) 
    print*,'Coefficients for SVD of density matrix'
    print*,s(1:8)
    
    allocate(dpe(ntot1))
    do i=1,ntot1
       dpe(i)=fulldens(i,i)
    end do
    !The starting states need to be played around with to get smooth starting states
    !first state can be the largest svd vector or some combination of the two largest
    ys(:,1)=dcmplx(u(:,1)+u(:,2),0.d0)
    ys(:,2)=dcmplx(u(:,2)-u(:,1),0.d0)
    
    !let the remainder be the second target state, in order to generate smooth/continuous initial states
    !one needs to assign the sign properly
    deallocate(s,a,vt)
    do j=1,100000
       !orthogonalize
       if ((j)/2*2==j) then
          ys(:,2)=ys(:,2)/norm(ys(:,2))
          ys(:,1)=ys(:,1)-dble(dot_product(ys(:,2),ys(:,1))*ys(:,2))
          ys(:,1)=ys(:,1)/norm(ys(:,1))
          u(:,1)=dble(ys(:,1))
          u(:,2)=dble(ys(:,2))
       else
          ys(:,1)=ys(:,1)/norm(ys(:,1))
          ys(:,2)=ys(:,2)-dble(dot_product(ys(:,2),ys(:,1))*ys(:,1))
          ys(:,2)=ys(:,2)/norm(ys(:,2))
          u(:,1)=dble(ys(:,1))
          u(:,2)=dble(ys(:,2))
       end if
       !adjust vectors to obtain correct density
       do i=1,ntot1
          if (dabs(dble(dconjg(ys(i,1))*ys(i,1))).le.dpe(i)) then
             ys(i,2)=dsign(dsqrt(dabs(max(0.d0,dpe(i)-dble(dconjg(ys(i,1))*ys(i,1))))),u(i,2))
          end if
          if (dabs(dble(dconjg(ys(i,2))*ys(i,2))).le.dpe(i)) then 
             ys(i,1)=dsign(dsqrt(dabs(max(0.d0,dpe(i)-dble(dconjg(ys(i,2))*ys(i,2))))),u(i,1))
          end if
       end do
    end do

  end subroutine svd_initial

  subroutine svd(Ntot1,fulldens,singvals)
    integer, intent(in) :: Ntot1
    real(8), intent(in) :: fulldens(ntot1,ntot1)
    real(8), intent(out) :: singvals(ntot1)
    integer :: lda,ldu,ldvt
    integer :: lwork,info
    character(1) :: Jobu,Jobvt
    real(8) :: qwork(1)
    real(8), allocatable :: s(:),A(:,:),u(:,:),vt(:,:),work(:)

    !Calculate SVD of full one body reduced density matrix 
    !The two vectors with the large coefficient are a good basis for for the starting vectors
    lda=ntot1
    ldu=ntot1
    ldvt=ntot1
    allocate(s(min(ntot1,ntot1)))
    allocate(A(lda,ntot1))
    allocate(u(ldu,ntot1))
    allocate(vt(ldvt,ntot1))
    Jobu='A'
    Jobvt='A'
    
    A=fulldens
    lwork=-1
    call dgesvd(Jobu,jobvt,ntot1,ntot1,A,lda,s,u,ldu,vt,ldvt,qwork,lwork,info)
    lwork=nint(qwork(1))
    allocate(work(lwork))
    call dgesvd(Jobu,jobvt,ntot1,ntot1,A,lda,s,u,ldu,vt,ldvt,work,lwork,info)
    deallocate(work) 
    print*,'Coefficients for SVD of density matrix'
    print*,s(1:5)
    singvals=s

   

  end subroutine svd

  subroutine groundstate_initial(sysparams,dpe,v1,ys,vksout,success,occupy,cutoff)
    use derivedtypes 
    use hamiltonian_mod
    type(systemparameters), intent(in) :: sysparams
    real(8), intent(in) :: dpe(:),v1(:)
    complex(8), intent(out) :: ys(:,:)
    real(8), intent(out) :: vksout(:)
    complex(8), allocatable :: wave(:,:)
    integer, intent(out) :: success
    real(8), optional, intent(in) :: occupy(:,:)
    real(8), optional, intent(in) :: cutoff
    integer :: ntot1
    real(8) :: ener,errorold,errornew,prefac,ecut
    real(8), allocatable :: T(:,:)
    integer :: i,j
    integer :: postrac(3,2),loc(1)
    real(8) :: occ(3,2)
    real(8), allocatable :: eighold(:,:,:)
    real(8), allocatable :: vks(:),dp(:),change(:)
    real(8), allocatable :: eigenstates(:,:)

    if (present(cutoff)) then
       ecut=cutoff
    else 
       ecut=5.d-13
    end if    

    !postrac(2,2) !number of states for each, number of particles
    postrac(1,1)=1
    postrac(2,1)=3
    postrac(3,1)=5
    postrac(1,2)=2
    postrac(2,2)=4
    postrac(3,2)=6
    success=0
    occ=0.d0
    if (sysparams%occupy_case==0) then
       occ(1,:)=1.d0
    elseif (sysparams%occupy_case==1) then
       occ(2,1)=1.d0
       occ(1,2)=1.d0
    elseif (sysparams%occupy_case==2) then
       occ(1,1)=1.d0
       occ(1,2)=1.0d0
       occ(2,2)=0.5d0
    elseif (sysparams%occupy_case==3) then
       occ(1,1)=3.5d0
       occ(2,1)=0.4d0
       occ(2,1)=0.1d0
       occ(1,2)=5.5d0
       occ(2,2)=0.0d0
       occ(3,2)=0.0d0
    elseif (sysparams%occupy_case==4) then
       occ(1,1)=1.0d0
       occ(2,1)=0.3d0
       occ(3,1)=0.0d0
       occ(1,2)=5.5d0
       occ(2,2)=1.0d0
       occ(3,2)=0.0d0
    elseif (sysparams%occupy_case==5) then
       occ(1,1)=5.d0
       occ(2,1)=2.0d0
       occ(3,1)=0.1d0
       occ(1,2)=5.d0
       occ(2,2)=2.0d0
       occ(3,2)=0.1d0
    elseif (sysparams%occupy_case==6) then
       occ(1,1)=4.0d0
       occ(2,1)=1.0d0
       occ(3,1)=0.1d0
       occ(1,2)=5.d0
       occ(2,2)=1.1d0
       occ(3,2)=1.1d0
    else
       occ(1,1)=1.5d0
       occ(2,1)=1.0d0
       occ(3,1)=0.5d0
       occ(1,2)=1.5d0
       occ(2,2)=1.0d0
       occ(3,2)=0.5d0
    end if
    if (present(occupy)) then
       occ(1:3,1:2)=occupy(1:3,1:2)
       print*,occ(1,:)
       print*,occ(2,:)
       print*,occ(3,:)

    end if
    
    errorold=2.d0

    
    ntot1=sysparams%ntot1
    
    allocate(eighold(ntot1,3,2))
    if (sysparams%nd==1) then
       allocate(eigenstates(ntot1,ntot1))
    else
       allocate(wave(ntot1,3))
    end if
    
    allocate(vks(ntot1),dp(ntot1))
    allocate(T(ntot1,ntot1))
    allocate(change(ntot1))
    T=sysparams%T
    !This section is if you want a completely stationary system
    vks =v1
    call random_number(vks)
    vks=vks*0.000d0+v1+10.d0
    prefac=0.05
    do i=1,20000
       evens=.true.
       odds=.false.
       call forcesym1(vks,ntot1)
       fas=1
       evens=.true.
       odds=.false.
       use2=.false.
       if (sysparams%nd>1) then
          call calceigenstates(sysparams,vks,25,3,wave)
          ys(:,1)=wave(:,1)
       else
          call calcalleigenstates(ntot1,T,vks,eigenstates)
          if (i.gt.1) then
          do j=1,3
             loc=maxloc(abs(matmul(eighold(:,j,1),eigenstates)))
             !print*,j,abs(matmul(eighold(:,j,1),eigenstates))
             postrac(j,1)=loc(1)
             loc=maxloc(abs(matmul(eighold(:,j,2),eigenstates)))
             postrac(j,2)=loc(1)
          end do
          end if
          !print*,postrac
          !pause 
          ys=0.d0
          do j=1,3
             evens=.true.
             odds=.false. 
             call forcesym1(eigenstates(:,postrac(j,1)),ntot1)
             evens=.false.
             odds=.true. 
             call forcesym1(eigenstates(:,postrac(j,2)),ntot1)
             eigenstates(:,postrac(j,2))=eigenstates(:,postrac(j,2))/&
                       sqrt(dot_product(eigenstates(:,postrac(j,2)),eigenstates(:,postrac(j,2))))
             eigenstates(:,postrac(j,1))=eigenstates(:,postrac(j,1))/&
                       sqrt(dot_product(eigenstates(:,postrac(j,1)),eigenstates(:,postrac(j,1))))
             eighold(:,j,1)=sign(1.d0,dot_product(eighold(:,j,1),eigenstates(:,postrac(j,1))))*eigenstates(:,postrac(j,1))
             eighold(:,j,2)=sign(1.d0,dot_product(eighold(:,j,2),eigenstates(:,postrac(j,2))))*eigenstates(:,postrac(j,2))
             ys(:,1)=ys(:,1)+occ(j,1)*eighold(:,j,1)
             ys(:,2)=ys(:,2)+occ(j,2)*eighold(:,j,2)
          end do
          if (sysparams%occupy_case==0.and.sysparams%singlet==1) then
             ys(:,2)=ys(:,1)
          end if
          do j=1,2
             ys(:,j)=ys(:,j)/sqrt(abs(dot_product(ys(:,j),ys(:,j))))
          end do
            
          !print*,matmul(transpose(ys),ys)
       end if
       fas=1
       evens=.false.
       odds=.true.
       use2=.false.
       if (sysparams%occupy_case==1) then
          use2=.true.
       end if
       
       if (sysparams%nd>1) then
          call calceigenstates(sysparams,vks,25,3,wave)
          ys(:,2)=wave(:,1)
       end if
       dp =0.d0
       ener=0.d0
       do j=1,2
          ener=ener+dble(dot_product(matmul(T,ys(:,j))+vks*ys(:,j),ys(:,j)))
          dp=dp+abs(dconjg(ys(:,j))*ys(:,j))
       end do
       change=-prefac*dlog(dp/dpe)
       change=-prefac*dlog(1.d0+(dp-dpe)/(dp+1.1d0*dpe))
       !print*,sign(min(abs(change),0.01d0),change)
       vks=vks-change!sign(min(abs(change),40.1d0),change)
       errornew=sqrt(dot_product(dpe-dp,dpe-dp))
       if ((i/100)*100==i) then
       print*,'iteration',i,'with error',errornew
       end if
          if (errornew.lt.errorold.and.(i/1000)*1000==i) then
          prefac=prefac*1.1d0
          end if
       if (errornew.lt.ecut) then
          success=1
          exit
       end if
       if (abs(errornew-errorold)/errornew.lt.1.d-14.and.errornew.gt.1.d-6) then
          call random_number(vks)
          vks=vks+v1
          prefac=0.1d0
       end if
       if (errorold.lt.errornew) then
          !stop
          if (prefac.gt.1.d-2) then
          prefac=prefac/1.1d0
          end if
          errorold=errornew
       else
          errorold=errornew
       end if
       vks=vks-(ener+15.d0)/2.d0
    end do
    vksout=vks
    if (sysparams%nd==1) then
       deallocate(eigenstates)
    end if
  end subroutine groundstate_initial

  !calculates all eigenstates of a 1D Hamiltonian with KEO matrix T and
  !diagonal potential v of length np1
  subroutine calcalleigenstates(np1,T,v,eigenstates)
     integer, intent(in) :: np1
     real(8), intent(in) :: T(:,:)
     real(8), intent(in) :: v(:)
     real(8), intent(out) :: eigenstates(np1,np1)
     character(1) :: JOBZ,UPLO
     real(8) :: eigenvals(np1)
     integer :: LDA
     real(8), allocatable :: work(:)
     integer :: liwork
     integer, allocatable :: iwork(:)
     real(8) :: qwork(1)
     integer :: qiwork(1)
     integer :: lwork,info
     integer :: i,j
       
     JOBZ='V'
     UPLO='U'
     LDA=np1
     do i=1,np1
        do j=1,np1
           eigenstates(i,j)=T(i,j)
        end do
        eigenstates(i,i)=eigenstates(i,i)+v(i)
     end do
     lwork=-1
     qiwork=-1     
     !call DSYEV( JOBZ, UPLO, Np1, eigenstates, LDA, eigenvals, qWORK, LWORK, INFO )
     call DSYEVD( JOBZ, UPLO, Np1, eigenstates, LDA, eigenvals, qWORK, LWORK, qiwork,liwork,INFO )
     lwork=nint(qwork(1))
     allocate(work(lwork))
     liwork=(qiwork(1))
     allocate(iwork(liwork))
     !call DSYEV( JOBZ, UPLO, Np1, eigenstates, LDA, eigenvals, WORK, LWORK, INFO )
     call DSYEVD( JOBZ, UPLO, Np1, eigenstates, LDA, eigenvals, WORK, LWORK, iwork,liwork,INFO )
     deallocate(work)
     deallocate(iwork)
     
  end subroutine calcalleigenstates
 
  
  !calculates L2-norm of a vector
  real(8) function norm(v)
    complex(8) :: v(:)
    norm=dsqrt(dabs(dble(dot_product(v,v))))
  end function norm

  subroutine initial_state_3d(np1,ntot,dx,xlattice,yh)
    integer, intent(in) :: np1,ntot
    real(8), intent(in) :: dx,xlattice(np1)
    complex(8), intent(out) :: yh(ntot)
    integer :: i,j,k,l
    real(8) :: pi
    pi=4.d0*datan(1.d0)

    l=0
    do i=1,np1
       do j=1,np1
          do k=1,np1
             l=l+1
             yh(l)=(dsqrt(dx)/dsqrt(dsqrt(pi)))**3*&
                  zexp(dcmplx(0.d0,-0.25d0*dsqrt(pi)*xlattice(i)))*&
                  zexp(dcmplx(0.d0,0.25d0)*dsqrt(pi)*xlattice(j))*&
                  dexp(-((xlattice(i)-0.0d0)**2+&
                  (xlattice(j)+0.0d0)**2+&
                  (xlattice(k)+0.25d0*sqrt(pi))**2)/2.d0)
          end do
       end do
    end do
    yh=yh/norm(yh)
  end subroutine initial_state_3d

  
  !calculates eigenstates of a system defined by sysparams with potential vks
  subroutine calceigenstates(sysparams,vks,niter,nev,yout,eigenvalues)
    use derivedtypes
    use sortmod
    use matmul_mod
    use hamiltonian_mod
    !np1: number of points in 1D grid points
    !npc: number of particles
    !nd: number of dimensions
    !niter: number of iterations, block size for either
    !       short iterative arnoldi or implicitly restarted arnoldi
    !ntsteps: If ntsteps=0, time-independent energy calculation
    !T:    Matrix of 1D 2nd derivative np1 x np1 matrix
    !vks:  Vector of potential at each grid point last dimension fast
    !      Lexigraphical order
    !dx    Grid spacing
    !dtin  Time step
    !yin   Starting vector before time step
    !neigs How many eigenvalues to calculate if ntsteps=0
    type(systemparameters), intent(in) :: sysparams
    integer, intent(in) :: niter
    real(8), intent(in) :: Vks(:)!,Vks(np1**npc**nd)
    integer, intent(in):: nev
    complex(8), optional, intent(out) :: yout(:,:)
    real(8), optional, intent(out) :: eigenvalues(nev)
    integer Np1,ntot,npc,nd
    real(8), allocatable :: T(:,:)!T(np1,np1)
    integer :: neigs
    real(8), allocatable :: veigs(:)
    real(8), allocatable :: wffunc(:,:)
    complex(8), allocatable :: vec(:),vecp(:),vecpp(:)
    integer :: i,j,k,ip
    integer :: s,c
    real(8), allocatable :: valwf(:,:,:)
    integer :: info,iter
    integer, allocatable :: ind(:,:)
    integer, allocatable :: pind(:,:,:)
    real(8) :: tol
    integer :: ido,ds,df,dsp,dfp,iv
    integer, allocatable :: iparam(:),ipntr(:)
    real(8), allocatable :: resid(:),workd(:),workl(:)
    real(8), allocatable :: qs(:,:)
    integer :: lworkl
    logical :: rvec
    character(1) :: howmny
    character(12) :: filename
    logical, allocatable :: select(:)
    real(8) :: sigmar,sigmai
    real(8), allocatable :: DR(:),DI(:),Zeig(:,:),workev(:)
    include 'mkl_blas.fi' 
    Np1=sysparams%np1
    ntot=sysparams%ntot
    npc=sysparams%npart
    nd=sysparams%nd
    allocate(T(np1,np1))
    T=sysparams%T
   
    if (sysparams%singlet==1) then
       singlet=.TRUE.
    else
       singlet=.False.
    end if
    if (sysparams%triplet==1) then
       triplet=.TRUE.
    else
       triplet=.False.
    end if 
    neigs=nev
    allocate(ind(ntot,3))
    allocate(pind(np1,np1,np1))
    !index vector, only works for 1 or two dimensional grid right now
    ip=0
    if (npc==2) then
       do i=1,np1
          do j=1,np1
             ip=ip+1
             ind(ip,1)=i
             ind(ip,2)=j
             iv=(np1+1)/2
             pind(i,j,iv)=ip
          end do
       end do
    elseif(npc==1.and.nd==1) then
       do i=1,np1
          ip=ip+1
          ind(ip,1)=i
          iv=(np1+1)/2
          pind(i,iv,iv)=ip
       end do
    else
       do i=1,np1
          do j=1,np1
             do k=1,np1
                ip=ip+1
                
                ind(ip,1)=i
                ind(ip,2)=j
                ind(ip,3)=k
                pind(i,j,k)=ip
             end do
          end do
       end do
    end if



    allocate(vecpp(ntot),vecp(ntot),vec(ntot))
    allocate(wffunc(ntot,Neigs))
    allocate(veigs(neigs+1))
    !force asymmetry under exchange of wavefunction
    !if two electrons
    if (npc==2) then
       call forcesymc(vecpp,Np1,Ntot)
    end if
    roswitch=0
    allocate(reord(ntot,3))
    
    

    
    iter=niter
    allocate(qs(ntot,niter))
    call random_number(qs(:,1))
    vecp=qs(:,1)
    vecp=vecp/dsqrt(dble(zdotc(ntot,vecp,1,vecp,1)))
    qs(:,1)=dble(vecp)
    if (npc==2) then
       call forcesym(qs(:,1),Np1,Ntot)
    end if
    if (npc==1) then
       call forcesym1(qs(:,1),Np1)
    end if
    
    !vecp=qs(:,1)
    !vecp=vecp/sqrt(zdotc(ntot,vecp,1,vecp,1))
    
    qs(:,1)=qs(:,1)/dsqrt(ddot(ntot,qs(:,1),1,qs(:,1),1))
    !ARPACK input parameters
    allocate(iparam(11))
    allocate(ipntr(14))
    lworkl=3*iter**2+6*iter
    allocate(workl(lworkl))
    allocate(resid(ntot))
    allocate(workd(3*ntot))
    iparam(1)=1
    iparam(3)=2*iter
    iparam(4)=1
    iparam(7)=1
    ipntr(1)=1
    tol=1.d-14!epsilon(1.d0)
    ido=0
    info=1
    resid=qs(:,1)
    workd(1:ntot)=qs(:,1)
    workd(ntot+1:2*ntot)=qs(:,1)
    workd(2*ntot+1:3*ntot)=qs(:,1)

    !print*,'step',tsteps
    arpack:do
       !call dnaupd(ido,'I',ntot,'SR',Neigs,tol,resid,niter,&
       !     qs(:,1:niter),ntot,iparam,ipntr,workd,workl,lworkl,info)
       call dsaupd(ido,'I',ntot,'SA',Neigs,tol,resid,niter,&
            qs(:,1:niter),ntot,iparam,ipntr,workd,workl,lworkl,info)

       !hamiltonian times vector section
       if (ido.eq.-1.or.ido.eq.1) then
          
          ds=ipntr(1)
          df=ipntr(1)+ntot-1
          if (npc==2) then
             call forcesym(workd(ds:df),Np1,Ntot)
          end if
          if (npc==1) then
             !call forcesym1(workd(ds:df),Np1)
          end if
          dsp=ipntr(2)
          dfp=ipntr(2)+ntot-1
          vecp=workd(ds:df)
       
          !Apply hamiltonian to vecp
          !vec=H*vecp
          if (npc==1) then
             if (nd==1) then
                vec=matmul(T,vecp)+vks*vecp
             elseif (nd==2) then
                vec=matmulac2(T,vecp,ntot)+vks*vecp
             elseif (nd==3) then
                !vec=matmulac(T,vecp,ntot)+vks*vecp
                vec=matmulac(T,vecp,ntot)
                call zaxpy(ntot,dcmplx(1.d0,0.d0),vtimeszv(vks,vecp,ntot),1,vec,1)
             else
                print*,'code not written for these parameters'
                print*,'stopping'
                stop
             end if
          elseif (npc==2) then
             if (nd==1) then
                vec=matmulac2(T,vecp,ntot) +vks*vecp
             else
                print*,'code not written for these parameters'
                print*,'stopping'
                stop
             end if
          elseif (npc==3) then
             if (nd==1) then
                vec=matmulac(T,vecp,ntot)+vks*vecp
             else
                print*,'code not written for these parameters'
                print*,'stopping'
                stop
             end if
          end if


          !ensure appropriate symmetry remains
          
          if (npc==2) then
             workd(dsp:dfp)=dble(vec)
             call forcesym(workd(dsp:dfp),Np1,Ntot)
          else
             workd(dsp:dfp)=dble(vec)
             if (npc==1) then
                call forcesym1(workd(dsp:dfp),Np1)
             end if
             workd(dsp:dfp)=workd(dsp:dfp)-2.d0*workd(ds:df)
          end if
       else
          !print*,'ido',ido
          exit arpack
       end if
    end do arpack

    

    !calculate time-independant eigenvalues/eigenvectors
    rvec=.True.
    howmny='A'
    !allocate(DR(Neigs+1),DI(Neigs+1),Zeig(ntot,Neigs+1),workev(3*iter))
    allocate(DR(Neigs),DI(Neigs),Zeig(ntot,Neigs),workev(3*iter))
    allocate(select(iter))
    select=.TRUE.
    sigmar=0.d0
    sigmai=0.d0
    zeig=qs(:,1:Neigs+1)
    !call dneupd2(rvec,howmny,select,dr,di,zeig,ntot,sigmar,sigmai,workev,&
    !     'I',ntot,'SR',Neigs,tol,resid,iter,qs(:,1:iter),ntot,iparam,&
    !     ipntr,workd,workl,lworkl,info)
    call dseupd(rvec,howmny,select,dr,zeig,ntot,sigmar,&
         'I',ntot,'SA',Neigs,tol,resid,iter,qs(:,1:iter),ntot,iparam,&
         ipntr,workd,workl,lworkl,info)
    if (npc==1) then
       veigs(1:neigs)=(dr(1:neigs)+2.d0)
       if (present(eigenvalues)) then
          eigenvalues=dr(1:nev)+2.d0
       end if
    else
       veigs=dr
       if (present(eigenvalues)) then
          eigenvalues=dr(1:nev)
       end if
       !print*,(di(1:Neigs))!*219474.63
    end if
    
    wffunc=zeig(:,1:neigs) 
    deallocate(DR,DI,Zeig,workev)
    deallocate(select)
    deallocate(qs)



    if (present(yout)) then
       !transfer eigenvectors to yout
       yout(:,1:nev)=wffunc(:,1:nev)
    end if
       
    
    
    deallocate(iparam)
    deallocate(ipntr)
    deallocate(workl)
    deallocate(resid)
    deallocate(workd)
    


    deallocate(reord)
    deallocate(vecpp,vecp,vec)
    deallocate(wffunc)
    deallocate(veigs)
    deallocate(ind)
    deallocate(pind)

  end subroutine calceigenstates


  
end Module initial_states
