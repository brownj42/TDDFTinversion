Module initial_states
  implicit none
contains

  !Generates the full system wavefunction for the examples in the !paper.
  subroutine initializefullsystem(sysparams,fullvals)
    use derivedtypes
    use hamiltonian_mod
    type(systemparameters), intent(in) :: sysparams
    type(fullvalues), intent(inout) :: fullvals
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
          call calceigenstates(sysparams,fullvals%v,30,1,fullvals%psi)
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
    call calceigenstates(sysparams1,sharedvals%v1,25,1,KSvals%phi(:,1))
    KSvals%phi(:,2)=0.d0
  end subroutine oneelectrongroundstate 
 
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
    type(fullvalues), optional, intent(in) :: fullvals
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
       
          call groundstate_initial(sysparams1part,dpe,sharedvals%v1,phi)
       
       
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

  subroutine groundstate_initial(sysparams,dpe,v1,ys)
    use derivedtypes 
    use hamiltonian_mod
    type(systemparameters), intent(in) :: sysparams
    real(8), intent(in) :: dpe(:),v1(:)
    complex(8), intent(out) :: ys(:,:)
    integer :: ntot1
    real(8) :: ener,errorold,errornew,prefac
    real(8), allocatable :: T(:,:)
    integer :: i,j
    integer :: postrac(3,2),loc(1)
    real(8) :: occ(3,2)
    real(8), allocatable :: eighold(:,:,:)
    real(8), allocatable :: vks(:),dp(:),change(:)
    real(8), allocatable :: eigenstates(:,:)
    
    !postrac(2,2) !number of states for each, number of particles
    postrac(1,1)=1
    postrac(2,1)=3
    postrac(3,1)=5
    postrac(1,2)=2
    postrac(2,2)=4
    postrac(3,2)=6
    occ=0.d0
    if (sysparams%occupy_case==0) then
       occ(1,:)=1.d0
    elseif (sysparams%occupy_case==1) then
       occ(2,1)=1.d0
       occ(1,2)=1.d0
    elseif (sysparams%occupy_case==2) then
       occ(1,1)=1.d0
       occ(2,1)=0.2d0
       occ(1,2)=1.0d0
    elseif (sysparams%occupy_case==3) then
       occ(1,1)=3.5d0
       occ(2,1)=0.3d0
       occ(2,1)=0.1d0
       occ(1,2)=5.5d0
       occ(2,2)=0.5d0
       occ(3,2)=0.1d0
    elseif (sysparams%occupy_case==4) then
       occ(1,1)=3.0d0
       occ(2,1)=0.7d0
       occ(3,1)=0.2d0
       occ(1,2)=7.0d0
       occ(2,2)=1.7d0
       occ(3,2)=0.2d0
    elseif (sysparams%occupy_case==5) then
       occ(1,1)=5.d0
       occ(2,1)=0.5d0
       occ(3,1)=0.5d0
       occ(1,2)=5.d0
       occ(2,2)=0.5d0
       occ(3,2)=0.5d0
    else
       occ(1,1)=4.d0
       occ(2,1)=4.0d0
       occ(3,1)=4.0d0
       occ(1,2)=9.d0
       occ(2,2)=1.0d0
       occ(3,2)=1.0d0
    end if
    
    errorold=2.d0
   
    ntot1=sysparams%ntot1
    allocate(eighold(ntot1,3,2))
    if (sysparams%nd==1) then
       allocate(eigenstates(ntot1,ntot1))
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
          call calceigenstates(sysparams,vks,25,2,ys(:,1))
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
          call calceigenstates(sysparams,vks,25,2,ys(:,2))
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
       if (errornew.lt.1.d-11) then
          exit
       end if
       if (abs(errornew-errorold)/errornew.lt.1.d-12.and.errornew.gt.1.d-6) then
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
end Module initial_states
