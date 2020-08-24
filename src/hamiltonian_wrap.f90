module propagate 
  implicit none
contains
    
  !subroutine calcynew(Np1,ntot,npc,nd,niter,T,vks,dx,dtin,yin,yout)
  subroutine advancewf(sysparams,sharedvals,niter,vks,yin,yout)
    use derivedtypes 
    use sortmod
    use matmul_mod
    use hamiltonian_mod
    use secondquant_mod
    use mapmod
    use nchoosekmod
    !np1: number of points in 1D grid points
    !npc: number of particles
    !nd: number of dimensions
    !niter: number of iterations, block size for either
    !       short iterative arnoldi or implicitly restarted arnoldi
    !T:    Matrix of 1D 2nd derivative np1 x np1 matrix
    !vks:  Vector of potential at each grid point last dimension fast
    !      Lexigraphical order
    !dx    Grid spacing
    !yout
    !nev How many eigenvalues to calculate
    type(systemparameters), intent(in) :: sysparams
    type(sharedvalues), intent(in) :: sharedvals
    integer, intent(in)  :: niter
    real(8), intent(in) :: Vks(:)!,Vks(np1**npc**nd)
    complex(8), intent(in) :: yin(:)
    complex(8), intent(out) :: yout(:)
    integer :: Np1,ntot,npc,nd
    real(8), allocatable :: T(:,:)!T(np1,np1)
    real(8) :: dtin
    real(8), allocatable :: wffunc(:,:)
    complex(8), allocatable :: vec(:),vecp(:),vecpp(:),q(:,:)
    complex(8), allocatable :: armat(:,:),vr(:,:),vl(:,:),w(:)
    integer :: i,j,k
    complex(8) :: exparmat(niter,niter)
    integer :: info,lwork,iter
    complex(8) :: qcwork(1)
    complex(8), allocatable :: cwork(:),vecs(:)
    real(8), allocatable :: rwork(:)
    integer, allocatable :: pos2(:,:),posrep(:,:)
    real(8) :: dtt, normh
    integer :: ido
    integer :: quantization
    integer :: sumn
    integer, allocatable :: plus(:),map(:,:)
    include 'mkl_blas.fi' 
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
    Np1=sysparams%np1
    ntot=sysparams%ntot
    npc=sysparams%npart
    nd=sysparams%nd
    allocate(T(np1,np1))
    T=sysparams%T
    dtin=sysparams%dt
    quantization=sysparams%quantization
    if (quantization==2.and.npc.gt.1) then
       do i=1,np1
          T(i,i)=T(i,i)+sharedvals%v1(i)
       end do
    end if

    if (npc.gt.1.and.quantization==2) then
       call makefocklist(np1,npc,posrep,pos2,ido)
       if (ido.ne.ntot) then
          stop
       end if
       allocate(plus(np1))
       call calcmap(np1,ntot,pos2,sumn,plus,map)
    end if

    
    dtt=dtin



    allocate(vecpp(ntot),vecp(ntot),vec(ntot),armat(niter,niter))
    allocate(wffunc(ntot,1))
    allocate(vl(1,niter),vr(niter,niter),w(niter))
    !force asymmetry under exchange of wavefunction
    !if two electrons
    if (npc==2.and.sysparams%quantization==1) then
       call forcesymc(vecpp,Np1,Ntot)
    end if
    roswitch=0
    allocate(reord(ntot,3))
    printit=.false.
    !ntsteps=1 means propagate dt
    !ntsteps=0 means time independent

    
    iter=niter
    armat=0.d0
    vecp=yin!q(:,1)
    allocate(q(ntot,niter+1))
    normh=dsqrt(dble(zdotc(ntot,vecp,1,vecp,1)))
    vecp=vecp/normh
    q(:,1)=vecp
    iter=1
    

    !print*,'step',tsteps
    arnoldi:do
       
       if (npc==2.and.sysparams%quantization==1) then
          call forcesymc(q(:,iter),Np1,Ntot)
       end if
       vecp=q(:,iter)
       ido=-1
       

       !hamiltonian times vector section
       if (ido.eq.-1.or.ido.eq.1) then
          
          !Apply hamiltonian to vecp
          !vec=H*vecp
          !if (nd*npc==1) then
          !   vec=matmul(T,vecp)+vks*vecp
          !elseif (nd*npc==2) then
          !   vec=matmulac2(T,vecp,ntot)+vks*vecp
          !elseif (nd*npc==3) then
          !   vec=matmulac(T,vecp,ntot)+vks*vecp
          !else
          !   print*,'choice of dimension and particle not coded yet'
          !   print*,'stopping'
          !   stop
          !end if
          if (quantization==1) then
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
          elseif (quantization==2) then
             if (npc==1) then
                vec=matmul(T,vecp)+vks*vecp
             else
                vec=onematvecmap(T,Vks,pos2,posrep,npc,np1,ntot,&
                    plus,sumn,map,vecp)
             end if
          end if

          
          

          !Build short iterative arnoldi matrix and vectors
          
          q(:,iter+1)=vec
          if (npc==2.and.sysparams%quantization==1) then
             call forcesymc(q(:,iter+1),Np1,Ntot)
          end if
          do j=1,iter
             armat(j,iter)=zdotc(ntot,q(:,j),1,q(:,iter+1),1)
             !q(:,iter+1)=q(:,iter+1)-q(:,j)*armat(j,iter)
             call zvplusab(q(:,iter+1),-armat(j,iter),q(:,j),ntot)
          end do
          if (npc==2.and.sysparams%quantization==1) then
             call forcesymc(q(:,iter+1),Np1,Ntot)
          end if
          if (iter+1.le.niter) then
             armat(iter+1,iter)=dsqrt(dble(zdotc(ntot,q(:,iter+1),1,q(:,iter+1),1)))
             q(:,iter+1)=q(:,iter+1)/armat(iter+1,iter)
             iter=iter+1
          else
             exit arnoldi
          end if
          
       else
          !print*,'ido',ido
          exit arnoldi
       end if
    end do arnoldi

    !calculate eigenvalues/eigenvectors of arnoldi matrix
    lwork=-1
    allocate(rwork(2*niter))
    call zgeev('N','V',niter,armat,niter,w,vl,1,vr,niter,qcwork,lwork,rwork,info)
    lwork=nint(real(qcwork(1)))
    allocate(cwork(lwork))
    call zgeev('N','V',niter,armat,niter,w,vl,1,vr,niter,cwork,lwork,rwork,info)
    deallocate(cwork,rwork)
    !Generate exponential of matrix times i*dt
    call expdtmatrix(niter,w,vr,dtt,exparmat)
    !propagate vector by dt
    allocate(vecs(niter))
    vecs=0.d0
    vecs(1)=1.d0
    vecs=matmul(exparmat,vecs)
    k=0
    do j=1,ntot
       !vecpp(j)=(dot_product(dconjg(q(j,1:niter)),vecs))
       vecpp(j)=zdotu(niter,q(j,1:niter),1,vecs,1)
    end do
    !ensure propagated vector has correct symmetry if two electrons
    if (npc==2.and.sysparams%quantization==1) then
       call forcesymc(vecpp,Np1,Ntot)
    end if
    vecpp=vecpp/dsqrt(dble(zdotc(ntot,vecpp,1,vecpp,1)))*normh
    deallocate(vecs)
    wffunc(1:ntot,1)=dble(vecpp)
    deallocate(q)
    yout(1:ntot)=vecpp(1:ntot)
    
 
    



    !generate 1D cuts of eigenvectors for neigs desired eigenvalues
    
    if (npc.gt.1.and.quantization==2) then
       deallocate(plus,map)
    end if


    deallocate(reord)
    deallocate(vecpp,vecp,vec,armat)
    deallocate(wffunc)
    deallocate(vl,vr,w)

  end subroutine advancewf 

  subroutine advancewftd(sysparams,sharedvals,niter,vin,vks,vksnew,dvks,dvksnew,yin,yout)
   use derivedtypes
   use hamiltonian_mod
   type(systemparameters), intent(in) :: sysparams
   type(sharedvalues), intent(in) :: sharedvals
   integer, intent(in)  :: niter
   real(8), intent(in) :: vin(:),Vks(:),vksnew(:)!,Vks(np1**npc**nd)
   real(8), intent(in), optional :: dvks(:),dvksnew(:)
   complex(8), intent(in) :: yin(:)
   complex(8), intent(out) :: yout(:)
   real(8), allocatable :: vksp(:)
   complex(8), allocatable :: phiminus(:),phiplus(:),f(:),phinew(:)
   complex(8), parameter :: j1=dcmplx(0.d0,1.d0)
   integer :: i,np1,nd,ntot
   real(8) :: dt

   np1=sysparams%Np1
   ntot=sysparams%ntot
   dt=sysparams%dt
   nd=sysparams%nd
   allocate(vksp(sysparams%ntot),phiminus(sysparams%ntot),phiplus(sysparams%ntot))
   allocate(phinew(sysparams%ntot),f(ntot))
   !do first part of propagation 
   if (present(dvks)) then
      f=calcf1(np1,nd,ntot,dt,sysparams%T,vin,vks,dvks,yin)
   else
      f=0.d0
   end if
   phiminus=yin-dt*j1/2.d0*vks*yin-f
   call advancewf(sysparams,sharedvals,niter,vin,phiminus,phiplus)
   
   !iterate solution
   if (present(dvksnew)) then
      f=calcf1(np1,nd,ntot,dt,sysparams%T,vin,vksnew,dvksnew,yin)
      yout=yin
      phinew=(phiplus+f)/(1.d0+dt*j1/2.d0*vksnew)
      td:do i=1,30
         if (abs(dot_product(yout-phinew,yout-phinew)).lt.1.d-32) then
            exit td
         else
            yout=phinew
            f=calcf1(np1,nd,ntot,dt,sysparams%T,vin,vksnew,dvksnew,yout)
            phinew=(phiplus+f)/(1.d0+dt*j1/2.d0*vksnew)
         end if
      end do td
   else
      yout=(phiplus+f)/(1.d0+dt*j1/2.d0*vksnew)
   end if

  end subroutine advancewftd


end module propagate 
