module secondquant_wrap
  use secondquant_mod
  implicit none
contains
  subroutine onerdm(sysparams,numorb,psi,rdm)
   use derivedtypes
   use nchoosekmod
   !function diagonerdm(pos2,psi,numorb,nfock) result(rdm)
   type(systemparameters), intent(in) :: sysparams
    integer, intent(in) :: numorb
    complex(8), intent(in) :: psi(:)
    complex(8), intent(out) :: rdm(numorb,numorb)
    integer :: i,o,op,fi(numorb),fk(numorb),k,nab, nfock
    real(8) :: si
    integer, ALLOCATABLE :: pos2(:,:), posrep(:,:)
    nfock=sysparams%ntot
    
    !nab=sysparams%npart
    if (sysparams%nalpha == 0 .and. sysparams%nbeta == 0) then
      call makefocklist(sysparams%np1,sysparams%npart,posrep,pos2,nfock)
    else 
      call makefocklist(sysparams%np1,sysparams%npart,posrep,pos2,nfock, sysparams%nalpha, sysparams%nbeta)
    end if
    nab=sum(pos2(:,1))
    rdm=dcmplx(0.d0,0.d0)
    do o=1,nfock
       do i=1,numorb
          fi=pos2(:,o)
          
          if (fi(i)==1) then
             fi(i)=0
             do op=1,nfock
                
                do k=1,numorb
                   fk(:)=pos2(:,op)
                   if (fk(k)==1) then
                      fk(k)=0
                      
                      if (sum(fi*fk)==nab-1) then
                         si=(-1.d0)**sum(fi(1:i-1))*(-1.d0)**sum(fk(1:k-1))
                         rdm(i,k)=rdm(i,k)+psi(o)*dconjg(psi(op))*si
                      end if
                   end if
                end do
             end do
          end if
       end do
    end do
  end subroutine onerdm

  subroutine diagonerdm(sysparams,numorb,psi,rdm)
    use derivedtypes
    use nchoosekmod
    !function diagonerdm(pos2,psi,numorb,nfock) result(rdm)
    type(systemparameters), intent(in) :: sysparams
    integer, intent(in) :: numorb
    integer :: nfock
    complex(8), intent(in) :: psi(:)
    real(8), intent(out) :: rdm(numorb)
    integer :: i,o,op,fi(numorb)
    integer, allocatable :: pos2(:,:),posrep(:,:)
    real(8) :: si
    nfock=sysparams%ntot
    
    !nab=sysparams%npart
    !call makefocklist(sysparams%np1,sysparams%npart,posrep,pos2,nfock)
    if (sysparams%nalpha == 0 .and. sysparams%nbeta == 0) then
      call makefocklist(sysparams%np1,sysparams%npart,posrep,pos2,nfock)
    else 
      call makefocklist(sysparams%np1,sysparams%npart,posrep,pos2,nfock, sysparams%nalpha, sysparams%nbeta)
    end if
    rdm=0.d0
    do o=1,nfock
       !print*,sysparams%pos2(:,o)
       fi=pos2(:,o)
       do i=1,numorb
          if (fi(i)==1) then
             !fi(i)=0
             do op=o,o
                !fk(:)=pos2(:,op)
                !do k=i,i
                   !if (fk(k)==1) then
                   !   fk(k)=0
                   !   if (sum(fi*fk)==nab-1) then
                         si=1.d0!(-1.d0)**sum(fi(1:i-1))*(-1.d0)**sum(fi(1:i-1))
                         rdm(i)=rdm(i)+dble(psi(o)*dconjg(psi(op)))*si
                   !   end if
                   !end if
                !end do
             end do
          end if
       end do
    end do
    !print*,rdm
    !pause
  end subroutine diagonerdm

  !subroutine diagdonerdm(pos2,posrep,psi,H1a,V2s,nab,numorb,nfock) result(rdm)
  subroutine diagdonerdm(sysparams,sharedvals,numorb,psi,rdm)
    use derivedtypes
    use nchoosekmod
    use mapmod
    type(systemparameters), intent(in) :: sysparams
    type(sharedvalues), intent(in) :: sharedvals
    integer, intent(in) :: numorb
    complex(8), intent(in) :: psi(:)
    real(8), intent(out) :: rdm(numorb)
    integer :: nfock
    !integer :: pos2(numorb,nfock),posrep(nab,nfock)
    !real(8) :: H1a(numorb,numorb),V2s(numorb**2)
    complex(8), allocatable :: dpsi(:)
    real(8), allocatable :: H1(:,:)
    integer :: i,o,op,fi(numorb)
    integer, allocatable :: pos2(:,:),posrep(:,:)
    integer, allocatable :: map(:,:),plus(:)
    integer :: sumn
    nfock=sysparams%ntot
    allocate(dpsi(nfock))
    allocate(H1(sysparams%np1,sysparams%np1))
    H1=sysparams%T
    do i=1,sysparams%np1
       H1(i,i)=H1(i,i)+sharedvals%v1(i)
    end do
    !call makefocklist(sysparams%np1,sysparams%npart,posrep,pos2,nfock)
    if (sysparams%nalpha == 0 .and. sysparams%nbeta == 0) then
      call makefocklist(sysparams%np1,sysparams%npart,posrep,pos2,nfock)
    else 
      call makefocklist(sysparams%np1,sysparams%npart,posrep,pos2,nfock, sysparams%nalpha, sysparams%nbeta)
    end if
    allocate(plus(sysparams%np1))
    call calcmap(sysparams%ntot1,sysparams%ntot,pos2,sumn,plus,map)
    !dpsi=dcmplx(0.d0,1.d0)*onematvec(H1,sharedvals%vinteract,pos2,&
    !     posrep,sysparams%npart,numorb,nfock,psi)
    dpsi=dcmplx(0.d0,1.d0)*onematvecmap(H1,sharedvals%vinteract,pos2,&
                posrep,sysparams%npart,numorb,nfock,&
                    plus,sumn,map,psi)
    deallocate(map,posrep)
    rdm=0.d0
    do o=1,nfock
       fi=pos2(:,o)
       do i=1,numorb
          if (fi(i)==1) then
             !fi(i)=0
             do op=o,o
                !fk(:)=pos2(:,op)
                !do k=i,i
                   !if (fk(k)==1) then
                   !   fk(k)=0
                      !if (sum(fi*fk)==nab-1) then
                         !si=1.d0!(-1.d0)*sum(fi(1:i-1))*(-1.d0)*sum(fk(1:k-1))
                         rdm(i)=rdm(i)-dble(dpsi(o)*dconjg(psi(op))+psi(o)*dconjg(dpsi(op)))
                      !end if
                   !end if
                !end do
             end do
          end if
       end do
    end do
  end subroutine diagdonerdm

  !function diagddonerdm(pos2,posrep,psi,H1a,V2s,nab,numorb,nfock) result(rdm)
  subroutine diagddonerdm(sysparams,sharedvals,numorb,psi,rdm)
    use derivedtypes
    use nchoosekmod
    use mapmod
    type(systemparameters), intent(in) :: sysparams
    type(sharedvalues), intent(in) :: sharedvals
    integer, intent(in) :: numorb
    complex(8), intent(in) :: psi(:)
    real(8), intent(out) :: rdm(numorb)
    complex(8), allocatable :: dpsi(:),ddpsi(:)
    real(8), allocatable :: H1(:,:)
    integer :: nfock,i,o,op,fi(numorb)
    integer, allocatable :: pos2(:,:),posrep(:,:)
    integer, allocatable :: map(:,:),plus(:)
    integer :: sumn
    nfock=sysparams%ntot
    allocate(dpsi(nfock),ddpsi(nfock))
    allocate(H1(sysparams%np1,sysparams%np1))
    H1=sysparams%T
    do i=1,sysparams%np1
       H1(i,i)=H1(i,i)+sharedvals%v1(i)
    end do
    !call makefocklist(sysparams%np1,sysparams%npart,posrep,pos2,nfock)
    if (sysparams%nalpha == 0 .and. sysparams%nbeta == 0) then
      call makefocklist(sysparams%np1,sysparams%npart,posrep,pos2,nfock)
    else 
      call makefocklist(sysparams%np1,sysparams%npart,posrep,pos2,nfock, sysparams%nalpha, sysparams%nbeta)
    end if
    allocate(plus(sysparams%np1))
    call calcmap(sysparams%ntot1,sysparams%ntot,pos2,sumn,plus,map)
    !dpsi=onematvec(H1,sharedvals%vinteract,pos2,posrep,&
    !     sysparams%npart,numorb,nfock,psi)
    !ddpsi=onematvec(H1,sharedvals%vinteract,pos2,posrep,&
    !     sysparams%npart,numorb,nfock,dpsi)
    dpsi=onematvecmap(H1,sharedvals%vinteract,pos2,&
                posrep,sysparams%npart,numorb,nfock,&
                    plus,sumn,map,psi)
    ddpsi=onematvecmap(H1,sharedvals%vinteract,pos2,&
                posrep,sysparams%npart,numorb,nfock,&
                    plus,sumn,map,dpsi)
    
    rdm=0.d0
    do o=1,nfock
       fi=pos2(:,o)
       do i=1,numorb
          if (fi(i)==1) then
             !fi(i)=0
             do op=o,o
                !fk(:)=pos2(:,op)
                !do k=i,i
                   !if (fk(k)==1) then
                   !   fk(k)=0
                      !if (sum(fi*fk)==nab-1) then
                         !si=1.d0!(-1.d0)*sum(fi(1:i-1))*(-1.d0)*sum(fk(1:k-1))
                rdm(i)=rdm(i)-dble(ddpsi(o)*dconjg(psi(op))+psi(o)*dconjg(ddpsi(op))-&
                     dconjg(dpsi(o))*dpsi(op)-dpsi(o)*dconjg(dpsi(op)))
                      !end if
                   !end if
                !end do
             end do
          end if
       end do
    end do
  end subroutine diagddonerdm

  subroutine calc_ground_state(sysparams,sharedvals,numorb,niter,psi)
   use derivedtypes
   use nchoosekmod
   use mapmod
   type(systemparameters), intent(in) :: sysparams
   type(sharedvalues), intent(in) :: sharedvals
   integer, intent(in) :: numorb, niter
   complex(8), intent(inout) :: psi(:)
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
   integer Np1,ntot,npc,nd
   real(8), allocatable :: T(:,:)!T(np1,np1)
   integer :: neigs, nev, nfock
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
   real(8), allocatable :: DR(:),DI(:),Zeig(:,:),workev(:), H1(:,:)
   integer, allocatable :: pos2(:,:),posrep(:,:)
   integer, allocatable :: map(:,:),plus(:)
   integer :: sumn
   include 'mkl_blas.fi' 
   Np1=sysparams%np1
   ntot=sysparams%ntot
   npc=sysparams%npart
   nd=sysparams%nd
   allocate(T(np1,np1))
   T=sysparams%T
   nev = 5
   !print*, "infunction"
   nfock=sysparams%ntot
   allocate(H1(sysparams%np1,sysparams%np1))
   H1=sysparams%T
   do i=1,sysparams%np1
      H1(i,i)=H1(i,i)+sharedvals%v1(i)
   end do
   !call makefocklist(sysparams%np1,sysparams%npart,posrep,pos2,nfock)
   if (sysparams%nalpha == 0 .and. sysparams%nbeta == 0) then
      call makefocklist(sysparams%np1,sysparams%npart,posrep,pos2,nfock)
    else 
      call makefocklist(sysparams%np1,sysparams%npart,posrep,pos2,nfock, sysparams%nalpha, sysparams%nbeta)
    end if
   allocate(plus(sysparams%np1))
   call calcmap(sysparams%ntot1,sysparams%ntot,pos2,sumn,plus,map)
   !dpsi=onematvec(H1,sharedvals%vinteract,pos2,posrep,&
   !     sysparams%npart,numorb,nfock,psi)
   !ddpsi=onematvec(H1,sharedvals%vinteract,pos2,posrep,&
   !     sysparams%npart,numorb,nfock,dpsi)
   neigs=nev



   allocate(vecpp(ntot),vecp(ntot),vec(ntot))
   allocate(wffunc(ntot,Neigs))
   allocate(veigs(neigs+1))
   
   

   
   iter=niter
   allocate(qs(ntot,niter))
   call random_number(qs(:,1))
   vecp=qs(:,1)
   vecp=vecp/dsqrt(dble(zdotc(ntot,vecp,1,vecp,1)))
   qs(:,1)=dble(vecp)
   
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
   info=0
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
         dsp=ipntr(2)
         dfp=ipntr(2)+ntot-1
         vecp=workd(ds:df)
      
         vec = onematvecmap(H1,sharedvals%vinteract,pos2,&
                posrep,sysparams%npart,numorb,nfock,&
                    plus,sumn,map,vecp)


         workd(dsp:dfp)=dble(vec)
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
   !call dneupd(rvec,howmny,select,dr,di,zeig,ntot,sigmar,sigmai,workev,&
   !     'I',ntot,'SR',Neigs,tol,resid,iter,qs(:,1:iter),ntot,iparam,&
   !     ipntr,workd,workl,lworkl,info)
   call dseupd(rvec,howmny,select,dr,zeig,ntot,sigmar,&
        'I',ntot,'SA',Neigs,tol,resid,iter,qs(:,1:iter),ntot,iparam,&
        ipntr,workd,workl,lworkl,info)
   print*, 'eig', dr
   
   wffunc=zeig(:,1:neigs) 
   deallocate(DR,DI,Zeig,workev)
   deallocate(select)
   deallocate(qs)

   psi(:)=wffunc(:,1)
      
   
   
   deallocate(iparam)
   deallocate(ipntr)
   deallocate(workl)
   deallocate(resid)
   deallocate(workd)
   



   deallocate(vecpp,vecp,vec)
   deallocate(wffunc)
   deallocate(veigs)

 end subroutine calc_ground_state


end module secondquant_wrap
