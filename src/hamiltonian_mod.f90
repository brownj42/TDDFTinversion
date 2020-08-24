Module hamiltonian_mod
  implicit none
  integer :: fas
  logical :: printit,evens,odds,use2
  logical :: singlet,triplet
contains

  !subroutine calcynew(Np1,ntot,npc,nd,niter,ntsteps,T,vks,dx,dtin,yin,yout,nev)
  !subroutine calceigenstates(Np1,ntot,npc,nd,niter,T,vks,dx,nev,yout)
  subroutine calceigenstates(sysparams,vks,niter,nev,yout)
    use derivedtypes
    use sortmod
    use matmul_mod
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
    !yout  If ntsteps=0, lowest energy triplet eigenvector
    !      if ntsteps>0  vector after dtin timestep
    !neigs How many eigenvalues to calculate if ntsteps=0
    type(systemparameters), intent(in) :: sysparams
    integer, intent(in) :: niter
    real(8), intent(in) :: Vks(:)!,Vks(np1**npc**nd)
    integer, intent(in):: nev
    complex(8), intent(out) :: yout(:)
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
       !print*,(dr(1:Neigs)+2.d0)
       !print*,(di(1:Neigs))!*219474.63
       !print*,(dr(1:Neigs)+0.d0)
    else
       veigs=dr
       print*,(dr(1:Neigs)+0.d0)
       !print*,(di(1:Neigs))!*219474.63
    end if
    
    wffunc=zeig(:,1:neigs) 
    deallocate(DR,DI,Zeig,workev)
    deallocate(select)
    deallocate(qs)
    !write(filename,'(A6,i2.2,A4)')'eigene',nd,'.dat'
    !open(unit=23,file=filename)
    !do i=1,Neigs
    !   write(23,'(f20.10,2x,f20.10)')veigs(i),veigs(i)-veigs(1)
    !end do
    !close(23)



    !generate 1D cuts of eigenvectors for neigs desired eigenvalues
    allocate(valwf(np1,nd,Neigs))
    valwf=0.d0
    do i=1,ntot
       do c=1,nd
          s=ind(i,c)
          valwf(s,c,:)=valwf(s,c,:)+real(wffunc(i,:)) 
       end do
    end do
    
    !transfer ground state eigenvector to yout
    if (npc==2) then
       call forcesym(wffunc(1:ntot,neigs),Np1,Ntot)
       do c=1,ntot
          yout(c)=dcmplx((wffunc(c,neigs)),0.d0)!+wffunc(c,2))/2.d0!sum(wffunc(c,1:Neigs))
       end do
    else
      
       if (nev==1) then
          yout=wffunc(:,1)
       else
          if (veigs(1).gt.veigs(2)) then
             if (use2) then
                yout=wffunc(:,1)
             else
                yout=wffunc(:,2)
             end if
          else
             if (use2) then
                yout=wffunc(:,2)
             else
                yout=wffunc(:,1)
             end if
          end if
       end if
      
    
    end if
    deallocate(iparam)
    deallocate(ipntr)
    deallocate(workl)
    deallocate(resid)
    deallocate(workd)
    deallocate(valwf)
    


    deallocate(reord)
    deallocate(vecpp,vecp,vec)
    deallocate(wffunc)
    deallocate(veigs)
    deallocate(ind)
    deallocate(pind)

  end subroutine calceigenstates


  
  !calculate Exp(-i dt A) 
  subroutine expdtmatrix(niter,w,vr,dt,exph)
    integer :: niter
    complex(8) :: w(niter),vr(niter,niter),exph(niter,niter)
    complex(8) :: vri(niter,niter),dw(niter,niter)
    real(8) :: dt
    integer :: info
    integer :: ipiv(niter)
    complex(8), allocatable :: work(:)
    complex(8) :: qwork(1)
    integer :: lwork
    integer :: i

    vri=vr
    call zgetrf(niter,niter,vri,niter,ipiv,info)
    lwork=-1
    call zgetri(niter,vri,niter,ipiv,qwork,lwork,info)
    lwork=nint(real(qwork(1)))
    allocate(work(lwork))
    call zgetri(niter,vri,niter,ipiv,work,lwork,info)
    deallocate(work)
    !print*,matmul(vri,vr)
    dw=0.d0
    do i=1,niter
       dw(i,i)=zexp(-dt*dcmplx(0.d0,1.d0)*w(i))
    end do
    exph=matmul(matmul(vr,dw),vri)
  end subroutine expdtmatrix




  real(8) function delta(i,j)
    integer :: i,j
    if (i==j) then
       delta=1.d0
    else
       delta=0.d0
    end if
  end function delta

  integer function pos(N,nb,nc)
    integer :: N
    integer :: nb(N),nc(N)
    integer :: phold
    integer :: i
    phold=1
    do i=1,N
       phold=phold+nb(i)*(nc(i)-1)
    end do
    pos=phold
  end function pos

  !force exchange anti-symmetry for real vector
  subroutine forcesym(a,N,N2)
    integer,intent(in)  :: N,N2
    real(8), intent(inout) :: a(n2)          
    integer :: i1,i2,pos1,pos2
    real(8) :: av
    !print*,n,n2,size(a)
    do i1=1,n
       do i2=1,n
          pos1=(i1-1)*n+i2
          pos2=(i2-1)*n+i1
          
          if (pos1==pos2) then
             if (triplet) then
                a(pos1)=0.d0
                a(pos2)=0.d0
             end if
             if (singlet) then
                av=((a(pos1))+(a(pos2)))/2.d0
                a(pos1)=sign(av,av)
                a(pos2)=sign(av,av)
             end if
          else
             if (triplet) then
                av=(a(pos1)-a(pos2))/(2.d0)
                a(pos1)=sign(av,a(pos1))
                a(pos2)=sign(av,-a(pos1))
             end if
             if (singlet) then
                av=((a(pos1))+(a(pos2)))/2.d0
                a(pos1)=sign(av,av)
                a(pos2)=sign(av,av)
             end if
          end if
       end do
    end do
  end subroutine forcesym

  !force even or odd symmetry for 1D problem with respect to center
  subroutine forcesym1(a,n)
    integer, intent(in) :: n
    real(8), intent(inout) :: a(n)
    integer :: i
    real(8) :: hold

    if (odds) then
       do i=1,n/2
          hold=(dabs(a(i))+dabs(a(n+1-i)))/2.d0
          a(i)=dsign(hold,a(i))
          a(n+1-i)=dsign(hold,-a(i))
       end do
       a((n+1)/2)=0.d0
    end if
    if (evens) then
       do i=1,n/2
          hold=(dabs(a(i))+dabs(a(n+1-i)))/2.d0
          a(i)=dsign(hold,a(i))
          a(n+1-i)=dsign(hold,a(i))
       end do
    end if
  end subroutine forcesym1

  !force anti-symmetric exchange for 2D complex vector
  subroutine forcesymc(a,N,N2)
    integer,intent(in)  :: N,N2
    complex(8), intent(inout) :: a(n2)          
    integer :: i1,i2,pos1,pos2
    complex(8) :: av
    !print*,n,n2,size(a)
    do i1=1,n
       do i2=1,n
          pos1=(i1-1)*n+i2
          pos2=(i2-1)*n+i1
          
          if (pos1==pos2) then
             if (triplet) then
                a(pos1)=0.d0
                a(pos2)=0.d0
             end if
             if (singlet) then
                av=((a(pos1))+(a(pos2)))/2.d0
                a(pos1)=dcmplx(sign(dble(av),dble(av)),sign(dimag(av),dimag(av)))
                a(pos2)=dcmplx(sign(dble(av),dble(av)),sign(dimag(av),dimag(av)))
             end if
          else
             if (triplet) then
                 av=(a(pos1)-a(pos2))/(2.d0)
                 a(pos1)=dcmplx(sign(dble(av),dble(av)),sign(dimag(av),dimag(av)))
                 a(pos2)=dcmplx(sign(dble(av),dble(-av)),sign(dimag(av),dimag(-av)))
             end if
             if (singlet) then
                 av=((a(pos1))+(a(pos2)))/2.d0
                 a(pos1)=dcmplx(sign(dble(av),dble(av)),sign(dimag(av),dimag(av)))
                 a(pos2)=dcmplx(sign(dble(av),dble(av)),sign(dimag(av),dimag(av)))
             end if
          end if
       end do
    end do
  end subroutine forcesymc

  function calcf1(np1,nd,ntot1,dt,T,vin,vks,dvks,phi) result(f)
    use matmul_mod
    integer :: np1,nd,ntot1
    real(8) :: dt,T(np1,np1),vin(ntot1),vks(ntot1),dvks(ntot1)
    complex(8) :: phi(ntot1)
    complex(8) :: f(ntot1)
    complex(8) :: ay(ntot1),nwvd(ntot1,0:1),Tnwvd0(ntot1)
    roswitch=0
    allocate(reord(ntot1,3))
    nwvd(:,0)=vks*phi !N
    if (nd==1) then
       Tnwvd0=matmul(T,nwvd(:,0))
       ay=-dcmplx(0.d0,1.d0)*(matmul(T,phi)+vks*phi+vin*phi)
    elseif(nd==3) then
       Tnwvd0=matmulac(T,nwvd(:,0),ntot1)
       ay=-dcmplx(0.d0,1.d0)*(matmulac(T,phi,ntot1)+vks*phi+vin*phi)
    end if
    deallocate(reord)
    nwvd(:,1)=dvks*phi 
    nwvd(:,1)=vks*ay+nwvd(:,1) !Nprime

    f=dt**2*dcmplx(1.d0,0.d0)/12.d0*(dcmplx(-1.d0,0.d0)*&
         (Tnwvd0+vin*nwvd(:,0))+dcmplx(0.d0,1.d0)*&
         nwvd(:,1))
  end function calcf1

  function calcvhar(ntot1,npart,ntot2,vinteract,dp) result(vhar)
     integer :: ntot1,npart,ntot2
     real(8) :: vinteract(ntot2),dp(ntot1)
     real(8) :: vhar(ntot1)
     integer :: i,j
     if (npart.gt.1) then
        do i=1,ntot1
           vhar(i)=0.d0
           do j=1,ntot1
              vhar(i)=vhar(i)+dp(j)*vinteract((j-1)*ntot1+i)
           end do
        end do
     else
        vhar=0.d0
     end if
     
   end function calcvhar
end Module hamiltonian_mod
