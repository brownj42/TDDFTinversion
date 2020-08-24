Module pot_invert_mod 
  use sortmod
  implicit none

contains
  
  

    
  !Orthogonalizes a set of vectors A
  subroutine orth(A,m,n)
    integer, intent(in) :: m,n
    real(8), intent(inout) :: A(m,n)
    integer :: k
    integer :: LDA
    real(8), allocatable :: tau(:),work(:),U(:,:),VT(:,:)
    integer, allocatable :: iwork(:)
    real(8) :: qwork(1)
    integer :: lwork,info,LDU,LDVT,ucol
    lda=m
    allocate(tau(min(m,n)))
    k=min(m,n)
    LDU=M
    UCOL=N
    LDVT=min(M,N)
    allocate(U(LDU,UCOL),VT(LDVT,N),iwork(8*min(M,N)))


    lwork=-1
    call dgesdd('S',M,N,A,LDA,tau,U,LDU,VT,LDVT,qwork,lwork,iwork,info)
    lwork=nint(qwork(1))
    allocate(work(lwork))
    call dgesdd('S',M,N,A,LDA,tau,U,LDU,VT,LDVT,work,lwork,iwork,info)
    deallocate(work)
    A=U
    deallocate(U,iwork,vt,tau)

  end subroutine orth
  

  !Apply's the force balance matrix to a vector Eq.(22)
  function matmulW(T,ry,iy,dy,x,n3) result(xp)
    use matmul_mod
    real(8), intent(in) :: T(:,:),ry(:),iy(:),dy(:),x(:)
    integer, intent(in) :: n3
    real(8) :: xp(n3)

    !apply real vectors r
    xp=-2.d0*ry*matmula(T,x*ry,n3)
    !apply imaginary vectors i
    xp=xp-2.d0*iy*matmula(T,x*iy,n3)
    xp=xp-dy*x
  end function matmulW
  






  !Subroutine for morse-penrose pseudoinverse with cutoff parameter
  subroutine pinv(Ah,m,n,cutoff,Apinv)
      integer, intent(in) :: m,n
      real(8), intent(in) :: Ah(m,n),cutoff
      real(8), intent(out) :: Apinv(n,m)
      real(8) :: a(m,n)
      integer :: i
      character(1) :: JOBU,JOBVT
      integer :: lda
      real(8), allocatable :: S(:),U(:,:),Sigma(:,:)
      integer :: LDU
      real(8), allocatable :: Vt(:,:)
      integer :: LDVT
      real(8), allocatable :: work(:)
      real(8) :: qwork(1)
      integer :: lwork
      integer :: info
      
      a=ah
      apinv=a
      JOBU='A'
      JOBVT='A'
      LDA=max(1,M)
      LDU=M
      LDVT=N
      allocate(U(M,M),sigma(M,N),VT(N,N),s(min(m,n)))
      lwork=-1
      call dgesvd(jobu,jobvt,m,n,A,lda,s,u,ldu,vt,ldvt,qwork,lwork,info)
      lwork=nint(qwork(1))
      allocate(work(lwork))
      call dgesvd(jobu,jobvt,m,n,A,lda,s,u,ldu,vt,ldvt,work,lwork,info)
      a=apinv
      deallocate(work)
      sigma=0.d0
      do i=1,min(m,n)
         if(abs(s(i)).ge.maxval(abs(s))*1.d0*abs(cutoff)) then       
             sigma(i,i)=1.d0/s(i)
         end if
      end do
      Apinv=matmul(matmul(transpose(Vt),transpose(sigma)),transpose(U))

  end subroutine pinv

  !Calculates inverted potential
  subroutine calcvks(niter,N,nd,N3,npart,T,pinv0minresqlp1,ddnx,yks,vks)
    use minresqlpdatamodule
    use minresqlpmodule
    use matmul_mod
    !lattice basis size one dimension, N^3, max number of iterations for
    !MINRES-QLP, number of particles
    integer, intent(in) :: N,nd,n3,niter,npart,pinv0minresqlp1 
    !Kinetic energy operator matrix, second derivative of density
    real(8), intent(in) :: T(:,:),ddnx(:)
    !current KS wavefunctions for npart particles
    complex(8), intent(in) :: yks(:,:)
    !Current V^KS potential
    real(8), intent(inout) :: vks(:)
    integer :: i,p,np
    real(8) :: maxit
    real(8), allocatable :: S(:),ry(:,:),iy(:,:),diagy(:),Mm1(:)
    real(8), allocatable :: Km(:,:),G(:,:),w(:)
    real(8) :: po(n3,1)
    complex(8), allocatable :: Gy(:,:)
    complex(8), allocatable :: ty(:),tty(:)
    include 'mkl_blas.fi' 
    allocate(ty(n3),tty(n3))
    allocate(s(n3),ry(n3,npart),iy(n3,npart),diagy(n3))
    allocate(MM1(n3))
    roswitch=0
    allocate(reord(n3,3))
    npartw=npart
    !calculate momentum-stress tensor and diagonal term of force-balance matrix
    s=0.d0
    diagy=0.d0
    do i=1,npart
       if (nd==1) then
          ty=matmul(T,yks(:,i))
          tty=matmul(T,ty)
          !tty=matmul(T2,yks(:,i))
       elseif(nd==3) then
          ty=matmulac(T,yks,n3)
          tty=matmulac(T,ty,n3)
       end if
       s=s-2.d0*(dble(ty)**2+dimag(ty)**2-&
            0.5d0*(dble(dconjg(yks(:,i))*tty)+dble(yks(:,i)*dconjg(tty))))
       diagy=diagy-1.d0*dble(dconjg(ty)*(yks(:,i))+dconjg(yks(:,i))*ty)
    end do
    s=ddnx+s
    deallocate(ty,tty)
    
    !project out constant vector from right side of equation
    !po=1.d0/dsqrt(dble(n3))
    !s=s-ddot(n3,s,1,po(:,1),1)*po(:,1)    
    
    
    !decompose KS wavefunctions into real and imaginary parts
    ry=dble(yks)
    iy=dimag(yks)

    !calculate diagonal elements for force-balance matrix for preconditioner
    do p=1,n3
       mm1(p)=0.d0
       do np=1,npart
          mm1(p)=mm1(p)-2.d0*(ry(p,np)**2+iy(p,np)**2)*(nd)!*5.d0!*T(1,1)
       end do
       mm1(p)=mm1(p)-diagy(p)
    end do
    maxit=maxval(dabs(mm1))
    mm1=mm1/maxit
    maxit=1.d0
    !generate preconditioner 
    p=0
    do i=1,n3
       p=p+1
       if (abs(mm1(p)).gt.1.d0*epsilon(maxit)) then
          mm1(p)=1.d0/dabs(mm1(p))
       else
          mm1(p)=1.d0*abs(epsilon(1.d0/(1.d0*epsilon(maxit))))
       end if
    end do
 
    if (pinv0minresqlp1==0.and.nd==1) then       
       allocate(G(n3,n3),Km(n3,n3),gy(n3,1),w(n3))
       !since KS wavefunctions are in one dimension, generate full force-balance
       !matrix and preconditionin matrix
       G=0.d0
       do i=1,npart
          gy(:,1)=yks(:,i)
          G=G+2.d0*dble(matmul(gy,dconjg(transpose(gy))))
       end do
       Km=-G*T
       G=0.d0
       do i=1,n
         Km(i,i)=Km(i,i)-diagy(i)
         G(i,i)=sqrt(mm1(i))
       end do
       G=-matmul(G,matmul(km,G))
       !km=g
       !lwork=-1    
       !call dsyev('V','U',n,km,n,w,qwork,lwork,info)
       !lwork=nint(qwork(1))
       !allocate(work(lwork))
       !call dsyev('V','U',n,km,n,w,work,lwork,info)
       !deallocate(work) 
       !print*,'K eigs'
       !print*,w(1:3)
       !print*,w(n-2:n)
    end if
      
    if (pinv0minresqlp1==1) then
       npartw=npart
        allocate(TW(n,n),ryw(n3,npartw),iyw(n3,npartw),dyw(n3),mW(n3))
        mW=abs(mm1)!sqrt(abs(mm1))
        TW=-T
        ryw=ry
        iyw=iy
        dyw=-diagy
    end if
    s=-s
    if (pinv0minresqlp1==1) then
       !use MINRES-QLP for inversion
       if (nd==1) then
          call minresqlp(n=n,aprod=aprodw,b=s,msolve=msolvew,itnlim=2*n,rtol=1.d-13,x=vks,trancond=5.d5)
       else
          
           call minresqlp(n=n3,aprod=aprodw3,b=s,shift=0.d0,msolve=msolvew,itnlim=1000,rtol=1.d-13,x=vks,trancond=5.d5)
       end if
       deallocate(TW,ryw,iyw,dyw,mW)
    end if
    
    if (pinv0minresqlp1==0) then       
       !use Morse-Penrose pseudoinverse for inversion
       if (sum(abs(G)).lt.1.d-62) then
          vks=0.d0
       else
       call pinv(G,n,n,1.d-14,km)
       vks=matmul(km,sqrt(mm1)*s)*sqrt(mm1)
       end if
       deallocate(G,Km,gy,w)
    end if
    deallocate(reord)
    deallocate(S,ry,iy,diagy)
    deallocate(MM1)
  end subroutine calcvks

  !Apply force-balance matrix to a vector
  subroutine AprodW(n,x,y)
    use minresqlpDataModule, only : iyw,ryw,Tw,ip,dp,dyw,npartw
    integer(ip), intent(in)  :: n
    real(dp),    intent(in)  :: x(n)
    real(dp),    intent(out) :: y(n)
    integer :: j
    real(dp) :: p(n)
    !p=1.d0/sqrt(dble(n))
    p=x!-dot_product(x,p)*p

    y=0.d0
    !for each particle
    do j=1,npartw
       !apply real part
       y=y-2.d0*ryW(:,j)*matmul(TW,p*ryW(:,j))
       !apply imaginary part
       y=y-2.d0*iyW(:,j)*matmul(TW,p*iyW(:,j))
    end do
    y=y-dyW*p
    !p=1.d0/sqrt(dble(n))
    !y=y-dot_product(y,p)*p
  end subroutine AprodW
  subroutine AprodW3(n,x,y)
    !use minresqlpDataModule
    use minresqlpDataModule
    use matmul_mod 
    integer(ip), intent(in)  :: n
    real(dp),    intent(in)  :: x(n)
    real(dp),    intent(out) :: y(n)
    integer :: j
    include 'mkl_blas.fi'
    
    y=0_dp
    !for each particle
    do j=1,npartw
       !apply real part
       !y=y-2_dp*ryW(:,j)*matmula(TW,x*ryW(:,j),n)
       call daxpy(n,-2.d0,vtimesv(ryW(:,j),matmula(TW,vtimesv(x,ryW(:,j),n),n),n),1,y,1)
       !apply imaginary part
       !y=y-2_dp*iyW(:,j)*matmula(TW,x*iyW(:,j),n)
       call daxpy(n,-2.d0,vtimesv(iyW(:,j),matmula(TW,vtimesv(x,iyW(:,j),n),n),n),1,y,1)
    end do
    !y=y-dyW*x
    call daxpy(n,-1.d0,vtimesv(dyW(:),x,n),1,y,1)
    
  end subroutine AprodW3


  !apply preconditioner to vector
  subroutine MsolveW(n,x,y)                   ! Solve M*y = x
    use minresqlpDataModule
    integer(ip), intent(in)    :: n
    real(dp),    intent(in)    :: x(n)
    real(dp),    intent(out)   :: y(n)
    y=Mw*x
  end subroutine Msolvew



   subroutine recenter_potential(np1,nd,ntot1,npart,enerh,T,phi,vks)
     use matmul_mod
     integer, intent(in) :: np1,nd,ntot1,npart
     real(8), intent(in) :: enerh,T(np1,np1)
     complex(8), intent(in) :: phi(ntot1,npart)
     real(8), intent(inout) :: vks(ntot1)
     complex(8) :: hphi(ntot1)
     integer :: i
     real(8) :: ener,vshifthold
     include 'mkl_blas.fi'
     roswitch=0
     allocate(reord(ntot1,3))
     ener=0.d0
     do i=1,npart
        if (nd==1) then
           hphi=matmul(T,phi(:,i))
           hphi=hphi+vks*phi(:,i)
        elseif (nd==2) then
           hphi=matmulac2(T,phi(:,i),ntot1)
           call zaxpy(ntot1,dcmplx(1.d0,0.d0),vtimeszv(vks,phi(:,i),ntot1),1,hphi,1)
        elseif (nd==3) then
           hphi=matmulac(T,phi(:,i),ntot1)
           call zaxpy(ntot1,dcmplx(1.d0,0.d0),vtimeszv(vks,phi(:,i),ntot1),1,hphi,1)
        end if
        ener=ener+dble(zdotc(ntot1,phi(:,i),1,hphi,1))
     end do
     deallocate(reord)
     vshifthold=(ener-enerh)/dble(npart)
     vks=vks-vshifthold
   end subroutine recenter_potential
   

   function vks_derivative(ntot1,dt,vks,vksh,vkshh,sw) result(dvks)
     integer :: ntot1,sw
     real(8) :: dt,vks(ntot1),vksh(ntot1),vkshh(ntot1)
     real(8) :: dvks(ntot1)
     integer :: i

     if (sw.gt.0) then
        !$OMP parallel shared(dvks,vks,vksh,ntot1) firstprivate(dt)
        !$OMP do schedule(static,4) 
        do i=1,ntot1
           dvks(i)=(vks(i)-vksh(i))/dt
        end do
        !$OMP end do nowait
        !$OMP end parallel
     else
        !$OMP parallel shared(dvks,vks,vksh,ntot1) firstprivate(dt)
        !$OMP do schedule(static,4) 
        do i=1,ntot1
           dvks(i)=(3.d0*vks(i)-4.d0*vksh(i)+vkshh(i))/dt/2.d0
        end do
        !$OMP end do nowait
        !$OMP end parallel
     end if
   end function vks_derivative

   subroutine restrict_vks(ntot1,dt,dvksmax,dvks,vks,vksh,vkshh)
     integer, intent(in) :: ntot1
     real(8), intent(in) :: dt,dvksmax
     real(8), intent(inout) :: dvks(:),vks(:),vksh(:),vkshh(:)
     integer :: i
   
    !$OMP parallel shared(dvks,vks,vksh,vkshh,ntot1) firstprivate(dt)
    !$OMP do schedule(static,4) 
     do i=1,ntot1
        if (abs(dvks(i)).gt.dvksmax.or.abs(vks(i)).gt.400.d0) then
           if (abs(vks(i)).lt.400.d0) then
              dvks(i)=sign(dvksmax,dvks(i))
              vks(i)=vksh(i)+dvks(i)*dt!+(vks(i-1)-vks(i-2))!(vks(i)+hdvks(i))/2.d0
              vkshh(i)=vksh(i)-dvks(i)*dt!0.d0
           else
              dvks(i)=0.d0
              vks(i)=vksh(i)!.d0
              !vksh(i)=0.d0
              vkshh(i)=vksh(i)!0.d0
           end if
           
        end if
     end do
    !$OMP end do nowait
    !$OMP end parallel
        
   end subroutine restrict_vks

   subroutine calcdvks(np1,nd,ntot1,npart,sw,enerh,dt,T,vin,phi,dvksmax,vkshh,vksh,vks,dvks)
     integer, intent(in) :: np1,nd,ntot1,npart,sw
     real(8), intent(in) :: enerh,dt,T(np1,np1),vin(ntot1),dvksmax
     complex(8), intent(in) :: phi(ntot1,npart)
     real(8), intent(inout) :: vksh(ntot1),vkshh(ntot1)
     real(8), intent(inout) :: vks(ntot1)
     real(8), intent(out) :: dvks(ntot1)
     integer :: i
     call recenter_potential(np1,nd,ntot1,npart,enerh,T,phi,vks)
     do i=1,20
        !vks=vks-vin
        call daxpy(ntot1,-1.d0,vin,1,vks,1)
        dvks=vks_derivative(ntot1,dt,vks,vksh,vkshh,sw)
        call restrict_vks(ntot1,dt,dvksmax,dvks,vks,vksh,vkshh)
        !vks=vks+vin
        call daxpy(ntot1,1.d0,vin,1,vks,1)
        call recenter_potential(np1,nd,ntot1,npart,enerh,T,phi,vks)
     end do
   end subroutine calcdvks

end Module pot_invert_mod 
