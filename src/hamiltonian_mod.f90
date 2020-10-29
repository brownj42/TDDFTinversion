Module hamiltonian_mod
  implicit none
  integer :: fas
  logical :: printit,evens,odds,use2
  logical :: singlet,triplet
contains

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
    elseif(nd==2) then
       Tnwvd0=matmulac2(T,nwvd(:,0),ntot1)
       ay=-dcmplx(0.d0,1.d0)*(matmulac2(T,phi,ntot1)+vks*phi+vin*phi)
    elseif(nd==3) then
       Tnwvd0=matmulac(T,nwvd(:,0),ntot1)
       ay=-dcmplx(0.d0,1.d0)*(matmulac(T,phi,ntot1)+vks*phi+vin*phi)
    else
       print*,'not implemented yet, stopping'
       stop
    end if
    deallocate(reord)
    nwvd(:,1)=dvks*phi 
    nwvd(:,1)=vks*ay+nwvd(:,1) !Nprime

    f=dt**2*dcmplx(1.d0,0.d0)/12.d0*(dcmplx(-1.d0,0.d0)*&
         (Tnwvd0+vin*nwvd(:,0))+dcmplx(0.d0,1.d0)*&
         nwvd(:,1))
  end function calcf1

  function calcf1c(np1,nd,ntot1,dt,T,vin,vks,dvks,phi) result(f)
    use matmul_mod
    integer :: np1,nd,ntot1
    real(8) :: dt,T(np1,np1),vin(ntot1),dvks(ntot1)
    complex(8) :: vks(ntot1)
    complex(8) :: phi(ntot1)
    complex(8) :: f(ntot1)
    complex(8) :: ay(ntot1),nwvd(ntot1,0:1),Tnwvd0(ntot1)
    roswitch=0
    allocate(reord(ntot1,3))
    nwvd(:,0)=vks*phi !N

    if (nd==1) then
       Tnwvd0=matmul(T,nwvd(:,0))
       ay=-dcmplx(0.d0,1.d0)*(matmul(T,phi)+vks*phi+vin*phi)
    elseif(nd==2) then
       Tnwvd0=matmulac2(T,nwvd(:,0),ntot1)
       ay=-dcmplx(0.d0,1.d0)*(matmulac2(T,phi,ntot1)+vks*phi+vin*phi)
    elseif(nd==3) then
       Tnwvd0=matmulac(T,nwvd(:,0),ntot1)
       ay=-dcmplx(0.d0,1.d0)*(matmulac(T,phi,ntot1)+vks*phi+vin*phi)
    else
       print*,'not implemented yet, stopping'
       stop
    end if
    deallocate(reord)
    nwvd(:,1)=dvks*phi
    nwvd(:,1)=vks*ay+nwvd(:,1) !Nprime

    f=dt**2*dcmplx(1.d0,0.d0)/12.d0*(dcmplx(-1.d0,0.d0)*&
         (Tnwvd0+vin*nwvd(:,0))+dcmplx(0.d0,1.d0)*&
         nwvd(:,1))
  end function calcf1c

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
