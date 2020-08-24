module density_mod
  implicit none
contains
  

  function fullwf_densitymatrix(ntot1,npart,yh) result(density)
    integer, intent(in) :: ntot1,npart
    complex(8), intent(in) :: yh(:)
    real(8) :: density(ntot1,ntot1)
    integer :: i,j,k
    density =0.d0
    if (npart==2) then
       do i=1,ntot1
          do k=1,ntot1
             do j=1,ntot1
                density(i,k)=density(i,k)+dble(2.d0*dconjg(yh((i-1)*ntot1+j))*yh((k-1)*ntot1+j))
             end do
          end do
       end do
    end if
  end function fullwf_densitymatrix

  function ks_density(ntot1,npart,y) result(density)
    integer,intent(in) :: ntot1,npart
    complex(8),intent(in) :: y(ntot1,npart)
    real(8) :: density(ntot1)
    integer :: i
    density =0.d0
    do i=1,npart
       density=density+dble(dconjg(y(:,i))*y(:,i))
    end do
  end function ks_density

 
  subroutine calcddnx3(N,n3,T,yex,vex,ddnx)
    use matmul_mod
    integer, intent(in) :: N,n3
    real(8), intent(in) :: T(N,N)
    complex(8), intent(in) :: yex(n3)
    real(8), intent(in) :: vex(n3)
    real(8), intent(out) :: ddnx(n3)
    complex(8), allocatable :: hy(:),hhy(:)
    roswitch=0
    allocate(reord(n3,3))
    allocate(hy(n3),hhy(n3))
    hy=matmulac(T,yex,n3)+vex*yex
    hhy=matmulac(T,hy,n3)+vex*hy
    ddnx=-1.d0*dble(hhy*dconjg(yex)+dconjg(hhy)*yex)+&
         2.d0*(dble(hy)**2+dimag(hy)**2)
    deallocatE(hy,hhy)
    deallocate(reord)
    
  end subroutine calcddnx3
  
  
  subroutine calcddnx1(N,T,yex,vex,ddnx)
    use matmul_mod
    integer, intent(in) :: N
    real(8), intent(in) :: T(N,N)
    complex(8), intent(in) :: yex(n)
    real(8), intent(in) :: vex(n)
    real(8), intent(out) :: ddnx(n)
    complex(8), allocatable :: hy(:),hhy(:)
    allocate(hy(n),hhy(n))
    hy=matmul(T,yex)+vex*yex
    hhy=matmul(T,hy)+vex*hy
    ddnx=-1.d0*dble(hhy*dconjg(yex)+dconjg(hhy)*yex)+2.d0*(dble(hy)**2+dimag(hy)**2)
  end subroutine calcddnx1
  
  
  !calculate 2nd derivative of density
  subroutine calcddnx2(N,n2,T,yex,vex,ddnx)
    use matmul_mod
    integer, intent(in) :: N,n2
    real(8), intent(in) :: T(N,N)
    complex(8), intent(in) :: yex(n2)
    real(8), intent(in) :: vex(n2)
    real(8), intent(out) :: ddnx(n2)
    complex(8), allocatable :: hy(:),hhy(:)
    !allocate index tracking vectors
    roswitch=0
    allocate(reord(n2,3))
    allocate(hy(n2),hhy(n2))
    !calculate diagonal elements of 2nd derivative of full density
    hy=matmulac2(T,yex,n2)+vex*yex
    hhy=matmulac2(T,hy,n2)+vex*hy
    ddnx=-1.d0*dble(hhy*dconjg(yex)+dconjg(hhy)*yex)+2.d0*(dble(hy)**2+dimag(hy)**2)
    deallocate(reord,hy,hhy)
  end subroutine calcddnx2
 
  !calculate 1st derivative of density 
  subroutine calcdnx2(N,n2,T,yex,vex,dnx)
    use matmul_mod
    integer, intent(in) :: N,n2
    real(8), intent(in) :: T(N,N)
    complex(8), intent(in) :: yex(n2)
    real(8), intent(in) :: vex(n2)
    real(8), intent(out) :: dnx(n2)
    complex(8), allocatable :: hy(:)
    !allocate index tracking vectors
    roswitch=0
    allocate(reord(n2,3))
    allocate(hy(n2))
    !calculate diagonal elements of 1st derivative of full density
    hy=-dcmplx(0.d0,1.d0)*(matmulac2(T,yex,n2)+vex*yex)
    dnx=dble(hy*dconjg(yex)+dconjg(hy)*yex)
    
    deallocate(reord,hy)
  end subroutine calcdnx2
  !calculate 1st derivative of density 
  subroutine calcdnx1(N,T,yex,vex,dnx)
    !use matmul_mod
    integer, intent(in) :: N
    real(8), intent(in) :: T(N,N)
    complex(8), intent(in) :: yex(n)
    real(8), intent(in) :: vex(n)
    real(8), intent(out) :: dnx(n)
    complex(8), allocatable :: hy(:)
    !allocate index tracking vectors
    allocate(hy(n))
    !calculate diagonal elements of 1st derivative of full density
    hy=-dcmplx(0.d0,1.d0)*(matmul(T,yex)+vex*yex)
    dnx=dble(hy*dconjg(yex)+dconjg(hy)*yex)
  end subroutine calcdnx1

  
  subroutine calcdnx3(N,n3,T,yex,dnx)
    use matmul_mod
    integer, intent(in) :: N,n3
    real(8), intent(in) :: T(N,N)
    complex(8), intent(in) :: yex(n3)
    real(8), intent(out) :: dnx(n3)
    complex(8), allocatable :: hy(:)
    roswitch=0
    allocate(reord(n3,3))
    allocate(hy(n3))
    hy=-dcmplx(0.d0,1.d0)*(matmulac(T,yex,n3))
    dnx=dble(hy*dconjg(yex)+dconjg(hy)*yex)
    deallocate(reord)
  end subroutine calcdnx3
   
end module density_mod
