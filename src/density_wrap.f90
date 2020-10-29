module density
  use density_mod
  implicit none
contains
  subroutine fullwf_density(sysparams,yh,density)
    use derivedtypes
    use secondquant_wrap 
    type(systemparameters), intent(in) :: sysparams
    complex(8), intent(in) :: yh(:)
    real(8), intent(out) :: density(sysparams%ntot1)
    integer :: ntot1,npart
    integer :: i,j
    ntot1=sysparams%ntot1
    npart=sysparams%npart
    if (sysparams%quantization==1) then
       density =0.d0
       if (npart==2) then
          do i=1,ntot1
             do j=1,ntot1
                density(i)=density(i)+dble(2.d0*dconjg(yh((i-1)*ntot1+j))*yh((i-1)*ntot1+j))
             end do
          end do
       else
          do i=1,ntot1
             density(i)=abs(dconjg(yh(i))*yh(i))
          end do
       end if
    else
       call diagonerdm(sysparams,sysparams%ntot1,yh,density)
    end if
  end subroutine fullwf_density
  
  
  subroutine calcddnx(sysparams,sharedvals,ntot1,psi,v,ddnx) 
    use derivedtypes 
    use secondquant_wrap
    type(systemparameters), intent(in) :: sysparams
    type(sharedvalues), intent(in) :: sharedvals 
    integer, intent(in) :: ntot1
    integer :: np1,nd,npart,ntot
    real(8), allocatable :: T(:,:)
    real(8), intent(in) :: v(:)
    complex(8), intent(in) :: psi(:)
    real(8), intent(out) :: ddnx(sysparams%ntot1)
    real(8), allocatable :: ddnx2(:)
    integer :: i
    np1=sysparams%np1
    nd=sysparams%nd
    npart=sysparams%npart
    ntot=sysparams%ntot
    if (sysparams%quantization==1) then
       allocate(T(np1,np1))
       T=sysparams%T
      
       if (npart==2.and.nd==1) then
          allocate(ddnx2(ntot))
          call calcddnx2(Np1,ntot,T,psi,v,ddnx2)
          !calculate reduced diagonal of 2nd derivative of density
          do i=1,np1
             ddnx(i)=2.d0*sum(ddnx2((i-1)*np1+1:i*np1))
          end do
          deallocate(ddnx2)
       elseif (npart==1.and.nd==1) then
          call calcddnx1(np1,t,psi,v,ddnx)
       elseif (npart==1.and.nd==3) then
          call calcddnx3(Np1,ntot1,T,psi,v,ddnx)
       end if
    else
       call diagddonerdm(sysparams,sharedvals,sysparams%ntot1,psi,ddnx)
    end if 
  end subroutine calcddnx
  subroutine calcdnx(sysparams,sharedvals,ntot1,psi,v,dnx)
    use derivedtypes 
    use secondquant_wrap
    type(systemparameters), intent(in) :: sysparams
    type(sharedvalues), intent(in) :: sharedvals
    integer, intent(in) :: ntot1
    integer :: np1,nd,npart,ntot
    real(8), allocatable :: T(:,:)
    real(8),intent(in) :: v(:)
    complex(8),intent(in) :: psi(:)
    real(8),intent(out) :: dnx(ntot1)
    real(8), allocatable :: dnx2(:)
    integer :: i
    
    nd=sysparams%nd
    np1=sysparams%np1
    npart=sysparams%npart
    ntot=sysparams%ntot
   
    
    allocate(T(np1,np1))
    T=sysparams%T
    if (sysparams%quantization==1) then
       if (nd==1.and.npart==2) then
          allocate(dnx2(ntot))
          call calcdnx2(np1,ntot,T,psi,v,dnx2)
          !calculate diagonal of reduced 1st derivative of density
          do i=1,np1
             
             dnx(i)=2.d0*sum(dnx2((i-1)*np1+1:i*np1))
          end do
          deallocate(dnx2)
       elseif (nd==1.and.npart==1) then
          call calcdnx1(np1,T,psi,v,dnx)
       elseif (nd==3.and.npart==1) then
          call calcdnx3(np1,ntot,T,psi,dnx)
       end if
    else
       call diagdonerdm(sysparams,sharedvals,sysparams%ntot1,psi,dnx)
    end if
  end subroutine calcdnx
  
  subroutine calcenergy(sysparams,sharedvals,psi,v,energy)
    use matmul_mod
    use derivedtypes
    type(systemparameters), intent(in) :: sysparams
    type(sharedvalues), intent(in) :: sharedvals
    complex(8),intent(in) :: psi(:)
    real(8),intent(in) :: v(:)
    real(8), intent(out) :: energy
    integer :: np1,nd,npart,ntot
    real(8), allocatable :: T(:,:)
    complex(8), allocatable :: hpsi(:)
    real(8), allocatable :: dnx2(:)
    integer :: i
    include 'mkl_blas.fi'
    nd=sysparams%nd
    np1=sysparams%np1
    npart=sysparams%npart
    ntot=sysparams%ntot
    allocate(hpsi(ntot))
    
    
    
    allocate(T(np1,np1))
    T=sysparams%T
    if (sysparams%quantization==1) then
       if (nd*npart==1) then
          hpsi=matmul(T,psi)+v*psi
       elseif (nd*npart==2) then
          allocate(reord(ntot,3))
          roswitch=0
          hpsi=matmulac2(T,psi,ntot)+v*psi
          deallocate(reord)
       elseif (nd*npart==3) then
          allocate(reord(ntot,3))
          roswitch=0
          hpsi=matmulac(T,psi,ntot)+v*psi
          deallocate(reord)
       else
          print*,'choice of parameters is not implemented yet'
          print*,'stopping now'
          stop
       end if
    else
       print*,'stopping: not implemented yet'
    end if
    energy=zdotc(ntot,hpsi,1,psi,1)
  end subroutine calcenergy
  
 
   
end module density
