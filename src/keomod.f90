module keomod
  implicit none
contains
  !returns sysparams%T with sinc DVR KEO
  subroutine buildkeo(sysparams)
    use derivedtypes 
    real(8) :: pi,bandlimit
    type(systemparameters) ,intent(inout) :: sysparams
    integer :: i,j
    if (sysparams%quantization==1) then
       pi=4.d0*datan(1.d0)
       bandlimit=dble(sysparams%np1-1)*pi/(sysparams%xmax-sysparams%xmin)
       do i=1,sysparams%np1
          do j=1,sysparams%np1
             if (i==j) then
                sysparams%T(i,j)=bandlimit**2/6.d0
             else
                sysparams%T(i,j)=(bandlimit/pi)**2*(-1.d0)**(j-i)/dble(j-i)**2
             end if
          end do
       end do
    else
       sysparams%T=0.d0
       do i=1,sysparams%np1
          do j=1,sysparams%np1
             if (i==j+1) then
                sysparams%T(i,j)=-1.0
             elseif (i+1==j) then
                sysparams%T(i,j)=-1.0
             end if
          end do
          sysparams%T(i,i)=0.d0
       end do
       sysparams%T(1,sysparams%np1)=-0.d0
       sysparams%T(sysparams%np1,1)=-0.d0
    end if
              
                
    call buildxlattice(sysparams%np1,sysparams%xmin,sysparams%xmax,sysparams%xlattice)
  end subroutine buildkeo

  subroutine buildxlattice(np1,xmin,xmax,xlattice)
    integer, intent(in)  :: np1
    real(8), intent(in)  :: xmin,xmax
    real(8) :: dx
    real(8), intent(inout) :: xlattice(np1)
    integer :: i
    dx=(xmax-xmin)/dble(np1-1)
    do i=1,np1
       xlattice(i)=(i-1)*dx+xmin
    end do
  end subroutine buildxlattice
end module keomod
