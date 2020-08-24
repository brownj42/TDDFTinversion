module outputdata_mod
  implicit none
contains
  subroutine printdata(np1,nd,ntot1,npart,xlattice,ct,vks,vint,v1,dp,dpe,phi,newfile)
    integer, intent(in) :: np1,nd,ntot1,npart
    real(8), intent(in) :: xlattice(np1),ct,vks(ntot1),vint(ntot1)
    real(8), intent(in) :: v1(ntot1)
    real(8), intent(in) :: dp(ntot1),dpe(ntot1)
    complex(8), intent(in) :: phi(ntot1,npart)
    integer, intent(inout) :: newfile
    complex(8) :: val
    integer :: i,p,nm
    logical :: exist
    
    inquire(file='times.dat',exist=exist)
    if (exist.and.newfile.ne.1) then
        open(unit=220,file='times.dat', status='old',position='append',action='write')
    else
        open(unit=220,file='times.dat',action='write')
    end if
    inquire(file='pots.dat',exist=exist)
    if (exist.and.newfile.ne.1) then
        open(unit=221,file='pots.dat', status='old',position='append',action='write')
    else
        open(unit=221,file='pots.dat',action='write')
    end if
    inquire(file='dense.dat',exist=exist)
    if (exist.and.newfile.ne.1) then
        open(unit=222,file='dense.dat', status='old',position='append',action='write')
    else
        open(unit=222,file='dense.dat',action='write')
    end if
    inquire(file='wave.dat',exist=exist)
    if (exist.and.newfile.ne.1) then
        open(unit=232,file='wave.dat', status='old',position='append',action='write')
    else
        open(unit=232,file='wave.dat',action='write')
    end if
    
    newfile=0    
    write(220,'(f20.10)') ct

    if (npart==2) then
       !do i=1,ntot1
       !   write(232,'(5f20.10)')xlattice(i),real(phi(i,1)),&
       !        aimag(phi(i,1)),real(phi(i,2)),aimag(phi(i,2))
       !end do
    end if
    if (nd==1) then
       do i=1,ntot1
          write(221,'(5f20.10)')xlattice(i),vks(i),vks(i)-vint(i)-v1(i),&
            vint(i),vks(i)-v1(i)
          write(222,*)xlattice(i),dpe(i),dp(i)
       end do
       do i=1,ntot1
          write(232,'(1f20.10)',advance='no') xlattice(i)
          do p=1,npart-1
             write(232,'(2f20.10)',advance='no')real(phi(i,p)),&
                  aimag(phi(i,p))
          end do
          p=npart
          write(232,'(2f20.10)')real(phi(i,p)),&
                  aimag(phi(i,p))
       end do 
    elseif(nd==3) then
       nm=(np1+1)/2
       do i=1,np1
          p=(nm-1)*np1**2+(nm-1)*np1
          write(221,'(5f20.10)')xlattice(i),vks(p+i),vks(p+i)-vint(p+i)-v1(p+i),&
            vint(p+i),vks(p+i)-v1(p+i)
          write(222,*)xlattice(i),dpe(p+i),dp(p+i)
       end do
       do i=1,np1
          val=sum(phi((i-1)*np1**2+1:(i-1)*np1**2+(np1-1)*np1+np1,1))
          write(232,'(3f20.10)')xlattice(i),real(val),aimag(val)
       end do
    end if
    write(221,*)' '
    write(221,*)' '
    write(222,*)' '
    write(222,*)' '
    write(232,*) ' '
    write(232,*) ' '
    close(220)
    close(221)
    close(222)
    close(232)
  end subroutine printdata
end module outputdata_mod
