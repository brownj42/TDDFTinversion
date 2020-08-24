Module readin_data_mod
use derivedtypes 

contains
  subroutine readparameters(readin,filename,sysparams)
    use keomod
    use nchoosekmod
    logical, intent(in)                   :: readin
    character(50), intent(in), optional   :: filename
    type(systemparameters), intent(inout) :: sysparams
    real(8) :: xmax,xmin
    integer :: np1,ntot1,npart
    integer :: jin,prev

    if (readin) then
       if (present(filename)) then
          open(unit=33,file=filename)
          !number of dimensions for each particle
          read(33,*) sysparams%nd
          !number of grid points for each particle in 1d
          read(33,*) sysparams%np1
          !number of particles
          read(33,*) sysparams%npart
          !number of grid points for each single particle
          sysparams%ntot1=sysparams%np1**(sysparams%nd)
          !min x for box
          read(33,*) sysparams%xmin
          !max x for box
          read(33,*) sysparams%xmax
          !grid spacing
          sysparams%dx=(sysparams%xmax-sysparams%xmin)/dble(sysparams%np1-1)
          !time step
          read(33,*) sysparams%dt
          sysparams%dth=sysparams%dt
          !choose pseudo inverse or MINRES-QLP
          read(33,*) sysparams%pinv0minresqlp1
          close(33)
          sysparams%quantization=1
          sysparams%dvks_sw=1
          sysparams%dvksh_sw=1
          
       else
          print*,'Must input filename of parameters file'
          print*,'Stopping'
          stop
       end if
       !build KEO and lattice
       np1=sysparams%np1
       xmin=sysparams%xmin
       xmax=sysparams%xmax
       sysparams%newfile=1
       sysparams%outputdata=1
       allocate(sysparams%T(np1,np1))
       allocate(sysparams%xlattice(np1))
       call buildkeo(sysparams)
       !sysparams%xlattice=buildxlattice(np1,xmin,xmax)
       ntot1=sysparams%ntot1
       npart=sysparams%npart
       !number of grid points in interaction potential
       sysparams%dvksmax=1000.d0 
       call fill_systemparameters(sysparams) 
      
    end if
  end subroutine readparameters
  
  subroutine alloc_arrays(sysparams)
     use derivedtypes 
     TYPE(systemparameters), INTENT(inout) :: sysparams
     integer :: np1
     np1=sysparams%np1
     ALLOCATE(sysparams%T(np1,np1),sysparams%xlattice(np1))
  end subroutine alloc_arrays
end module readin_data_mod
