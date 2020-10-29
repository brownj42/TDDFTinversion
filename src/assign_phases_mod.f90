Module assign_phases_mod
  use matmul_mod
  implicit none
    include 'mkl_blas.fi' 
Contains
  !function that assigns phases such that 1st derivative of density matches
  function assign_phases(psii,dn_aim,T,np1,nd,N,npart,tol) result(psif)
    use lsqrdatamodule
    !N=number of lattice points,npart is number of particles
    integer :: N,npart,np1,nd
    !Initial KS wavefunctions
    complex(8) :: psii(N,npart)
    !first derivative of density aim
    real(8) :: dn_aim(N)
    !2nd derivative matrix, tolerance of assigned phases
    real(8) :: T(np1,np1),tol
    integer :: ctr,j
    complex(8) :: j1
    real(8) :: new_error,abspsii(N,Npart)
    complex(8) :: psif(n,npart)
    real(8) :: new_phases(n*npart)
    
    j1=dcmplx(0.d0,1.d0)
    roswitch=0
    allocate(reord(n,3)) 
    !Transfer npart to global variable for LSQR
    npartw=npart
    print*,'Assigning phases'
    new_error=dnrm2(n,dn_aim-derivative_of_density(psii,T,np1,nd,n,npart),1)/dble(n)
    print*,'Initial phases error',new_error
    do j=1,npart
       !psif(:,j)=psii(:,j)
       call zcopy(n,psii(:,j),1,psif(:,j),1)
    end do
    !take newton steps to assign phases 
    if (new_error.gt.1.d-14) then
       do j=1,npart
          !abspsii(:,j)=abs(psii(:,j))
          abspsii(:,j)=absvec(n,psii(:,j))
       end do
       newton:do ctr=1,14
          new_phases=newton_step(psif,dn_aim,T,np1,nd,n,npart)
          do j=1,npart
             psif(:,j)=abspsii(:,j)*zexp(j1*new_phases((j-1)*n+1:j*n))
          end do
          new_error=dnrm2(n,dn_aim-derivative_of_density(psif,T,np1,nd,n,npart),1)/dble(n)
          if (new_error.lt.tol) then
             exit newton
          end if
       end do newton
    end if
    print*,'Final phases error',new_error
    deallocate(reord)
  end function assign_phases
 
  !function returns abs values of double complex vector
  function absvec(n,vec) result(avec)
    integer :: n
    complex(8) :: vec(n)
    real(8) :: avec(n)
    integer :: i
    !OMP Parallel shared(avec,vec)
    !OMP DO schedule(static,4) firstprivate(n)
    do i=1,n
       avec(i)=abs(vec(i))
    end do
    !OMP end do nowait
    !OMP end parallel
  end function absvec

  !calculates L2-norm of a vector
  real(8) function norm(v)
    include 'mkl_blas.fi' 
    complex(8) :: v(:)
    integer :: ntot
    ntot=size(v)
    norm=dsqrt(abs(zdotc(ntot,v,1,v,1)))
  end function norm
  real(8) function normd(v)
    include 'mkl_blas.fi' 
    real(8) :: v(:)
    integer :: ntot
    ntot=size(v)
    normd=dsqrt(ddot(ntot,v,1,v,1))
  end function normd

  !Calculates phase angle for each component of complex vectos
  function angle(v,n,npart) result(ph)
    !length of vector, number of vectors (particles)
    integer :: n,npart
    !matrix of vectors
    complex(8) :: v(n,npart)
    real(8) :: ph(n*npart)
    integer :: i,p,ps
    !stack in single vector
    !calculate phase angles
   
    do p=1,npart
       ps=(p-1)*n
       !OMP parallel shared(ph,v) firstprivate(npart,p,n,ps) default(private) 
       !OMP do schedule(static,4)
       do i=1,n
          ph(ps+i)=datan2(dimag(v(i,p)),dble(v(i,p)))
       end do
       !OMP end do nowait
       !OMP end parallel
    end do
  end function angle
  
  !function that takes step in direction of improved phases
  function newton_step(y,dn_aim,T,np1,nd,n,npart) result(phases)
    !use derivedtypes
    use lsqrdatamodule
    use pot_invert_mod 
    use lsqrmodule
    use minresqlpmodule
    !size of total grid, number of particles
    integer :: np1,nd,n,npart 
    !KS wavefunctions
    complex(8) :: y(n,npart)
    !target first derivative of density
    real(8) :: dn_aim(n)
    !2nd spatial derivative matrix
    real(8) :: T(np1,np1),Tdiag
    real(8), allocatable :: mT(:,:)
    real(8) :: phases(n*npart)
    real(8) :: gp(n)
    real(8) :: epsmaxit,epsepsmaxit
    real(8) :: step(n*npart)
    real(8) :: rgp(n),ry(n,npart),iy(n,npart),dy(n,npart),Mm1(n*npart),maxit,mms(n)
    real(dp) :: se(n)
    real(8), allocatable :: G(:,:),Km(:,:)
    complex(8) :: ty(n)
    integer :: sn3,p,ps,i,nnp
    integer :: istop, itn,np,loc(1,npart)
    real(dp) ::  Anorm, Acond, rnorm, Arnorm, xnorm,phasefix(npart),phaseshift
    sn3=size(T,1)
    allocate(mT(sn3,sn3))
    mT=-T
    !calculate phases of KS system wavefunctions
    phases=angle(y,n,npart)
    !fix phases at highest density point for each KS state
    do i=1,npart
       !loc(:,i)=maxloc(abs(y(:,i)))
       loc(1,i)=izamax(n,y(:,i),1)
       phasefix(i)=phases((i-1)*n+loc(1,i))
    end do
    
    !real and imaginary parts for each vector
    ry=dble(y)
    iy=dimag(y)
    nnp=n*npart
    !obtain right hand side vector
    gp=derivative_of_density(y,T,np1,nd,n,npart)!-dn_aim
    call daxpy(n,-1.d0,dn_aim,1,gp,1) 
    !rgp=gp
    call dcopy(n,gp,1,rgp,1)
    !obtain first derivative of density for each particle in KS system
    do i=1,npart
       if (nd==1) then
           ty=matmul(mT,y(:,i))
       elseif (nd==3) then
           ty=matmulac(mT,y(:,i),n)
       end if
       !OMP parallel shared(ty,y,dy) firstprivate(n)
       !OMP do schedule(static,4)
       do p=1,n
          dy(p,i)=-1.d0*dble(dconjg(ty(p))*(y(p,i))+dconjg(y(p,i))*ty(p))
       end do
       !OMP end do nowait
       !OMP end parallel
    end do

    !obtain preconditioners
    p=0
    Tdiag=T(1,1)
    do np=1,npart
       ps=(np-1)*n
       !OMP parallel shared(mm1,ry,iy,dy) firstprivate(n,nd,Tdiag,np,ps)
       !OMP do schedule(static,4)
       do i=1,n
          mm1(ps+i)=-2.d0*(ry(i,np)**2+iy(i,np)**2)*(nd)-dy(i,np)
       end do
       !OMP end do nowait
       !OMP end parallel
    end do
    !OMP parrallel shared(mms) firstprivate(n)
    !OMP do
    do i=1,n
    mms(i)=0.d0
    end do
    !OMP end do
    !OMP end parallel
    do np=1,npart
       !OMP parallel shared(mms,ry,iy,dy) firstprivate(n,nd,Tdiag,np)
       !OMP do schedule(static,4)
       do i=1,n
          mms(i)=mms(i)-2.d0*(ry(i,np)**2+iy(i,np)**2)*(nd)!*Tdiag)
          mms(i)=mms(i)-dy(i,np)
       end do
       !OMP end do nowait
       !OMP end parallel
    end do
    !maxit=maxval(dabs(mm1))
    p=idamax(nnp,mm1,1)
    maxit=dabs(mm1(p))
    p=0
    epsmaxit=epsilon(maxit)
    epsepsmaxit=dsqrt(1.d0*dabs(epsilon(1.d0/epsmaxit)))
    !OMP parallel shared(mm1) firstprivate(epsepsmaxit,epsmaxit)
    !OMP do schedule(static,4)
    do p=1,nnp
       if (abs(mm1(p)).gt.epsmaxit) then
          !mm1(p)=dsqrt(1.d0/(dabs(mm1(p))))
          mm1(p)=dsqrt(1.d0/(dabs(mm1(p))))
       else
          !mm1(p)=dsqrt(1.d0*dabs(epsilon(1.d0/epsilon(maxit))))
          mm1(p)=epsepsmaxit
       end if
    end do
    !OMP end do nowait
    !OMP end maxit
    !maxit=maxval(dabs(mms))
    p=idamax(n,mms,1)
    maxit=dabs(mms(p))
    epsmaxit=epsilon(maxit)
    epsepsmaxit=dsqrt(1.d0*dabs(epsilon(1.d0/epsmaxit)))
    !OMP parallel shared(mms) firstprivate(epsepsmaxit,epsmaxit)
    !OMP do schedule(static,4)
    do p=1,n
       if (abs(mms(p)).gt.epsmaxit) then
          mms(p)=dsqrt(1.d0/(dabs(mms(p))))
       else
          !mms(p)=sqrt(1.d0*dabs(epsilon(1.d0/epsilon(maxit))))
          mms(p)=epsepsmaxit
       end if
    end do
    !OMP end do nowait
    !OMP end maxit

    !use LSQR program to obtain new phases
    allocate(TW(sn3,sn3),ryw(n,npartw),iyw(n,npartw),dyw(n,npartw),mWs(n),mW(nnp))
    !mW=mm1
    call dcopy(nnp,mm1,1,mW,1)
    !mWs=mms
    call dcopy(n,mms,1,mWs,1)
    TW=mT
    do p=1,npartw
       !ryw=ry
       call dcopy(n,ry(:,p),1,ryw(:,p),1)
       !iyw=iy
       call dcopy(n,iy(:,p),1,iyw(:,p),1)
       !dyw=dy 
       call dcopy(n,dy(:,p),1,dyw(:,p),1)
    end do
    !rgp=rgp*mWs
    rgp=vtimesv(rgp,mWs,n)
    open(111, file='LSQR.txt')
    if (nd==1) then
       
       if (npart==1) then
          allocate(G(n,n),km(n,n))
          G=0.d0
          !print*,'balyb'
          G=2.d0*dble(matmul(y(:,1:1),dconjg(transpose(y(:,1:1)))))
          Km=-G*mT
          G=0.d0
          do i=1,n
             Km(i,i)=Km(i,i)-dy(i,1)
             G(i,i)=mws(i)
          end do
          G=matmul(G,matmul(Km,G))
          !print*,'balyb'
          call pinv(G,n,n,1.d-15,km)
          step=matmul(km,rgp)
          deallocate(G,km)
       else
          call LSQR(n,nnp,Aprodw1, Aprodw2, rgp, 0.d0, .false.,   &
                     step, se,                                    &
                     1.d-15,1.d-15, 1.d5, 2*nnp,111,              &
                     istop, itn, Anorm, Acond, rnorm, Arnorm, xnorm )
       end if
    else if (nd==3) then
       call LSQR(n,nnp,Aprodw1_3d, Aprodw2_3d, rgp, 0.d0, .false.,       &
                     step, se,                                         &
                     1.d-15,1.d-15, 1.d5, 40,111,              &
                     istop, itn, Anorm, Acond, rnorm, Arnorm, xnorm )

    end if
    close(111)
    step=vtimesv(step,mw,nnp)
    if (normd(step).gt.1d-2) then
            step=step/normd(step)*1d-2
    end if
    !phases=phases-step*0.5d0
    call daxpy(nnp,-1.d0,step,1,phases,1)
    do p=1,npart
       ps=(p-1)*n
       phaseshift=phases(ps+loc(1,p))-phasefix(p)
       !$OMP Parallel shared(phases) firstprivate(phaseshift,ps,n)
       !$OMP DO schedule(static,4)
       do i=1,n
          phases(ps+i)=phases(ps+i)-phaseshift
       end do
       !$OMP END do nowait
       !$OMP end parallel
    end do
    deallocate(mT)
    deallocate(TW,ryw,iyw,dyw,mW,mws)
    
  end function newton_step

  subroutine MsolveW_3d(n,x,y)                   ! Solve M*y = x
    use lsqrDataModule
    integer, intent(in)    :: n
    real(dp),    intent(in)    :: x(n)
    real(dp),    intent(out)   :: y(n)
    !y=Mw*x
    y=vtimesv(Mw,x,n)
  end subroutine Msolvew_3d

  !apply K for each particle to obtain N x Npart length vector
  subroutine AprodW1(n,nnp,x,y)
    use lsqrDataModule
    integer, intent(in)  :: n,nnp
    real(dp),    intent(in)  :: x(nnp)
    real(dp),    intent(inout) :: y(n)
    real(dp) :: xp(nnp)
    integer :: j,js,jf
    xp=mW*x

    do j=1,npartw
       js=(j-1)*n+1
       jf=j*n
       y=y-mws*2.d0*ryW(:,j)*matmul(TW,xp(js:jf)*ryW(:,j))
       y=y-mws*2.d0*iyW(:,j)*matmul(TW,xp(js:jf)*iyW(:,j))
       y=y-mws*dyW(:,j)*xp((j-1)*n+1:j*n)
    end do
  end subroutine AprodW1
  subroutine AprodW1_3d(n,nnp,x,y)
    use lsqrDataModule
    integer, intent(in)  :: n,nnp
    real(dp),    intent(in)  :: x(nnp)
    real(dp),    intent(inout) :: y(n)
    real(dp) :: xp(nnp),vec(n),vech(n)
    integer :: j,js,jf
    !xp=mW*x
    xp=vtimesv(mW,x,nnp)

    do j=1,npartw
       js=(j-1)*n+1
       jf=j*n
       !y=y-2.d0*ryW(:,j)*matmula(TW,xp(js:jf)*ryW(:,j),n)
       !y=y-2.d0*iyW(:,j)*matmula(TW,xp(js:jf)*iyW(:,j),n)
       !y=y-dyW(:,j)*xp((j-1)*n+1:j*n)
       vec=vtimesv(ryW(:,j),matmula(TW,vtimesv(xp(js:jf),ryW(:,j),n),n),n)
       vech=vtimesv(iyW(:,j),matmula(TW,vtimesv(xp(js:jf),iyW(:,j),n),n),n)
       call daxpy(n,1.d0,vech,1,vec,1)
       vech=vtimesv(dyW(:,j),xp(js:jf),n)
       call daxpy(n,0.5d0,vech,1,vec,1)
       vech=vtimesv(mws,vec,n)
       call daxpy(n,-2.d0,vech,1,y,1)
    end do
  end subroutine AprodW1_3d


  !apply K for each particle to length N vector to obtain length  N vector
  subroutine AprodW2(n,nnp,x,y)
    use lsqrDataModule
    integer, intent(in)  :: n,nnp
    real(dp),    intent(inout)  :: x(nnp)
    real(dp),    intent(in) :: y(n)
    real(dp) :: xp(nnp)
    integer :: j,js,jf
    xp=x
    x=0.d0
    do j=1,npartw
       js=(j-1)*n+1
       jf=j*n
       x(js:jf)=x(js:jf)-2.d0*ryW(:,j)*matmul(TW,mws*y*ryW(:,j))
       x(js:jf)=x(js:jf)-2.d0*iyW(:,j)*matmul(TW,mws*y*iyW(:,j))
       x(js:jf)=x(js:jf)-dyW(:,j)*y*mws
    end do
    x=xp+mw*x
  end subroutine AprodW2
  subroutine AprodW2_3d(n,nnp,x,y)
    use lsqrDataModule
    integer, intent(in)  :: n,nnp
    real(dp),    intent(inout)  :: x(nnp)
    real(dp),    intent(in) :: y(n)
    real(dp) :: xp(nnp),yp(n)
    integer :: j,js,jf
    yp=vtimesv(mws,y,n)
    !xp=x
    call dcopy(n,x,1,xp,1)
    x=0.d0
    do j=1,npartw
       js=(j-1)*n+1
       jf=j*n
       !x(js:jf)=x(js:jf)-2.d0*vtimesv(ryW(:,j),matmula(TW,vtimesv(y,ryW(:,j),n),n),n)
       !x(js:jf)=x(js:jf)-2.d0*vtimesv(iyW(:,j),matmula(TW,vtimesv(y,iyW(:,j),n),n),n)
       !x(js:jf)=x(js:jf)-vtimesv(dyW(:,j),y,n)
       call daxpy(n,-2.d0,vtimesv(ryW(:,j),matmula(TW,vtimesv(yp,ryW(:,j),n),n),n),1,x(js:jf),1)
       call daxpy(n,-2.d0,vtimesv(iyW(:,j),matmula(TW,vtimesv(yp,iyW(:,j),n),n),n),1,x(js:jf),1)
       call daxpy(n,-1.d0,vtimesv(dyW(:,j),yp,n),1,x(js:jf),1)
    end do
    x=xp+vtimesv(mw,x,nnp)
  end subroutine AprodW2_3d


  !calculates derivative of density of KS system y with n3 points and npart
  !particles, T is 2nd spatial derivative matrix.
  function derivative_of_density(y,T,np1,nd,n3,npart) result(yout)
    use pot_invert_mod 
    integer :: n3,npart,np1,nd
    complex(8) :: y(n3,npart)
    real(8) :: T(np1,np1)
    real(8) :: yout(n3)
    integer :: i
    yout=0.d0
    do i=1,npart
       if (nd==1) then
          yout=yout+dble(dcmplx(0.d0,1.d0)*(y(:,i)*matmul(T,dconjg(y(:,i)))-&
            matmul(T,y(:,i))*dconjg(y(:,i))))
       elseif(nd==3) then
          yout=yout+dble(dcmplx(0.d0,1.d0)*(y(:,i)*matmulac(T,dconjg(y(:,i)),n3)-&
            matmulac(T,y(:,i),n3)*dconjg(y(:,i))))
       end if
    end do
  end function derivative_of_density
  
end Module assign_phases_mod
