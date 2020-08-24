
  !This is implicitly restarted preconditioned GMRES method
  !Does not work very well to invert the potential
  subroutine gmres3(s3n,n,A,Mm1,ry,iy,dy,b,max_iterations,threshold,x,maxi)
    integer :: s3n,n
    real(8), intent(in) :: A(s3n,s3n),b(n),ry(n),iy(n),dy(n),threshold,Mm1(n)
    integer, intent(in) :: max_iterations,maxi
    real(8), intent(inout) :: x(n)
    integer :: m
    real(8) :: beta
    real(8), allocatable :: r(:)
    real(8), allocatable :: e(:),Q(:,:),H(:,:),hh(:,:),em(:,:),ihh(:,:),hharm(:,:)
    real(8) :: r_norm,b_norm,error
    integer, allocatable :: ipos(:)
    real(8), allocatable :: Gin(:,:),P(:,:),c(:),d(:),vadd(:),mQ(:,:),hhold(:,:)
    real(8), allocatable :: y(:),yin(:,:)
    integer :: info,k,nk,kstart,km,i,j,kmh
    Character :: Jobvl,Jobvr
    integer :: LWORK,blah
    real(8) :: qwork(1)
    real(8), allocatable :: work(:),wr(:),wi(:),Vr(:,:),VL(:,:)
    logical, allocatable :: bwork(:)
    integer :: sdim
    kstart=1
    km=2
    kmh=2
    m=max_iterations
    nk=m-1

    
    allocate(r(n))
    allocate(Q(n,m),H(m,nk),hhold(m,nk),e(m),yin(m,1),em(m-1,1))
    allocate(ihh(nk,nk),hh(nk,nk))
    allocate(y(nk),c(m),d(nk),vadd(m))
    allocate(bwork(nk))
    Jobvl='N'
    Jobvr='V'
    allocate(wr(nk),wi(nk),vl(nk,nk),vr(nk,nk))
    allocate(hharm(nk,nk),ipos(nk))
    
    
    !r=b-matmulW(A,ry,iy,dy,x,n)
    call aprodw(n,x,r)
    r=b-r
    
    !pause
    b_norm=sqrt(dot_product(b,b))
    r_norm=sqrt(dot_product(r,r))
    error=r_norm/b_norm
    Q=0.d0
    Q(:,1)=r/r_norm
    beta=r_norm
    c=0.d0
    c(1)=beta
    
    H=0.d0
    if (b_norm.gt.1.d-14) then
       kstart=1
       do k=kstart,m-1
          call arnoldi(A,Mm1,ry,iy,dy,Q,k,H(1:k+1,k))
       end do
       !c=e
       y=matmul(transpose(H(1:nk+1,1:nk)),c(1:nk+1))
       hh=matmul(transpose(H(1:nk+1,1:nk)),H(1:nk+1,1:nk))
       call pinv(hh(1:nk,1:nk),nk,nk,threshold,ihh)
       d=matmul(ihh,y(1:nk))
       !hhold=h
       !yin(:,1)=c
       !lwork=-1
       !call dgels('N',m,nk,1,hhold,m,yin,m,qwork,lwork,info)
       !lwork=nint(qwork(1))
       !allocate(work(lwork))
       !call dgels('N',m,nk,1,hhold,m,yin,m,work,lwork,info)
       !deallocate(work)
       !d=yin(1:nk,1)
       !yin(:,1)=c
       !print*,matmul(H,y)-c
       x=x+Mm1*matmul(Q(:,1:nk),d(1:nk))
       !r=b-matmulW(A,ry,iy,dy,x,n)
       call aprodw(n,x,r)
       r=b-r
       b_norm=sqrt(dot_product(b,b))
       r_norm=sqrt(dot_product(r,r)) 
       !print*,'after',r_norm
       if (r_norm.lt.1.d-13) then
          return
       end if
       
       hh=transpose(H(1:nk,1:nk))
       call pinv(hh,nk,nk,threshold,ihh)
       !call pinv(hh,nk,nk,1.d-16,ihh)
       em=0.d0
       em(nk,1)=1.d0
       Hharm=H(1:nk,1:nk)+matmul(matmul(ihh,em),transpose(em))
       lwork=-1
       call dgeev(JOBVL,JOBVR,nk,hharm,nk,wr,wi,vl,nk,vr,nk,qwork,lwork,info)
       lwork=nint(qwork(1))
       allocate(work(lwork))
       call dgeev(JOBVL,JOBVR,nk,hharm,nk,wr,wi,vl,nk,vr,nk,work,lwork,info)
       deallocate(work)
       do i=1,nk
          e(i)=sqrt(wr(i)**2+wi(i)**2)
       end do
       call IIInsertionSort(e(1:nk),ipos,nk)
       !print*,'ebot',wr(ipos(1))
       !print*,'etop',wr(ipos(km)),wi(ipos(km+1))
       if (dabs(wi(ipos(km))).gt.1.d-15) then
          if (dabs(wi(ipos(km+1))+wi(ipos(km))).lt.1.d-15) then
             km=km+1
          end if
       end if
       allocate(wrh(km))
       kpos=km
       do i=1,km
           wrh(i)=wr(ipos(i))
       end do
       lwork=-1
       Hharm=H(1:nk,1:nk)+matmul(matmul(ihh,em),transpose(em))
       call dgees(JOBVR,'S',selectr,nk,hharm,nk,sdim,wr,wi,vr,nk,qwork,lwork,bwork,info)
       lwork=nint(qwork(1))
       allocate(work(lwork))
       call dgees(JOBVR,'S',selectr,nk,hharm,nk,sdim,wr,wi,vr,nk,work,lwork,bwork,info)
       deallocate(work)
       deallocate(wrh)
       scoob:do blah=1,maxi
          allocate(Gin(nk,km))
          do i=1,km
             !print*,e(i),wr(ipos(i)),wi(ipos(i))
             !Gin(:,i)=VR(:,ipos(i))
             !print*,e(i),sqrt(wr(i)**2+wi(i)**2)
             Gin(:,i)=VR(:,i)
          end do
          !print*,'here'
          !call orth(Gin(1:nk,1:km),nk,km)
          allocate(P(m,km+1))
          P(1:nk,1:km)=Gin(1:nk,1:km)
          deallocate(Gin)
          P(m,1:km)=0.d0
          !call orth(P(:,1:km),m,km)
          vadd(1:m)=c(1:m)-matmul(h(1:m,1:nk),d(1:nk))
          vadd=vadd/sqrt(dot_product(vadd,vadd))
          do j=1,km
             vadd=vadd-dot_product(vadd,P(:,j))/dot_product(P(:,j),P(:,j))*P(:,j)
          end do
          P(:,km+1)=vadd/sqrt(dot_product(vadd,vadd))
          !print*,matmul(transpose(P),P)
          
          !call orth(P,m,km+1)
          
          H(1:km+1,1:km)=matmul(matmul(transpose(P),H),P(1:nk,1:km))
          H(km+2:m,:)=0.d0
          H(:,km+1:nk)=0.d0
          Q(:,1:km+1)=matmul(Q,P)
          Q(:,km+2:m)=0.d0
          allocate(mQ(n,km))
          mq=q(:,1:km)
          call orth(mq,n,km)
          do j=1,km
             Q(:,km+1)=Q(:,km+1)-dot_product(Q(:,km+1),mQ(:,j))/dot_product(mQ(:,j),mQ(:,j))*mQ(:,j)
          end do
          !Q(:,km+1)=Q(:,km+1)/sqrt(dot_product(Q(:,km+1),Q(:,km+1)))
          !print*,matmul(transpose(Q(:,1:km+1)),Q(:,km+1))
          do j=1,km
             mq(:,j)=mm1*q(:,j)
             !print*,maxval(abs(matmulw(A,ry,iy,dy,mQ(:,j),n)-matmul(Q(:,1:km+1),H(1:km+1,j))))
          end do
          deallocate(mq)
          deallocate(P)
          !print*,matmul(transpose(Q(:,1:km+1)),Q(:,1:km+1))
          !call orth(Q(:,1:km+1),n,km+1)
          kstart=km+1
          do k=kstart,m-1
             call arnoldi(A,Mm1,ry,iy,dy,Q,k,H(1:k+1,k))
          end do
          km=kmh
          beta=H(m,nk)
          c=matmul(transpose(Q),r)
          y=matmul(transpose(H(1:nk+1,1:nk)),c(1:nk+1))
          hh=matmul(transpose(H(1:nk+1,1:nk)),H(1:nk+1,1:nk))
          !do k=1,nk
          !   hh(k,k)=hh(k,k)+1.d-3
          !   end do
          call pinv(hh(1:nk,1:nk),nk,nk,threshold,ihh)
          d=matmul(ihh,y(1:nk))
       !hhold=h
       !yin(:,1)=c
       !lwork=-1
       !call dgels('N',m,nk,1,hhold,m,yin,m,qwork,lwork,info)
       !lwork=nint(qwork(1))
       !allocate(work(lwork))
       !call dgels('N',m,nk,1,hhold,m,yin,m,work,lwork,info)
       !deallocate(work)
       !d=yin(1:nk,1)
          !print*,'norm d',sqrt(dot_product(d,d))
          x=x+Mm1*matmul(Q(:,1:nk),d(1:nk))
          !r=b-matmulW(A,ry,iy,dy,x,n)
          call aprodw(n,x,r)
          r=b-r
          b_norm=sqrt(dot_product(b,b))
          r_norm=sqrt(dot_product(r,r)) 
          !print*,'after',r_norm
          if (r_norm.lt.1.d-13) then
             exit scoob
          end if
          hh=transpose(H(1:nk,1:nk))
          call pinv(hh,nk,nk,threshold,ihh)
          em=0.d0
          em(nk,1)=1.d0
          Hharm=H(1:nk,1:nk)+matmul(matmul(ihh,em),transpose(em))
          lwork=-1
          call dgeev(JOBVL,JOBVR,nk,hharm,nk,wr,wi,vl,nk,vr,nk,qwork,lwork,info)
          lwork=nint(qwork(1))
          allocate(work(lwork))
          call dgeev(JOBVL,JOBVR,nk,hharm,nk,wr,wi,vl,nk,vr,nk,work,lwork,info)
          deallocate(work)
          do i=1,nk
             e(i)=sqrt(wr(i)**2+wi(i)**2)
          end do
          call IIInsertionSort(e(1:nk),ipos,nk)
          print*,e(1),e(nk)
          if (dabs(wi(ipos(km))).gt.1.d-15) then
             if (dabs(wi(ipos(km+1))+wi(ipos(km))).lt.1.d-15) then
                km=km+1
             end if
          end if
       allocate(wrh(km))
       kpos=km
       do i=1,km
           wrh(i)=wr(ipos(i))
       end do
       Hharm=H(1:nk,1:nk)+matmul(matmul(ihh,em),transpose(em))
       lwork=-1
       call dgees(JOBVR,'S',selectr,nk,hharm,nk,sdim,wr,wi,vr,nk,qwork,lwork,bwork,info)
       lwork=nint(qwork(1))
       allocate(work(lwork))
       call dgees(JOBVR,'S',selectr,nk,hharm,nk,sdim,wr,wi,vr,nk,work,lwork,bwork,info)
       deallocate(work)
       deallocate(wrh)
       !print*,km,sdim
       end do scoob
    end if
       
    
    
    
  end subroutine gmres3


  !Preconditioned conjugate gradient method with constant vector projected out. 
  !K matrix is not always positive definite although tends to be
  !not used but left here
  subroutine pconjgrad(n,Mm1,b,max_iterations,threshold,x)
    integer :: n
    real(8), intent(in) :: b(n),threshold,Mm1(n)
    integer, intent(in) :: max_iterations
    real(8), intent(inout) :: x(n)
    real(8), allocatable :: p(:),z(:),ap(:),r(:),xh(:),ph(:)
    real(8) :: rnew,rold,rnorm,alpha,elow
    integer :: k

    elow=1.d0
    allocate(p(n),z(n),ap(n),r(n),xh(n),ph(n))
    ap=1.d0/sqrt(dble(n))
     ph=ap
    !ap=ap/sqrt(dot_product(ap,ap))
    x=x-dot_product(x,ap)*ap
    !r=b-matmulW(A,ry,iy,dy,x,n)
    call aprodw(n,x,r)
    r=b-r
    r=r-dot_product(r,ap)*ap
    z=Mm1*r
    rold=dot_product(r,z)
    print*,'rold',rold
    p=z
    ph=ph/mm1
    ph=ph/sqrt(dot_product(ph,ph))
    screw:do k=1,max_iterations!*30
       !ap=matmulw(A,ry,iy,dy,p,n)
       p=p-dot_product(ph,p)*ph
       call aprodw(n,p,ap)
       !print*,dot_product(ap,ph)
       !pause
       alpha=rold/dot_product(p,ap)
       x=x+alpha*p
       r=r-alpha*ap
       rnorm=dot_product(r,r)
       !print*,'rnorm',k,rnorm
       !er(k)=sqrt(rnorm)
       !print*,sqrt(rnorm),rold
       if (sqrt(rnorm).lt.threshold.or.k==Max_Iterations) then
          !print*,'exited at',k
          exit screw
       end if
       if (sqrt(rold).lt.0.1*sqrt(rnorm)) then
          if (sqrt(rnorm).lt.1.d-20) then
             
          exit screw
          !print*,'stalled at',k
          end if
       end if
       if (sqrt(rnorm).lt.elow) then
          elow=sqrt(rnorm)
          xh=x
       end if
       if (sqrt(rnorm).gt.1.d4*elow) then
          x=xh
          exit screw
       end if
       z=Mm1*r
       rnew=dot_product(z,r)
       p=z+(rnew/rold)*p
       rold=rnew
       !ap=1.d0
       !ap=ap/sqrt(dot_product(r,r))
       !x=x-dot_product(x,ap)*ap
       !ap=1.d0/sqrt(dble(n))
       !ap=ap/sqrt(dot_product(ap,ap))
       !x=x-dot_product(x,ap)*ap
    end do screw
    !print*,x
    
  end subroutine pconjgrad

    !one arnoldi iteration  
  subroutine arnoldi(A,MM1,ry,iy,dy,Q,k,H)
    real(8), intent(in) :: A(:,:),ry(:),iy(:),dy(:),Mm1(:)
    integer, intent(in) :: k
    real(8), intent(inout) :: Q(:,:),H(:)
    real(8), allocatable :: qvec(:),qvecp(:)
    integer :: i,n3
    n3=size(Q,1)
    allocate(qvec(n3),qvecp(n3))
    qvecp=Mm1*Q(:,k)
    !qvec=matmulW(A,ry,iy,dy,Qvecp,n3)
    call aprodw(n3,qvecp,qvec)
    !qvec=matmul(A,Q(:,k))
    do i=1,k
       h(i)=dot_product(qvec,Q(:,i))
       qvec=qvec-h(i)*Q(:,i)
    end do
    h(k+1)=sqrt(dot_product(qvec,qvec))
    Q(:,k+1)=qvec/h(k+1)
    deallocate(qvec)
  end subroutine arnoldi

  !Givens rotation for GMRES
  subroutine apply_givens_rotation(h,cs,sn,k)
    real(8), intent(inout) :: h(:),cs(:),sn(:)
    integer, intent(in) :: k
    real(8) :: temp,cs_k,sn_k,t
    integer :: i

    do i=1,k-1
       temp=cs(i)*h(i)+sn(i)*h(i+1)
       h(i+1)=-sn(i)*h(i)+cs(i)*h(i+1)
       h(i)=temp
    end do

    if (abs(h(k)).lt.1.d-16) then
       cs_k=0.d0
       sn_k=1.d0
    else
       t=sqrt(h(k)**2+h(k+1)**2)
       cs_k=abs(h(k))/t
       sn_k=cs_k*h(k+1)/h(k)
    end if

    h(k)=cs_k*h(k)+sn_k*h(k+1)
    cs(k)=cs_k
    sn(k)=sn_k
    h(k+1)=0.d0
  end subroutine apply_givens_rotation

  
  !Rational function extrapolation for errors. Not currently used 
  subroutine rzextr(iest,imax,x,y,xest,yest,dy,nuse)
    integer :: iest,nuse,imax
    real(8) :: xest,yest,dy,y(imax),R(nuse,0:nuse)
    real(8) :: x(imax),Ra,Rb,Rc,xin(nuse)
    real(8) :: bot
    integer :: i,m,m1

    x(iest)=xest
    y(iest)=yest
    if (iest.le.3) then
       dy=yest
    elseif (abs(y(iest)-y(iest-1)).le.1.d-13) then
       dy=abs(y(iest)-y(iest-1))
    else
       m1=min(iest,nuse)
       do m=1,m1
          xin(m)=x(iest-m1+m)
          R(m,0)=y(iest-m1+m)
       end do
       do m=1,m1-1
         do i=m1-m,1,-1
            bot=xin(i)/xin(i+m)
            Ra=R(i,m-1)
            Rb=R(i+1,m-1)
            if (m.gt.2) then
               Rc=R(i+1,m-2)
            else
               Rc=0.d0
            end if
            R(i,m)=Rb+(Rb-Ra)/(bot*(1.d0-(Rb-Ra)/(Rb-Rc))-1.d0)
         end do
       end do
       yest=R(1,m1-1)
       dy=abs(yest-R(2,m1-2))
    end if
  end subroutine rzextr
