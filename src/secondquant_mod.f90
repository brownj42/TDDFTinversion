module secondquant_mod
  implicit none
contains
  
  integer function fcheck(pos2,pp,p,i,j,k,l,nbasis)
    integer :: pos2(:,:),pp,p,i,j,k,l,nbasis
    integer :: b
    fcheck=0
    do b=1,nbasis
       if (b.ne.i.and.b.ne.j.and.b.ne.k.and.b.ne.l) then
          fcheck=fcheck+abs(pos2(b,pp)-pos2(b,p))
       end if
    end do
  end function fcheck

  integer function fcheck1(pos2,pp,p,i,k,nbasis)
    integer :: pos2(:,:),pp,p,i,k,nbasis
    integer :: b
    fcheck1=0
    do b=1,nbasis
       if (b.ne.i.and.b.ne.k) then
          fcheck1=fcheck1+abs(pos2(b,pp)-pos2(b,p))
       end if
    end do
  end function fcheck1

  real(8) function getsign(pos2p,i,k, length)
    integer :: length
    integer :: pos2p(length),i,k
    integer :: su
    integer :: vals(length)
    vals = pos2p
    su=0
    if (k.gt.1) then
    su=sum(vals(1:k-1))
    end if
    vals(k)=0
    if (i.gt.1) then
    su=su+sum(vals(1:i-1))
    end if
    vals(k)=1
    getsign=(-1.d0)**su
  end function getsign

  real(8) function getsign2(pos2p,i,j,k,l,length)
    integer :: length
    integer :: pos2p(length),i,j,k,l
    integer :: su
    integer :: vals(length)
    vals = pos2p
    su=0
    if (l.gt.1) then
    su=sum(vals(1:l-1))
    end if
    vals(l)=0
    if (k.gt.1) then
    su=su+sum(vals(1:k-1))
    end if
    vals(k)=0
    if (j.gt.1) then
    su=su+sum(vals(1:j-1))
    end if
    vals(j)=1
    if (i.gt.1) then
    su=su+sum(vals(1:i-1))
    end if
    vals(k)=1
    vals(l)=1
    vals(j)=0
    getsign2=(-1.d0)**su
  end function getsign2

  function onematvec(H1a,V2s,pos2,posrep,nab,nbasis,numfci,v) result(vp)
    integer :: nbasis,numfci,nab
    real(8) :: H1a(nbasis,nbasis),V2s(nbasis*nbasis)
    integer :: pos2(nbasis,numfci),pos2p(nbasis)
    integer :: posrep(nab,numfci)
    complex(8) :: v(numfci),vp(numfci),su,vv
    integer :: pp,p,ip,i,k,kp,lp,l,j,pij,pkl
    real(8) :: si
             
    do pp=1,numfci
       su=dcmplx(0.d0,0.d0)
       do p=1,numfci
          vv=v(p)
          pos2p=pos2(:,p)
          do ip=1,nab
             i=posrep(ip,pp)
             do kp=1,nab
                k=posrep(kp,p)
                if (fcheck1(pos2,pp,p,i,k,nbasis)==0) then
                   si=getsign(pos2p,i,k,nbasis)
                   su=su+si*H1a(i,k)*vv
                end if
             end do
          end do
          do ip=1,nab
             i=posrep(ip,pp)
             do j=1,nbasis
                if (pos2(j,pp)==1) then
                   pij=(i-1)*nbasis+j
                   do kp=1,nab
                      k=posrep(kp,p)
                      do lp=1,nab
                         l=posrep(lp,p)
                         if (k.ne.l) then
                            pkl=(k-1)*nbasis+l
                            if (i==k.and.j==l) then
                               if (fcheck(pos2,pp,p,i,j,k,l,nbasis)==0) then
                                  si=getsign2(pos2p,i,j,l,k,nbasis)
                                  su=su+si*V2s(pij)*vv
                               end if
                            end if
                         end if
                      end do
                   end do
                end if
             end do
          end do
       end do
       vp(pp)=dble(su)
    end do
    !print*,vp
    !pause 

  end function onematvec

  function onematvecmap(H1a,V2s,pos2,posrep,nab,nbasis,numfci,plus,sumn,map,v) result(vp)
    use mapmod
    integer :: nbasis,numfci,nab
    real(8) :: H1a(nbasis,nbasis),V2s(nbasis*nbasis)
    integer :: pos2(nbasis,numfci),pos2pp(nbasis)
    integer :: posrep(nab,numfci),posrepp(nab)
    integer :: plus(nbasis),sumn,map(sumn,0:numfci)
    complex(8) :: v(numfci),vp(numfci),su
    integer :: pp,p,ip,jp,i,k,l,j,pij,ks,kf,ls,lf
    real(8) :: si
    

    
    !$OMP parallel shared(vp,pos2,posrep) firstprivate(H1a,V2s,numfci,v,sumn,&
    !$OMP plus,map,nab,nbasis) default(private)
    !$OMP do schedule(static,4)
    
    do pp=1,numfci
       su=dcmplx(0.d0,0.d0)
       pos2pp=pos2(:,pp)
       posrepp=posrep(:,pp)
       do ip=1,nab
          i=posrepp(ip)
          pos2pp(i)=0
          if (i.le.nbasis/2) then
            ks = 1
            kf = nbasis/2
          else
            ks = nbasis/2+1
            kf = nbasis
          end if
          do k=ks,kf
          !do k=max(i-1,1),min(i+1,nbasis)
             if (pos2pp(k)==0) then
                pos2pp(k)=1
                p=findpos(nbasis,sumn,numfci,plus,pos2pp,map)
                si=getsign(pos2pp,i,k,nbasis)
                su=su+si*H1a(i,k)*v(p)
                pos2pp(k)=0
             end if
          end do
          pos2pp(i)=1
       end do
       !diagonal 2-body term as is required
       pos2pp = pos2(:, pp)
       posrepp=posrep(:,pp)
       do ip=1,nab
          i=posrepp(ip)
          do jp=1,nab
             j = posrepp(jp)
             if (ip.ne.jp) then ! .and. pos2(j, pp)==1) then
                !j=posrepp(jp)
                !if (i .le. nbasis/2) then
                !  ks = 1
                !  kf = nbasis/2
                !else
                !  ks = nbasis/2+1
                !  kf = nbasis
                !end if
                !if (j .le. nbasis/2) then
                !  ls = 1
                !  lf = nbasis/2
                !else
                !  ls = nbasis/2+1
                !  lf = nbasis
                !end if
                !do k = ks,kf
                !  do l = ls,lf
                  !k=i
                  !l=j
                  !   if (k==i .and. l==j) then
                        pij=(i-1)*nbasis+j
                        si=getsign2(pos2pp,i,j,j,i,nbasis)
                        su=su+si*V2s(pij)*v(pp)/2
                  !   end if
                !  end do
                !end do
                !k=i
                !l=j
                ! up, up
                !if (i .le. nbasis/2 .and. j .le. nbasis/2) then
                !  ! k = i, l = j
                !  pij=(i-1)*nbasis+l
                !  si=getsign2(pos2pp,i,j,l,k,nbasis)
                !  su=su+si*V2s(pij)*v(pp)/2
                !  ! k = j, l = i
                !  pij=(i-1)*nbasis+k
                !  si=getsign2(pos2pp,i,j,k,l,nbasis)
                !  su=su+si*V2s(pij)*v(pp)/2
                ! down, down
                !elseif (i .gt. nbasis/2 .and. j .gt. nbasis/2) then
                  ! k = i, l = j
                !  pij=(i-1)*nbasis+l
                !  si=getsign2(pos2pp,i,j,l,k,nbasis)
                !  su=su+si*V2s(pij)*v(pp)/2
                !  ! k = j, l = i
                !  pij=(i-1)*nbasis+k
                !  si=getsign2(pos2pp,i,j,k,l,nbasis)
                !  su=su+si*V2s(pij)*v(pp)/2
                ! up, down down up
                !elseif (i .le. nbasis/2 .and. j .gt. nbasis/2) then
                !  ! k = i, l = j
                !  pij=(i-1)*nbasis+l
                !  si=getsign2(pos2pp,i,j,l,k,nbasis)
                !  su=su+si*V2s(pij)*v(pp)/2
                !elseif (i .gt. nbasis/2 .and. j .le. nbasis/2) then
                !  ! k = i, l = j
                !  pij=(i-1)*nbasis+l
                !  si=getsign2(pos2pp,i,j,l,k,nbasis)
                !  su=su+si*V2s(pij)*v(pp)/2
                !end if

             end if
          end do
       end do
       vp(pp)=su
    end do
    !$OMP end do  
    !$OMP end parallel
    

  end function onematvecmap

  function calcf1_sq(nab,np1,ntot,dt,T,v2s,vks,dvks,phi,pos2,posrep,plus,sumn,map) result(f)
   integer, intent(in) :: nab,np1, ntot
   real(8) :: dt,T(np1,np1), v2s(np1*np1), vks(ntot),dvks(ntot)
   complex(8) :: phi(ntot)
   complex(8) :: f(ntot)
   integer :: pos2(np1,ntot)
   integer :: posrep(nab,ntot),posrepp(nab)
   integer :: sumn
   integer :: plus(np1),map(sumn,0:ntot)
   complex(8) :: ay(ntot),nwvd(ntot,0:1),Tnwvd0(ntot)

   nwvd(:,0)=vks*phi !N

   Tnwvd0=onematvecmap(T,V2s,pos2,posrep,nab,np1,ntot,plus,sumn,map,nwvd(:,0))
   ay=-dcmplx(0.d0,1.d0)*(onematvecmap(T,V2s,pos2,posrep,nab,np1,ntot,plus,sumn,map,phi)+vks*phi)

   nwvd(:,1)=dvks*phi 
   nwvd(:,1)=vks*ay+nwvd(:,1) !Nprime

   f=dt**2*dcmplx(1.d0,0.d0)/12.d0*(dcmplx(-1.d0,0.d0)*&
        (Tnwvd0)+dcmplx(0.d0,1.d0)*&
        nwvd(:,1))
 end function calcf1_sq
  
end module secondquant_mod
