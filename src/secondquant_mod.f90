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

  real(8) function getsign(pos2p,i,k)
    integer :: pos2p(:),i,k
    integer :: su
 
    su=0
    if (k.gt.1) then
    su=sum(pos2p(1:k-1))
    end if
    pos2p(k)=0
    if (i.gt.1) then
    su=su+sum(pos2p(1:i-1))
    end if
    pos2p(k)=1
    getsign=(-1.d0)**su
  end function getsign

  real(8) function getsign2(pos2p,i,j,k,l)
    integer :: pos2p(:),i,j,k,l
    integer :: su
    su=0
    if (l.gt.1) then
    su=sum(pos2p(1:l-1))
    end if
    pos2p(l)=0
    if (k.gt.1) then
    su=su+sum(pos2p(1:k-1))
    end if
    pos2p(k)=0
    if (j.gt.1) then
    su=su+sum(pos2p(1:j-1))
    end if
    pos2p(j)=1
    if (i.gt.1) then
    su=su+sum(pos2p(1:i-1))
    end if
    pos2p(k)=1
    pos2p(l)=1
    pos2p(j)=0
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
                   si=getsign(pos2p,i,k)
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
                                  si=getsign2(pos2p,i,j,l,k)
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
       vp(pp)=su
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
    integer :: pp,p,ip,jp,i,k,l,j,pij
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
          do k=1,nbasis
          !do k=max(i-1,1),min(i+1,nbasis)
             if (pos2pp(k)==0) then
                pos2pp(k)=1
                p=findpos(nbasis,sumn,numfci,plus,pos2pp,map)
                si=getsign(pos2pp,i,k)
                su=su+si*H1a(i,k)*v(p)
                pos2pp(k)=0
             end if
          end do
          pos2pp(i)=1
       end do
       !diagonal 2-body term as is required
       do ip=1,nab
          i=posrepp(ip)
          do jp=1,nab
             if (ip.ne.jp) then
                j=posrepp(jp)
                pij=(i-1)*nbasis+j
                k=i
                l=j
                si=getsign2(pos2pp,i,j,l,k)
                su=su+si*V2s(pij)*v(pp)
             end if
          end do
       end do
       vp(pp)=su
    end do
    !$OMP end do  
    !$OMP end parallel
    

  end function onematvecmap

  
end module secondquant_mod
