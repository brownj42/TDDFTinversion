module secondquant_wrap
  use secondquant_mod
  implicit none
contains
  subroutine onerdm(pos2,psi,numorb,nfock,rdm)
    integer, intent(in) :: numorb,nfock
    integer, intent(in) :: pos2(numorb,nfock)
    complex(8), intent(in) :: psi(nfock)
    complex(8), intent(out) :: rdm(numorb,numorb)
    integer :: i,o,op,fi(numorb),fk(numorb),k,nab
    real(8) :: si
    nab=sum(pos2(:,1))
    rdm=dcmplx(0.d0,0.d0)
    do o=1,nfock
       do i=1,numorb
          fi=pos2(:,o)
          
          if (fi(i)==1) then
             fi(i)=0
             do op=1,nfock
                
                do k=1,numorb
                   fk(:)=pos2(:,op)
                   if (fk(k)==1) then
                      fk(k)=0
                      
                      if (sum(fi*fk)==nab-1) then
                         si=(-1.d0)**sum(fi(1:i-1))*(-1.d0)**sum(fk(1:k-1))
                         rdm(i,k)=rdm(i,k)+psi(o)*dconjg(psi(op))*si
                      end if
                   end if
                end do
             end do
          end if
       end do
    end do
  end subroutine onerdm

  subroutine diagonerdm(sysparams,numorb,psi,rdm)
    use derivedtypes
    use nchoosekmod
    !function diagonerdm(pos2,psi,numorb,nfock) result(rdm)
    type(systemparameters), intent(in) :: sysparams
    integer, intent(in) :: numorb
    integer :: nfock
    complex(8), intent(in) :: psi(:)
    real(8), intent(out) :: rdm(numorb)
    integer :: i,o,op,fi(numorb)
    integer, allocatable :: pos2(:,:),posrep(:,:)
    real(8) :: si
    nfock=sysparams%ntot
    
    !nab=sysparams%npart
    call makefocklist(sysparams%np1,sysparams%npart,posrep,pos2,nfock)
    rdm=0.d0
    do o=1,nfock
       !print*,sysparams%pos2(:,o)
       fi=pos2(:,o)
       do i=1,numorb
          if (fi(i)==1) then
             !fi(i)=0
             do op=o,o
                !fk(:)=pos2(:,op)
                !do k=i,i
                   !if (fk(k)==1) then
                   !   fk(k)=0
                   !   if (sum(fi*fk)==nab-1) then
                         si=1.d0!(-1.d0)**sum(fi(1:i-1))*(-1.d0)**sum(fi(1:i-1))
                         rdm(i)=rdm(i)+dble(psi(o)*dconjg(psi(op)))*si
                   !   end if
                   !end if
                !end do
             end do
          end if
       end do
    end do
    !print*,rdm
    !pause
  end subroutine diagonerdm

  !subroutine diagdonerdm(pos2,posrep,psi,H1a,V2s,nab,numorb,nfock) result(rdm)
  subroutine diagdonerdm(sysparams,sharedvals,numorb,psi,rdm)
    use derivedtypes
    use nchoosekmod
    use mapmod
    type(systemparameters), intent(in) :: sysparams
    type(sharedvalues), intent(in) :: sharedvals
    integer, intent(in) :: numorb
    complex(8), intent(in) :: psi(:)
    real(8), intent(out) :: rdm(numorb)
    integer :: nfock
    !integer :: pos2(numorb,nfock),posrep(nab,nfock)
    !real(8) :: H1a(numorb,numorb),V2s(numorb**2)
    complex(8), allocatable :: dpsi(:)
    real(8), allocatable :: H1(:,:)
    integer :: i,o,op,fi(numorb)
    integer, allocatable :: pos2(:,:),posrep(:,:)
    integer, allocatable :: map(:,:),plus(:)
    integer :: sumn
    nfock=sysparams%ntot
    allocate(dpsi(nfock))
    allocate(H1(sysparams%np1,sysparams%np1))
    H1=sysparams%T
    do i=1,sysparams%np1
       H1(i,i)=H1(i,i)+sharedvals%v1(i)
    end do
    call makefocklist(sysparams%np1,sysparams%npart,posrep,pos2,nfock)
    allocate(plus(sysparams%np1))
    call calcmap(sysparams%ntot1,sysparams%ntot,pos2,sumn,plus,map)
    !dpsi=dcmplx(0.d0,1.d0)*onematvec(H1,sharedvals%vinteract,pos2,&
    !     posrep,sysparams%npart,numorb,nfock,psi)
    dpsi=dcmplx(0.d0,1.d0)*onematvecmap(H1,sharedvals%vinteract,pos2,&
                posrep,sysparams%npart,numorb,nfock,&
                    plus,sumn,map,psi)
    deallocate(map,posrep)
    rdm=0.d0
    do o=1,nfock
       fi=pos2(:,o)
       do i=1,numorb
          if (fi(i)==1) then
             !fi(i)=0
             do op=o,o
                !fk(:)=pos2(:,op)
                !do k=i,i
                   !if (fk(k)==1) then
                   !   fk(k)=0
                      !if (sum(fi*fk)==nab-1) then
                         !si=1.d0!(-1.d0)*sum(fi(1:i-1))*(-1.d0)*sum(fk(1:k-1))
                         rdm(i)=rdm(i)-dble(dpsi(o)*dconjg(psi(op))+psi(o)*dconjg(dpsi(op)))
                      !end if
                   !end if
                !end do
             end do
          end if
       end do
    end do
  end subroutine diagdonerdm

  !function diagddonerdm(pos2,posrep,psi,H1a,V2s,nab,numorb,nfock) result(rdm)
  subroutine diagddonerdm(sysparams,sharedvals,numorb,psi,rdm)
    use derivedtypes
    use nchoosekmod
    use mapmod
    type(systemparameters), intent(in) :: sysparams
    type(sharedvalues), intent(in) :: sharedvals
    integer, intent(in) :: numorb
    complex(8), intent(in) :: psi(:)
    real(8), intent(out) :: rdm(numorb)
    complex(8), allocatable :: dpsi(:),ddpsi(:)
    real(8), allocatable :: H1(:,:)
    integer :: nfock,i,o,op,fi(numorb)
    integer, allocatable :: pos2(:,:),posrep(:,:)
    integer, allocatable :: map(:,:),plus(:)
    integer :: sumn
    nfock=sysparams%ntot
    allocate(dpsi(nfock),ddpsi(nfock))
    allocate(H1(sysparams%np1,sysparams%np1))
    H1=sysparams%T
    do i=1,sysparams%np1
       H1(i,i)=H1(i,i)+sharedvals%v1(i)
    end do
    call makefocklist(sysparams%np1,sysparams%npart,posrep,pos2,nfock)
    allocate(plus(sysparams%np1))
    call calcmap(sysparams%ntot1,sysparams%ntot,pos2,sumn,plus,map)
    !dpsi=onematvec(H1,sharedvals%vinteract,pos2,posrep,&
    !     sysparams%npart,numorb,nfock,psi)
    !ddpsi=onematvec(H1,sharedvals%vinteract,pos2,posrep,&
    !     sysparams%npart,numorb,nfock,dpsi)
    dpsi=onematvecmap(H1,sharedvals%vinteract,pos2,&
                posrep,sysparams%npart,numorb,nfock,&
                    plus,sumn,map,psi)
    ddpsi=onematvecmap(H1,sharedvals%vinteract,pos2,&
                posrep,sysparams%npart,numorb,nfock,&
                    plus,sumn,map,dpsi)
    
    rdm=0.d0
    do o=1,nfock
       fi=pos2(:,o)
       do i=1,numorb
          if (fi(i)==1) then
             !fi(i)=0
             do op=o,o
                !fk(:)=pos2(:,op)
                !do k=i,i
                   !if (fk(k)==1) then
                   !   fk(k)=0
                      !if (sum(fi*fk)==nab-1) then
                         !si=1.d0!(-1.d0)*sum(fi(1:i-1))*(-1.d0)*sum(fk(1:k-1))
                rdm(i)=rdm(i)-dble(ddpsi(o)*dconjg(psi(op))+psi(o)*dconjg(ddpsi(op))-&
                     dconjg(dpsi(o))*dpsi(op)-dpsi(o)*dconjg(dpsi(op)))
                      !end if
                   !end if
                !end do
             end do
          end if
       end do
    end do
  end subroutine diagddonerdm

end module secondquant_wrap
