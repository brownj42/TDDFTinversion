module mapmod
  implicit none
contains
  subroutine calcmap(nd,ntotk,ind,sumn,plus,map)
    integer, intent(in) :: ntotk
    integer, intent(in) :: nd
    integer, intent(in) :: ind(nd,ntotk)
    integer, intent(out) :: sumn
    integer, intent(out) :: plus(nd)
    integer, allocatable, intent(out) :: map(:,:)
    integer :: c(nd),cu(0:nd),ci(1:nd)
    integer :: b,i
    sumn=nd*2
    allocate(map(sumn,0:ntotk))
    plus(1)=1
    do i=2,nd
       plus(i)=plus(i-1)+2
    end do
    map=0
    do i=1,ntotk
       ci=ind(:,i)
       cu(0)=1
       c(1:nd) = ci(1:nd)+plus
       do b=1,nd
          cu(b)=map(c(b),cu(b-1))
          if (cu(b)==0) then
             map(c(b),cu(b-1))=i
             cu(b)=i
          end if
       end do
    end do


  end subroutine calcmap

  integer function findpos(nd,sumn,ntotk,plus,ind,map)
    integer, intent(in) :: nd,sumn,ntotk
    integer, intent(in) :: plus(nd),ind(nd)
    integer, intent(in) :: map(sumn,0:ntotk)
    integer :: l,cin,c(nd)
    cin=1
    c=ind+plus
    do l=1,nd
       cin=map(c(l),cin)
    end do
    findpos=cin
  end function findpos
end module mapmod
