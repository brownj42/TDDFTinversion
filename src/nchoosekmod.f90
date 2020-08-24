module nchoosekmod
  implicit none
contains
  subroutine makefocklist(M,N,plist,flist,tot)
    integer, intent(in) :: M,N
    integer, allocatable :: plist(:,:),flist(:,:)
    integer, intent(out) :: tot
    integer:: the_list(M),str_builder(N),jin,prev,l
    do l=1,M
       the_list(l)=l
    end do
    jin=0
    prev=0 
    tot=0
    call countnchoosek(M,N,jin,prev,tot)
    allocate(plist(N,tot))
    jin=0
    prev=0 
    tot=0
    call buildnchoosek(the_list,N,str_builder,plist,jin,prev,tot)
    allocate(flist(M,tot))
    flist=0
    do l=1,tot
       do jin=1,N
          flist(plist(jin,l),l)=1
       end do
    end do
  end subroutine makefocklist

  recursive subroutine countnchoosek(size_the_list,kin,j,prev,tot)
    integer,intent(in) :: size_the_list !the number of elements from which to choose
    integer,intent(in) :: kin !the total number of elements to choose
    integer :: tot !the total number of nchoosek
    integer :: j !the current element of str_builder we are adding to
    integer :: prev !the previous element position added
    integer :: i !counter
    if (j==kin) then
       tot=tot+1
       return
    else
       do i=prev+1,size_the_list
          j=j+1
          call countnchoosek(size_the_list,kin,j,i,tot)
          j=j-1
       end do
    end if
  end subroutine countnchoosek
  
  recursive subroutine buildnchoosek(the_list,kin,str_builder,list,j,prev,tot)
    integer,intent(in) :: the_list(:) !the list of elements from which to choose
    integer,intent(in) :: kin !the total number of elements to choose
    integer :: tot !the total number of nchoosek
    integer :: list(:,:),str_builder(kin) !the build list
    integer :: j !the current element of str_builder we are adding to
    integer :: prev !the previous element position added
    integer :: i !counter
    if (j==kin) then
       tot=tot+1
       list(:,tot)=str_builder
       return
    else
       do i=prev+1,size(the_list)
          j=j+1
          str_builder(j)=the_list(i)
          call buildnchoosek(the_list,kin,str_builder,list,j,i,tot)
          j=j-1
       end do
    end if
  end subroutine buildnchoosek
end module nchoosekmod
