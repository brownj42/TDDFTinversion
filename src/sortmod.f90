!Module that contains sorting algorithms
Module SORTMOD
Implicit None

Contains
  Subroutine rInsertionSort(A,n)
  Real(8), dimension(:) :: A
  Integer, intent(in) :: n
  Real(8) :: key
  Integer :: i,  j
  do j = 1,n
     key = A(j)
     i = j-1

     do 
        if ((i<=0)) exit
        if((A(i)<key))exit
           A(i+1) = A(i)
           i = i-1
     end do
     
     A(i+1) = key

  End do
  End Subroutine rInsertionSort
  Subroutine cInsertionSort(A,A2,n)
  real(8), dimension(:) :: A
  complex(8), dimension(:,:) :: A2
  Integer, intent(in) :: n
  real(8) :: key
  complex(8) :: key2(n)
  Integer :: i,j
  do j = 1,n
     key = A(j)
     key2=A2(:,j)
     i = j-1

     do 
        if ((i<=0)) exit
        if((A(i)<key))exit
           A(i+1) = A(i)
           A2(:,i+1)=A2(:,i)
           i = i-1
     end do
     
     A(i+1) = key
     A2(:,i+1)= key2

  End do
  End Subroutine cInsertionSort
  Subroutine InsertionSort(A,A2,n)
  Real(8), dimension(:) :: A
  real(8), dimension(:,:) :: A2
  Integer, intent(in) :: n
  Real(8) :: key
  Real(8) :: key2(n)
  Integer :: i,j
  do j = 1,n
     key = A(j)
     key2=A2(:,j)
     i = j-1

     do 
        if ((i<=0)) exit
        if((A(i)<key))exit
           A(i+1) = A(i)
           A2(:,i+1)=A2(:,i)
           i = i-1
     end do
     
     A(i+1) = key
     A2(:,i+1)= key2

  End do
  End Subroutine InsertionSort
  
  Subroutine IInsertionSort(A,ipos,n)
  Integer, intent(in) :: n
  Real(8) :: A(1:n)
  integer :: ipos(1:n)
  Real(8) :: key
  integer :: key3
  Integer :: i,j
  do j = 1,n
     key = A(j)
     key3=ipos(j)
     i = j-1

     do 
        if ((i<=0)) exit
        if((A(i)<key))exit
           A(i+1) = A(i)
           ipos(i+1)=ipos(i)
           i = i-1
     end do
     
     A(i+1) = key
     ipos(i+1)=key3

  End do
  END Subroutine IInsertionSort

  Subroutine IIInsertionSort(A,ipos,n)
  integer, intent(in) :: n
  real(8) :: A(1:n)
  integer :: ipos(1:n)
  real(8) :: key
  integer :: key3
  Integer :: i,j
  do j=1,n
  ipos(j)=j
  end do
  do j = n,1,-1
     key = A(j)
     key3=ipos(j)
     i = j+1

     do 
        if (i.gt.n) exit
        if((A(i).ge.key)) exit
           A(i-1) = A(i)
           ipos(i-1)=ipos(i)
           i = i+1
     end do
     
     A(i-1) = key
     ipos(i-1)=key3

  End do
  END Subroutine IIInsertionSort
  

  Subroutine Swap(A,B)
    Real(8), intent(inout) :: A, B
    Real(8) :: s
    
    s = A
    A = B
    B = s
  End Subroutine Swap


  Subroutine BubbleSort(B,m)
    Real(8), dimension(:) :: B
    Integer :: m,i,j
   
    do i=2,m
      do j=m,(i + 1),-1
          If (B(j)<B(j-1)) then
          call Swap(B(j),B(j-1))
          End IF
      End do
    End do
  End Subroutine BubbleSort


  real(8) function fact(i)
    integer :: i,j
    real(8) :: hold
    if (i.le.0) then
       fact=1.d0
    else
       hold=1.d0
       do j=1,i
          hold=hold*dble(j)
       end do
       fact=hold
    end if
  end function fact
End Module SORTMOD
