module matmul_mod
  implicit none
integer, allocatable :: reord(:,:)
integer :: roswitch
contains
  function vtimesv(v1,v2,n3) result(v3)
    integer :: n3
    real(8),intent(in) :: v1(n3),v2(n3)
    real(8) :: v3(n3)
    integer :: i
    !$OMP parallel shared(v1,v2,v3) firstprivate(n3)
    !$OMP do schedule(static,4) 
    do i=1,n3
       v3(i)=v1(i)*v2(i)
    end do
    !$OMP end do nowait
    !$OMP end parallel
  end function vtimesv
  function vtimeszv(v1,v2,n3) result(v3)
    integer :: n3
    real(8), intent(in) :: v1(n3)
    complex(8),intent(in) :: v2(n3)
    complex(8) :: v3(n3)
    integer :: i
    !$OMP parallel shared(v1,v2,v3) firstprivate(n3)
    !$OMP do schedule(static,4) 
    do i=1,n3
       v3(i)=v1(i)*v2(i)
    end do
    !$OMP end do nowait
    !$OMP end parallel
  end function vtimeszv

  subroutine vplusab(v,a,b,n3)
    integer, intent(in) :: n3
    real(8), intent(in) :: b(n3),a
    real(8), intent(out) :: v(n3)
    integer :: i
    !$OMP parallel shared(v,b) firstprivate(a,n3)
    !$OMP do  
    do i=1,n3
       v(i)=v(i)+a*b(i)
    end do
    !$OMP end do nowait
    !$OMP end parallel
  end subroutine vplusab
   
  subroutine zvplusab(v,a,b,n3)
    integer, intent(in) :: n3
    complex(8), intent(in) :: b(n3),a
    complex(8), intent(out) :: v(n3)
    integer :: i
    !$OMP parallel shared(v,b) firstprivate(a,n3)
    !$OMP do
    do i=1,n3
       v(i)=v(i)+a*b(i)
    end do
    !$OMP end do nowait
    !$OMP end parallel
  end subroutine zvplusab

  function matmula(A,x,n3) result(cvec)
    real(8) :: A(:,:),x(n3)
    integer :: n3
    real(8) :: cvec(n3)
    real(8) :: cvecp(n3),xp(n3)
    integer :: n,n2,i1,i2,i3,i,pos,posp
    n=size(A,1)
    n2=n**2
    cvec=0.d0
    cvecp=0.d0
    if (roswitch==0) then
       roswitch=1
       do i1=1,n
          do i2=1,n
             do i3=1,n
                pos=(i1-1)*n2+(i2-1)*n+i3
                posp=(i2-1)*n2+(i3-1)*n+i1
                reord(pos,1)=posp
                posp=(i3-1)*n2+(i1-1)*n+i2
                reord(pos,2)=posp
                posp=(i1-1)*n2+(i2-1)*n+i3
                reord(pos,3)=posp
             end do
          end do
       end do
    end if
    !reorder to 2,3,1 order from 1,2,3 order
    do i=1,n3
       !cvecp(i)=cvec(reord(i,3))
       xp(reord(i,1))=x(i)
    end do
    !$OMP parallel shared(cvecp,xp) firstprivate(A,n,n2) default(private)
    !$OMP do collapse(3) schedule(static,4)
    !!$OMP do  
    do i1=1,n 
       do i2=1,n 
          do i3=1,n
             cvecp((i1-1)*n2+(i2-1)*n+i3)=dot_product(A(:,i3),xp((i1-1)*n2+(i2-1)*n+1:(i1-1)*n2+(i2-1)*n+n))
          end do
       end do
    end do
    !$OMP end do nowait
    !$OMP end parallel
    do i=1,n3
       cvec(i)=cvecp(reord(i,1))
    end do
    do i=1,n3
       cvecp(reord(i,2))=cvec(i)
       xp(reord(i,2))=x(i)
    end do
    !$OMP parallel shared(cvecp,xp) firstprivate(A,n,n2) default(private)
    !$OMP do collapse(3) schedule(static,4)
    !!$OMP do 
    do i1=1,n
       do i2=1,n
          do i3=1,n
             cvecp((i1-1)*n2+(i2-1)*n+i3)=cvecp((i1-1)*n2+(i2-1)*n+i3)+&
                  dot_product(A(:,i3),xp((i1-1)*n2+(i2-1)*n+1:(i1-1)*n2+(i2-1)*n+n))
          end do
       end do
    end do
    !$OMP end do nowait
    !$OMP end parallel
    do i=1,n3
       cvec(i)=cvecp(reord(i,2))
    end do
    !$OMP parallel shared(cvec,x) firstprivate(A,n,n2) default(private)
    !$OMP do collapse(3) schedule(static,4) 
    !!$OMP do 
    do i1=1,n
       do i2=1,n
          do i3=1,n
             cvec((i1-1)*n2+(i2-1)*n+i3)=cvec((i1-1)*n2+(i2-1)*n+i3)+&
                  dot_product(A(:,i3),x((i1-1)*n2+(i2-1)*n+1:(i1-1)*n2+(i2-1)*n+n))
          end do
       end do
    end do
    !$OMP end do nowait
    !$OMP end parallel
    
  end function matmula
 
  function matmulac(A,x,n3) result(cvec)
    real(8) :: A(:,:)
    integer :: n3
    complex(8) :: x(n3)
    complex(8) :: cvec(n3)
    complex(8) :: cvecp(n3),xp(n3)
    integer :: n,n2,i1,i2,i3,i,pos,posp
    n=size(A,1)
    n2=n**2
    cvec=0.d0
    cvecp=0.d0
    if (roswitch==0) then
       roswitch=1
       do i1=1,n
          do i2=1,n
             do i3=1,n
                pos=(i1-1)*n2+(i2-1)*n+i3
                posp=(i2-1)*n2+(i3-1)*n+i1
                reord(pos,1)=posp
                posp=(i3-1)*n2+(i1-1)*n+i2
                reord(pos,2)=posp
                posp=(i1-1)*n2+(i2-1)*n+i3
                reord(pos,3)=posp
             end do
          end do
       end do
    end if
    !reorder to 2,3,1 order from 1,2,3 order
    do i=1,n**3
       !cvecp(i)=cvec(reord(i,3))
       xp(reord(i,1))=x(i)
    end do
    !$OMP parallel shared(cvecp,xp) firstprivate(A,n,n2) default(private)
    !$OMP do
    do i1=1,n 
       do i2=1,n 
          do i3=1,n
             !do i=1,n
                !cvec((i1-1)*n**2+(i2-1)*n+i3)=&
                !     cvec((i1-1)*n**2+(i2-1)*n+i3)+&
                !     x((i-1)*n**2+(i2-1)*n+i3)*A(i1,i)
                cvecp((i1-1)*n2+(i2-1)*n+i3)=dot_product(A(:,i3),xp((i1-1)*n2+(i2-1)*n+1:(i1-1)*n2+(i2-1)*n+n))
             !end do
          end do
       end do
    end do
    !$OMP end do nowait
    !$OMP end parallel
    do i=1,n**3
       cvec(i)=cvecp(reord(i,1))
    end do
    do i=1,n**3
       cvecp(reord(i,2))=cvec(i)
       xp(reord(i,2))=x(i)
    end do
    !$OMP parallel shared(cvecp,xp) firstprivate(A,n,n2) default(private)
    !$OMP do
    do i1=1,n
       do i2=1,n
          do i3=1,n
             !do i=1,n
                !cvec((i1-1)*n**2+(i2-1)*n+i3)=cvec((i1-1)*n**2+(i2-1)*n+i3)+&
                !     x((i1-1)*n**2+(i-1)*n+i3)*A(i2,i)
             cvecp((i1-1)*n2+(i2-1)*n+i3)=cvecp((i1-1)*n2+(i2-1)*n+i3)+&
                  dot_product(A(:,i3),xp((i1-1)*n2+(i2-1)*n+1:(i1-1)*n2+(i2-1)*n+n))
             !end do
          end do
       end do
    end do
    !$OMP end do nowait
    !$OMP end parallel
    do i=1,n**3
       cvec(i)=cvecp(reord(i,2))
    end do
    !$OMP parallel shared(cvec,x) firstprivate(A,n,n2) default(private)
    !$OMP do
    do i1=1,n
       do i2=1,n
          do i3=1,n
             !do i=1,n
             cvec((i1-1)*n2+(i2-1)*n+i3)=cvec((i1-1)*n2+(i2-1)*n+i3)+&
                  dot_product(A(:,i3),x((i1-1)*n2+(i2-1)*n+1:(i1-1)*n2+(i2-1)*n+n))
                     !x((i1-1)*n**2+(i2-1)*n+i)*A(i3,i)
             !end do
          end do
       end do
    end do
    !$OMP end do nowait
    !$OMP end parallel
    
  end function matmulac


    !apply full dimensional T to dimensions individually for 2D problem
  function matmulac2(A,x,n3) result(cvec)
    real(8) :: A(:,:)
    integer :: n3
    complex(8) :: x(n3)
    complex(8) :: cvec(n3)
    complex(8) :: cvecp(n3),xp(n3)
    integer :: n,i1,i2,i,pos,posp
    n=size(A,1)
    cvec=0.d0
    cvecp=0.d0
    if (roswitch==0) then
       roswitch=1
       do i1=1,n
          do i2=1,n
             !do i3=1,n
                pos=(i1-1)*n+i2
                posp=(i2-1)*n+i1
                reord(pos,1)=posp
                posp=(i1-1)*n+i2
                reord(pos,2)=posp
                !posp=(i1-1)*n**2+(i2-1)*n+i3
                !reord(pos,3)=posp
             !end do
          end do
       end do
    end if
    !reorder to 2,1 order from 1,2 order
    do i=1,n**2
       !cvecp(i)=cvec(reord(i,3))
       xp(reord(i,1))=x(i)
    end do

    !apply matrix A to 2nd dimension of transformed vector
    !$OMP parallel shared(cvecp,xp) firstprivate(A,n) default(private)
    !$OMP do
    do i1=1,n 
       do i2=1,n 
          
          cvecp((i1-1)*n+i2)=dot_product(A(i2,:),xp((i1-1)*n+1:(i1-1)*n+n))
          
       end do
    end do
    
    !transform back to normal order
    !$OMP end do nowait
    !$OMP end parallel
    do i=1,n**2
       cvec(i)=cvecp(reord(i,1))
    end do

    !apply A to second dimension
    !$OMP parallel shared(cvec,x) firstprivate(A,n) default(private)
    !$OMP do
    do i1=1,n 
       do i2=1,n 
          
          cvec((i1-1)*n+i2)=cvec((i1-1)*n+i2)+dot_product(A(i2,:),x((i1-1)*n+1:(i1-1)*n+n))
          
       end do
    end do
    !$OMP end do nowait
    !$OMP end parallel
    
  end function matmulac2

    function diag_vec(v)
        implicit none 
        real(8) :: v(:)
        real(8), allocatable  :: diag_vec(:,:)
        integer :: i,s,f
        s=lbound(v,dim=1)
        f=ubound(v,dim=1)
        allocate(diag_vec(s:f,s:f))
        diag_vec=0.d0
        do i=s,f
        diag_vec(i,i)=v(i)
        end do
   end function diag_vec 
   function diag_mat(M)
     implicit none 
     real(8) :: M(:,:)
     real(8), allocatable  :: diag_mat(:)
     integer :: i,s,f
     s=lbound(M,dim=1)
     f=ubound(M,dim=1)
     allocate(diag_mat(lbound(M,dim=1):ubound(M,dim=1)))
     do i=s,f
        diag_mat(i)=M(i,i)
     end do
   end function diag_mat
  
end module matmul_mod
