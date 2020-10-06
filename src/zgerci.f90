
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

pure subroutine zgerci(m,n,x,y,ld,a)
implicit none
! arguments
integer, intent(in) :: m,n
complex(8), intent(in) :: x(m),y(n)
integer, intent(in) :: ld
complex(8), intent(inout) :: a(ld,*)
! local variables
integer i,j
! numbers less than eps are considered to be zero
real(8), parameter :: eps=1.d-10
real(8) a1,b1
complex(8) z1
do j=1,n
  z1=y(j)
  if (abs(dble(z1)).gt.eps) then
    if (abs(aimag(z1)).gt.eps) then
! complex prefactor
      do i=1,m
        a(i,j)=a(i,j)+z1*x(i)
      end do
    else
! real prefactor
      a1=dble(z1)
      do i=1,m
        a(i,j)=a(i,j)+a1*x(i)
      end do
    end if
  else if (abs(aimag(z1)).gt.eps) then
! imaginary prefactor
    b1=aimag(z1)
    do i=1,m
      a(i,j)=a(i,j)+b1*cmplx(-aimag(x(i)),dble(x(i)),8)
    end do
  end if
end do
return
end subroutine

