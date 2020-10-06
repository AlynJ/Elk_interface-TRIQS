
! Copyright (C) 2015 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

real(8) function fintgt(m,n,x,f)
implicit none
! arguments
integer, intent(in) :: m,n
real(8), intent(in) :: x(n),f(n)
! local variables
integer i
real(8) x0,x1,x2
! external functions
real(8) splint
external splint
if (n.le.0) then
  write(*,*)
  write(*,'("Error(fintgt): invalid number of points : ",I8)') n
  write(*,*)
  stop
end if
select case(m)
! low accuracy trapezoidal integration
case(-3)
  fintgt=0.d0
  do i=1,n-1
    fintgt=fintgt+0.5d0*(x(i+1)-x(i))*(f(i+1)+f(i))
  end do
  return
case(-2)
! medium accuracy Simpson integration
  fintgt=0.d0
  do i=1,n-2
    x0=x(i)
    x1=x(i+1)
    x2=x(i+2)
    fintgt=fintgt &
     +(x0-x1)*(f(i+2)*(x0-x1)**2+f(i+1)*(x2-x0)*(x0+2.d0*x1-3.d0*x2) &
     +f(i)*(x2-x1)*(2.d0*x0+x1-3.d0*x2))/(6.d0*(x0-x2)*(x1-x2))
  end do
  x0=x(n)
  x1=x(n-1)
  x2=x(n-2)
  fintgt=fintgt+(x1-x0)*(f(n-2)*(x1-x0)**2+f(n)*(x1-x2)*(3.d0*x2-x1-2.d0*x0) &
   +f(n-1)*(x0-x2)*(3.d0*x2-2.d0*x1-x0))/(6.d0*(x2-x1)*(x2-x0))
  return
case(-1)
! high accuracy integration from spline interpolation
  fintgt=splint(n,x,f)
  return
case default
  write(*,*)
  write(*,'("Error(fintgt): undefined m : ",I8)') m
  write(*,*)
  stop
end select
end function

