
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine olpaa(ias,ngp,apwalm,ld,o)
use modmain
implicit none
! arguments
integer, intent(in) :: ias,ngp
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw)
integer, intent(in) :: ld
complex(8), intent(inout) :: o(*)
! local variables
integer is,l,m,lm,io
! automatic arrays
complex(8) x(ngp)
is=idxis(ias)
lm=0
do l=0,lmaxapw
  do m=-l,l
    lm=lm+1
    do io=1,apword(l,is)
      x(1:ngp)=conjg(apwalm(1:ngp,io,lm))
      call zher('U',ngp,1.d0,x,1,o,ld)
    end do
  end do
end do
return
end subroutine

