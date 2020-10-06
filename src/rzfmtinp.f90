
! Copyright (C) 2014 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

complex(8) function rzfmtinp(nr,nri,r3,rfmt,zfmt)
use modmain
implicit none
! arguments
integer, intent(in) :: nr,nri
real(8), intent(in) :: r3(nr)
real(8), intent(in) :: rfmt(*)
complex(8), intent(in) :: zfmt(*)
! local variables
integer n,ir,i
real(8) t1
complex(8) z1
! automatic arrays
real(8) fr1(nr),fr2(nr)
! external functions
real(8) splint
external splint
n=lmmaxi-1
i=1
do ir=1,nri
  z1=dot_product(rfmt(i:i+n),zfmt(i:i+n))
  fr1(ir)=dble(z1); fr2(ir)=aimag(z1)
  i=i+lmmaxi
end do
t1=dble(lmmaxi)/dble(lmmaxo)
n=lmmaxo-1
do ir=nri+1,nr
  z1=t1*dot_product(rfmt(i:i+n),zfmt(i:i+n))
  fr1(ir)=dble(z1); fr2(ir)=aimag(z1)
  i=i+lmmaxo
end do
! integrate over r
t1=fourpi/dble(3*lmmaxi)
rzfmtinp=t1*cmplx(splint(nr,r3,fr1),splint(nr,r3,fr2),8)
return
end function

