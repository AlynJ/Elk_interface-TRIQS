
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

real(8) function rfint(rfmt,rfir)
use modmain
implicit none
! arguments
real(8), intent(in) :: rfmt(npmtmax,natmtot),rfir(ngtot)
! local variables
integer is,ias,i
integer nr,nri,ir
real(8) t1
! automatic arrays
real(8) fr(nrmtmax)
! external functions
real(8) splint
external splint
! interstitial contribution
rfint=dot_product(rfir(:),cfunir(:))
rfint=rfint*omega/dble(ngtot)
! muffin-tin contribution
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmti(is)
  i=1
  do ir=1,nri
    fr(ir)=rfmt(i,ias)*rlmt(ir,2,is)
    i=i+lmmaxi
  end do
  do ir=nri+1,nr
    fr(ir)=rfmt(i,ias)*rlmt(ir,2,is)
    i=i+lmmaxo
  end do
  t1=splint(nr,rlmt(:,1,is),fr)
  rfint=rfint+fourpi*y00*t1
end do
return
end function

