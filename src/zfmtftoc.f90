
! Copyright (C) 2016 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine zfmtftoc(nr,nri,zfmt,zfcmt)
use modmain
implicit none
! arguments
integer, intent(in) :: nr,nri
complex(8), intent(in) :: zfmt(*)
complex(8), intent(out) :: zfcmt(*)
! local variables
integer ir,i,j,n
i=1
j=1
n=lmmaxi*lradstp
do ir=1,nri,lradstp
  call zcopy(lmmaxi,zfmt(i),1,zfcmt(j),1)
  i=i+n
  j=j+lmmaxi
end do
i=i+(lradstp-1)*(lmmaxo-lmmaxi)
n=lmmaxo*lradstp
do ir=nri+lradstp,nr,lradstp
  call zcopy(lmmaxo,zfmt(i),1,zfcmt(j),1)
  i=i+n
  j=j+lmmaxo
end do
return
end subroutine

