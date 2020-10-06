!Written by A. D. N. James

subroutine wanolp(ik,ias,nst,l,lmmax,u,wfmt,wantemp)
!Calculates the overlap matrix elements between the wavfunctions and (apw)
!radial function
use modmain
implicit none
!arguments
integer, intent(in) :: ik
integer, intent(in) :: ias
integer, intent(in) :: nst
integer, intent(in) :: l
integer, intent(in) :: lmmax
real(8), intent(in) :: u(nrmtmax)
complex(8), intent(in) :: wfmt(npcmtmax,natmtot,nspinor,nst)
complex(8), intent(out) :: wantemp(lmmax,nst,nspinor)
! local variables
integer is,nrc,nrci,nrco,npc,npci,ispn
integer idx,lm,m,ist,nr
real(8) t1,t2
real(8), allocatable :: r(:)
! automatic arrays
real(8) fr(npcmtmax), imfr(npcmtmax), wave(npcmtmax)
! external functions
real(8) fintgt
external fintgt
!correlated orb ias index
is=idxis(ias)
if (lradstp.eq.1) then
  nrc=nrmt(is)
  nrci=nrmti(is)
  npc=npmt(is)
  npci=npmti(is)
!radial mesh for integral
  allocate(r(nrc))
  r(:)=rsp(1:nrc,is)
else 
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  npc=npcmt(is)
  npci=npcmti(is)
  allocate(r(nrc))
  r(:)=rcmt(1:nrc,is)
end if
nr=nrmt(is)
nrco=nrc-nrci
!Calculate the temporary Wannier matrix elements
do ispn=1,nspinor
  do ist=1,nst
    idx=1
    do m=-l,l
      lm=idxlm(l,m)
      fr(:)=0.d0
      imfr(:)=0.d0
! calcualte overlap for the real part
      wave(:)=dble(wfmt(:,ias,ispn,ist))
      call ovlp(l,lm,is,npci,nr,nrci,nrc,wave,u,fr)
! calcualte overlap for the imaginary part
      wave(:)=aimag(wfmt(:,ias,ispn,ist))
      call ovlp(l,lm,is,npci,nr,nrci,nrc,wave,u,imfr)
! integrate over radius to produce temporary Wannier projector matrix elements
      t1=0.d0
      t2=0.d0
      t1=fintgt(-1,nrc,r(1:nrc),fr)
      t2=fintgt(-1,nrc,r(1:nrc),imfr)
! put real and imaginary elements into array   
      wantemp(idx,ist,ispn)=cmplx(t1,t2,8)
      idx=idx+1
    end do
  end do
end do
if(allocated(r)) deallocate(r)
return

contains 

subroutine ovlp(l,lm,is,npci,nr,nrci,nrc,wf,u,fr)
use modmain
implicit none
!arguments
integer, intent(in) :: l
integer, intent(in) :: lm
integer, intent(in) :: is
integer, intent(in) :: npci
integer, intent(in) :: nr
integer, intent(in) :: nrci
integer, intent(in) :: nrc
real(8), intent(in) :: wf(*)
real(8), intent(in) :: u(nr)
real(8), intent(out) :: fr(nrc)
!local
integer lmmt,i,ir,iir,iro,ldi,ldo

ldi=lmmaxi
ldo=lmmaxo
i=npci+lm
iro=nrmti(is)+lradstp
!Generate temporary projectors from wfmt
if (l.le.lmaxi) then
  lmmt=lm
  iir=1
  do ir=1,nrci          
    if (abs(dble(wf(lmmt))).gt.1.d-14) then
      fr(ir)=rlmt(iir,2,is)*u(iir)*wf(lmmt)
     endif
     lmmt=lmmt+ldi
     iir=iir+lradstp
  end do
end if
iir=iro
lmmt=i
do ir=nrci+1,nrc          
  if (abs(dble(wf(lmmt))).gt.1.d-14) then
    fr(ir)=rlmt(iir,2,is)*u(iir)*wf(lmmt)
  endif
  lmmt=lmmt+ldo
  iir=iir+lradstp
end do
return
end subroutine

end subroutine
!EOC

