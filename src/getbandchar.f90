
! Copyright (C) 2002-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!rewritten by A. D. N. James

subroutine getbandchar(tspndg,tlmdg,lmmax,bc)
use modmain
use modmpi
use modomp
implicit none
! arguments
logical, intent(in) :: tspndg,tlmdg
integer, intent(in) :: lmmax
real(4), intent(out) :: bc(lmmax,nspinor,natmtot,nstsv,nkpt)


! local variables
logical tsqaz
integer nsk(3),ik,jk,ist,iw,ld
integer nsd,ispn,jspn,is,ia,ias
integer l0,l1,l,lm,nthd
real(8) dw,th,sps(2),vl(3),vc(3)
real(8) v1(3),v2(3),v3(3),t1
complex(8) su2(2,2),b(2,2),c(2,2)
character(256) fname
! low precision for band/spin character array saves memory
real(8), allocatable :: elm(:,:)
complex(8), allocatable :: ulm(:,:,:),a(:,:)
complex(8), allocatable :: dmat(:,:,:,:,:)
complex(8), allocatable :: apwalm(:,:,:,:,:)
complex(8), allocatable :: evecfv(:,:,:),evecsv(:,:)
ld=lmmax*nspinor
if (dosssum) then
  nsd=1
else
  nsd=nspinor
end if
if (dosmsum) then
  l0=0; l1=lmaxdos
else
  l0=1; l1=lmmax
end if
! generate unitary matrices which convert the (l,m) basis into the irreducible
! representation basis of the symmetry group at each atomic site
if (lmirep) then
  allocate(elm(lmmax,natmtot))
  allocate(ulm(lmmax,lmmax,natmtot))
  call genlmirep(lmaxdos,lmmax,elm,ulm)
end if
! compute the SU(2) operator used for rotating the density matrix to the
! desired spin-quantisation axis
v1(:)=sqados(:)
t1=sqrt(v1(1)**2+v1(2)**2+v1(3)**2)
if (t1.le.epslat) then
  write(*,*)
  write(*,'("Error(dos): spin-quantisation axis (sqados) has zero length")')
  write(*,*)
  stop
end if
v1(:)=v1(:)/t1
if (v1(3).ge.1.d0-epslat) then
  tsqaz=.true.
else
  tsqaz=.false.
  v2(1:2)=0.d0
  v2(3)=1.d0
  call r3cross(v1,v2,v3)
! note that the spin-quantisation axis is rotated, so the density matrix should
! be rotated in the opposite direction
  th=-acos(v1(3))
  call axangsu2(v3,th,su2)
end if

! begin parallel loop over k-points
call holdthd(nkpt,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(apwalm,evecfv,evecsv,dmat,a) &
!$OMP PRIVATE(jk,ispn,jspn,vl,vc) &
!$OMP PRIVATE(is,ia,ias,ist,lm,b,c,t1) &
!$OMP NUM_THREADS(nthd)
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
allocate(evecfv(nmatmax,nstfv,nspnfv),evecsv(nstsv,nstsv))
allocate(dmat(lmmax,nspinor,lmmax,nspinor,nstsv),a(lmmax,lmmax))
!$OMP DO
do ik=1,nkpt
! loop over first-variational spins
  do ispn=1,nspnfv
    vl(:)=vkl(:,ik)
    vc(:)=vkc(:,ik)
! spin-spiral case
    if (spinsprl) then
      if (ispn.eq.1) then
        vl(:)=vl(:)+0.5d0*vqlss(:)
        vc(:)=vc(:)+0.5d0*vqcss(:)
      else
        vl(:)=vl(:)-0.5d0*vqlss(:)
        vc(:)=vc(:)-0.5d0*vqcss(:)
      end if
    end if
! find the matching coefficients
    call match(ngk(ispn,ik),vgkc(:,:,ispn,ik),gkc(:,ispn,ik), &
     sfacgk(:,:,ispn,ik),apwalm(:,:,:,:,ispn))
  end do
! get the eigenvectors from file for non-reduced k-point
  call getevecfv('.OUT',0,vkl(:,ik),vgkl(:,:,:,ik),evecfv)
  call getevecsv('.OUT',0,vkl(:,ik),evecsv)
  do is=1,nspecies
    do ia=1,natoms(is)
      ias=idxas(ia,is)
! generate the density matrix
      call gendmatk(.false.,.false.,0,lmaxdos,ias,ngk(:,ik),apwalm,evecfv, &
       evecsv,lmmax,dmat)
! convert (l,m) part to an irreducible representation if required
      if (lmirep) then
        do ist=1,nstsv
          do ispn=1,nspinor
            do jspn=1,nspinor
              call zgemm('N','N',lmmax,lmmax,lmmax,zone,ulm(:,:,ias),lmmax, &
               dmat(:,ispn,1,jspn,ist),ld,zzero,a,lmmax)
              call zgemm('N','C',lmmax,lmmax,lmmax,zone,a,lmmax,ulm(:,:,ias), &
               lmmax,zzero,dmat(:,ispn,1,jspn,ist),ld)
            end do
          end do
        end do
      end if
! spin rotate the density matrices to desired spin-quantisation axis
      if (spinpol.and.(.not.tsqaz)) then
        do ist=1,nstsv
          do lm=1,lmmax
            b(:,:)=dmat(lm,:,lm,:,ist)
            call z2mm(su2,b,c)
            call z2mmct(c,su2,b)
            dmat(lm,:,lm,:,ist)=b(:,:)
          end do
        end do
      end if
! determine the band characters from the density matrix
    do ist=1,nstsv
      do ispn=1,nspinor
          do lm=1,lmmax
            t1=dble(dmat(lm,ispn,lm,ispn,ist))
            bc(lm,ispn,ias,ist,ik)=real(t1)
          end do
        end do
      end do
    end do
  end do
end do
!$OMP END DO
deallocate(apwalm,evecfv,evecsv,dmat,a)
!$OMP END PARALLEL
call freethd(nthd)
if (lmirep) deallocate(elm,ulm)
return
end subroutine
!EOC

