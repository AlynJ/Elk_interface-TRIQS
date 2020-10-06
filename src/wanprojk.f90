!Written by A. D. N. James

subroutine wanprojk(ik,idx,ld,nproj,mst,nst,subulm,sublm,symlm,wanprj)
!this routine generates the wannier projectors 
use modmain

implicit none
!arguments
integer, intent(in) :: ik
integer, intent(in) :: idx(nstsv,nkpt)
integer, intent(in) :: ld
integer, intent(in) :: nproj
integer, intent(in) :: mst
integer, intent(in) :: nst 
complex(8), intent(in) :: subulm(ld,ld,norb,nproj)
integer, intent(in) :: sublm(ld,iorb,norb,nproj)
complex(8), intent(in) :: symlm(ld,ld,norb,nproj)
complex(8), intent(inout) :: wanprj(ld,mst,nspinor,norb,nproj) 
!local
integer ist,ispn,lm,is,nr,ikp,iorb,ia,ias,l,lmmax
integer ilo !temp
real(8) t1
complex(8), allocatable :: evecfv(:,:,:),evecsv(:,:)
real(8), allocatable :: evalfv(:,:)
integer, allocatable :: igpig(:,:) 
complex(8), allocatable :: wfmt(:,:,:,:),wfir(:,:,:)
complex(8), allocatable :: wantemp(:,:,:,:,:)
complex(8), allocatable :: z(:),z1(:,:),z2(:,:)
real(8), allocatable :: u(:)
integer ngp(nspnfv)

!k index for reading in eigenvectors in genwfsvpwan
ikp=0
if((task.eq.807).or.(task.eq.820)) ikp=ik
! this generates the second variational wavefunctions in the spherical harmoincs basis and the 
! intertitial wavefunctions are in real space.
allocate(igpig(ngkmax,nspnfv))
allocate(wfmt(npcmtmax,natmtot,nspinor,nst),wfir(ngtc,nspinor,nst))
call genwfsvpwan(.true.,.false.,ikp,nst,idx(1:nst,ik),ngdc,igfc,vkl(:,ik),ngp,igpig,wfmt,ngtc, &
       wfir)
!intersitial wavefunction not used in projection
deallocate(igpig,wfir)
! generate the overlap between the APW radial function and the KS
! second variational wfmt. These are the temporary Wannier projectors
! allocate the temporary wannier projectors
allocate(wantemp(ld,nst,nspinor,norb,nproj))
allocate(u(nrmtmax))
wantemp(:,:,:,:,:)=0.d0
do iorb=1,norb
  is=orb(iorb,1)
  l=orb(iorb,2)
  lmmax=2*l+1
  nr=nrmt(is)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    u(:)=0.d0
  !choose the local function for wannier projection to be the (first)
  !apw radial function
    u(1:nr)=apwfr(1:nr,1,1,l,ias)
    call wanolp(ik,ias,nst,l,lmmax,u,wfmt,wantemp(:,:,:,iorb,ia)) 
!transform temporary projector for equivalent atoms by symlm symmetry matrix
    allocate(z1(lmmax,nst))
    allocate(z2(lmmax,lmmax))
    z2(:,:)=symlm(1:lmmax,1:lmmax,iorb,ia)
    do ispn=1,nspinor
  !Transforming to equivalent atom
      call zgemm('N','N',lmmax,nst,lmmax,zone,z2,lmmax, &
              wantemp(1:lmmax,:,ispn,iorb,ia),lmmax,zzero,z1,lmmax)
      wantemp(1:lmmax,:,ispn,iorb,ia)=z1(:,:)
    enddo
    deallocate(z1,z2)
! rotate temporary projectors into different basis
    if (lmirep) then
      allocate(z2(lmmax,lmmax))
      z2(:,:)=subulm(1:lmmax,1:lmmax,iorb,ia)
      allocate(z(lmmax))
      do ispn=1,nspinor
        do ist=1,nst
          call zgemv('N',lmmax,lmmax,zone,z2,lmmax, &
               wantemp(1:lmmax,ist,ispn,iorb,ia),1,zzero,z,1)
          wantemp(1:lmmax,ist,ispn,iorb,ia)=z(:)
        end do
      end do
      deallocate(z,z2)
    end if
  enddo
enddo
deallocate(u,wfmt)
! calculate the "true" wannier projectors
wanprj(:,:,:,:,:)=0.d0 
call wannorm(nst,ld,nproj,sublm,wantemp, &
         wanprj(:,1:nst,:,:,:))
deallocate(wantemp)
!rotate back into spherical harmonics 
if (lmirep) then
  do iorb=1,norb
    is=orb(iorb,1)
!rotate projectors into spherical harmonic basis
    l=orb(iorb,2)
    lmmax=2*l+1
    do ia=1,natoms(is)
      allocate(z(lmmax),z2(lmmax,lmmax))
      z2(:,:)=subulm(1:lmmax,1:lmmax,iorb,ia)
      do ispn=1,nspinor
        do ist=1,nst
          call zgemv('C',lmmax,lmmax,zone,z2,lmmax, &
               wanprj(1:lmmax,ist,ispn,iorb,ia),1,zzero,z,1)
          wanprj(1:lmmax,ist,ispn,iorb,ia)=z(:)
        end do
      end do
      deallocate(z,z2)
    end do
  end do
endif
return
end subroutine
