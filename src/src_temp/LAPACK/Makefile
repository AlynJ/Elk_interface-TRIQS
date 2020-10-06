
AR = ar

include ../../make.inc

#-------------------------------------------------------------------------------
# Suffix Rules
#-------------------------------------------------------------------------------
.SUFFIXES: .o .f
.f.o:
	$(F77) $(F77_OPTS) -c $<

#-------------------------------------------------------------------------------
# File dependencies
#-------------------------------------------------------------------------------
SRC = \
 dbdsqr.f  dgebd2.f  dgebrd.f  dgelq2.f  dgelqf.f  dgeqr2.f  dgeqrf.f \
 dgesvd.f  dgesv.f   dgetrf2.f dgetrf.f  dgetri.f  dgetrs.f  disnan.f \
 dlabad.f  dlabrd.f  dlacpy.f  dladiv.f  dlae2.f   dlaebz.f  dlaed0.f \
 dlaed1.f  dlaed2.f  dlaed3.f  dlaed4.f  dlaed5.f  dlaed6.f  dlaed7.f \
 dlaed8.f  dlaed9.f  dlaeda.f  dlaev2.f  dlagtf.f  dlagts.f  dlaisnan.f \
 dlamch.f  dlamrg.f  dlange.f  dlansp.f  dlanst.f  dlansy.f  dlapy2.f \
 dlapy3.f  dlarfb.f  dlarf.f   dlarfg.f  dlarft.f  dlarnv.f  dlartg.f \
 dlaruv.f  dlas2.f   dlascl.f  dlaset.f  dlasq1.f  dlasq2.f  dlasq3.f \
 dlasq4.f  dlasq5.f  dlasq6.f  dlasr.f   dlasrt.f  dlassq.f  dlasv2.f \
 dlaswp.f  dlatrd.f  dopgtr.f  dopmtr.f  dorg2l.f  dorg2r.f  dorgbr.f \
 dorgl2.f  dorglq.f  dorgql.f  dorgqr.f  dorgtr.f  dorm2l.f  dorm2r.f \
 dormbr.f  dorml2.f  dormlq.f  dormql.f  dormqr.f  dormtr.f  dpotrf2.f \
 dpotrf.f  dspev.f   dspevx.f  dsptrd.f  dstebz.f  dstedc.f  dstein.f \
 dsteqr.f  dsterf.f  dsyev.f   dsyevx.f  dsygs2.f  dsygst.f  dsygvx.f \
 dsytd2.f  dsytrd.f  dtrti2.f  dtrtri.f  dzsum1.f  ieeeck.f  iladlc.f \
 iladlr.f  ilaenv.f  ilazlc.f  ilazlr.f  iparmq.f  izmax1.f  xerbla.f \
 zbdsqr.f  zgbsv.f   zgbtf2.f  zgbtrf.f  zgbtrs.f  zgebak.f  zgebal.f \
 zgebd2.f  zgebrd.f  zgees.f   zgeev.f   zgehd2.f  zgehrd.f  zgelq2.f \
 zgelqf.f  zgeqr2.f  zgeqrf.f  zgesvd.f  zgesv.f   zgetrf2.f zgetrf.f \
 zgetri.f  zgetrs.f  zheevd.f  zheev.f   zheevx.f  zhegs2.f  zhegst.f \
 zhegv.f   zhegvx.f  zhetd2.f  zhetrd.f  zhpevx.f  zhpgst.f  zhpgvx.f \
 zhptrd.f  zhseqr.f  zlabrd.f  zlacgv.f  zlacn2.f  zlacpy.f  zlacrm.f \
 zladiv.f  zlaed0.f  zlaed7.f  zlaed8.f  zlahqr.f  zlahr2.f  zlange.f \
 zlanhe.f  zlanhp.f  zlaqr0.f  zlaqr1.f  zlaqr2.f  zlaqr3.f  zlaqr4.f \
 zlaqr5.f  zlarfb.f  zlarf.f   zlarfg.f  zlarft.f  zlartg.f  zlascl.f \
 zlaset.f  zlasr.f   zlassq.f  zlaswp.f  zlatrd.f  zlatrs.f  zpotrf2.f \
 zpotrf.f  zpptrf.f  zrot.f    zstedc.f  zstein.f  zsteqr.f  ztrevc3.f \
 ztrexc.f  ztrsen.f  ztrsyl.f  ztrti2.f  ztrtri.f  zung2l.f  zung2r.f \
 zungbr.f  zunghr.f  zungl2.f  zunglq.f  zungql.f  zungqr.f  zungtr.f \
 zunm2l.f  zunm2r.f  zunmbr.f  zunmhr.f  zunml2.f  zunmlq.f  zunmql.f \
 zunmqr.f  zunmtr.f  zupgtr.f  zupmtr.f  iparam2stage.f

OBJ = $(SRC:.f=.o)

lapack:	$(OBJ)
	$(AR) -rc lapack.a $(OBJ)

clean:
	rm -f *.o *.mod *~ *.a ifc* *.gcno gmon.out

ls:
	ls -x --tabsize=0 --width=80 *.f
