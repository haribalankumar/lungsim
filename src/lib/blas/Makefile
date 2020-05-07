# -*- Makefile -*-

############################################################################
#
#  Program:         blas
#
#  Module:          Makefile
#
############################################################################

BUILDING_BLAS = true

include ../Makefile.inc

#-----------------------------------------------------------------------------
# End of initialisation: set targets and rules.
#-----------------------------------------------------------------------------

#
#  The name of the libraries to be created/linked to
#
BLASLIB   = $(LIBDIR)/libblas$(SUFFIX).a
LIB = $(BLASLIB)

#
#  The archiver and the flag(s) to use when building archive (library)
#  Convert between the rest of the worlds (and make's) default names
#  and those used in the SuperLU library.
#
ARCH         = $(AR)
ARCHFLAGS    = $(AFLAGS)
FORTRAN	     = $(FC)
LOADER       = $(CC)
LOADOPTS     =

#######################################################################
#  This is the makefile to create a library for the BLAS.
#  The files are grouped as follows:
#
#       SBLAS1 -- Single precision real BLAS routines
#       CBLAS1 -- Single precision complex BLAS routines
#       DBLAS1 -- Double precision real BLAS routines
#       ZBLAS1 -- Double precision complex BLAS routines
#
#       CB1AUX -- Real BLAS routines called by complex routines
#       ZB1AUX -- D.P. real BLAS routines called by d.p. complex
#                 routines
#
#      ALLBLAS -- Auxiliary routines for Level 2 and 3 BLAS
#
#       SBLAS2 -- Single precision real BLAS2 routines
#       CBLAS2 -- Single precision complex BLAS2 routines
#       DBLAS2 -- Double precision real BLAS2 routines
#       ZBLAS2 -- Double precision complex BLAS2 routines
#
#       SBLAS3 -- Single precision real BLAS3 routines
#       CBLAS3 -- Single precision complex BLAS3 routines
#       DBLAS3 -- Double precision real BLAS3 routines
#       ZBLAS3 -- Double precision complex BLAS3 routines
#
#  The library can be set up to include routines for any combination
#  of the four precisions.  To specify the objects to include in the build,
#  enter make variable OBJ to the objects desired.  Some examples:
#       make OBJ='$(SINGLEOBJ)'
#       make OBJ='$(ALLAUX) $(DZLAUX) $(DLASRC) $(ZLASRC)'
#       make OBJ='$(ALLOBJ)'
#  Alternatively, the command
#       make
#  without any arguments creates a library of just the double precision.
#
#  To remove the object files after the library is created, enter
#       make clean
#  To force the source files to be recompiled, enter, for example,
#       make single FRC=FRC
#
#---------------------------------------------------------------------
#
#  Edward Anderson, University of Tennessee
#  March 26, 1990
#  Susan Ostrouchov, Last updated September 30, 1994
#
#######################################################################

#---------------------------------------------------------
#  Comment out the next 6 definitions if you already have
#  the Level 1 BLAS.
#---------------------------------------------------------
SBLAS1 = isamax.o sasum.o saxpy.o scopy.o sdot.o snrm2.o \
	srot.o srotg.o sscal.o sswap.o

CBLAS1 = scasum.o scnrm2.o icamax.o caxpy.o ccopy.o \
	cdotc.o cdotu.o csscal.o crotg.o cscal.o cswap.o

DBLAS1 = idamax.o dasum.o daxpy.o dcopy.o ddot.o dnrm2.o \
	drot.o drotg.o dscal.o dswap.o

ZBLAS1 = dcabs1.o dzasum.o dznrm2.o izamax.o zaxpy.o zcopy.o \
	zdotc.o zdotu.o zdscal.o zrotg.o zscal.o zswap.o

CB1AUX = isamax.o sasum.o saxpy.o scopy.o snrm2.o sscal.o

ZB1AUX = idamax.o dasum.o daxpy.o dcopy.o dnrm2.o dscal.o

#---------------------------------------------------------------------
#  The following line defines auxiliary routines needed by both the
#  Level 2 and Level 3 BLAS.  Comment it out only if you already have
#  both the Level 2 and 3 BLAS.
#---------------------------------------------------------------------
ALLBLAS  = lsame.o xerbla.o

#---------------------------------------------------------
#  Comment out the next 4 definitions if you already have
#  the Level 2 BLAS.
#---------------------------------------------------------
SBLAS2 = sgemv.o sgbmv.o ssymv.o ssbmv.o sspmv.o \
	strmv.o stbmv.o stpmv.o strsv.o stbsv.o stpsv.o \
	sger.o ssyr.o sspr.o ssyr2.o sspr2.o

CBLAS2 = cgemv.o cgbmv.o chemv.o chbmv.o chpmv.o \
	ctrmv.o ctbmv.o ctpmv.o ctrsv.o ctbsv.o ctpsv.o \
	cgerc.o cgeru.o cher.o chpr.o cher2.o chpr2.o

DBLAS2 = dgemv.o dgbmv.o dsymv.o dsbmv.o dspmv.o \
	dtrmv.o dtbmv.o dtpmv.o dtrsv.o dtbsv.o dtpsv.o \
	dger.o dsyr.o dspr.o dsyr2.o dspr2.o

ZBLAS2 = zgemv.o zgbmv.o zhemv.o zhbmv.o zhpmv.o \
	ztrmv.o ztbmv.o ztpmv.o ztrsv.o ztbsv.o ztpsv.o \
	zgerc.o zgeru.o zher.o zhpr.o zher2.o zhpr2.o

#---------------------------------------------------------
#  Comment out the next 4 definitions if you already have
#  the Level 3 BLAS.
#---------------------------------------------------------
SBLAS3 = sgemm.o ssymm.o ssyrk.o ssyr2k.o strmm.o strsm.o 

CBLAS3 = cgemm.o csymm.o csyrk.o csyr2k.o ctrmm.o ctrsm.o \
	chemm.o cherk.o cher2k.o

DBLAS3 = dgemm.o dsymm.o dsyrk.o dsyr2k.o dtrmm.o dtrsm.o

ZBLAS3 = zgemm.o zsymm.o zsyrk.o zsyr2k.o ztrmm.o ztrsm.o \
	zhemm.o zherk.o zher2k.o


ALLOBJ = $(SBLAS1) $(DBLAS1) \
	$(CBLAS1) $(CB1AUX) $(ZBLAS1) $(ZB1AUX) \
	$(ALLBLAS) $(SBLAS2) $(SBLAS3) $(DBLAS2) $(DBLAS3) \
	$(CBLAS2) $(CBLAS3) $(ZBLAS2) $(ZBLAS3)
SINGLEOBJ = $(SBLAS1) $(ALLBLAS) $(SBLAS2) $(SBLAS3)
DOUBLEOBJ = $(DBLAS1) $(ALLBLAS) $(DBLAS2) $(DBLAS3)
COMPLEXOBJ = $(CBLAS1) $(CB1AUX) $(ALLBLAS) $(CBLAS2) $(CBLAS3)
COMPLEX16OBJ = $(ZBLAS1) $(ZB1AUX) $(ALLBLAS) $(ZBLAS2) $(ZBLAS3)

# Change OBJ if you want to include other precisions.
OBJ = $(DOUBLEOBJ)

# single complex double complex16 cannot be made in parallel because they all
# build into the same archive.
$(LIB): $(OBJ)
	$(ARCH) $(ARCHFLAGS) $(LIB) $(OBJ)
	$(RANLIB) $(LIB)

.INTERMEDIATE: $(OBJ)

$(OBJ): $(FRC)

FRC:
	@FRC=$(FRC)

clean:
	rm -rf $(CLEAN_OBJS)
