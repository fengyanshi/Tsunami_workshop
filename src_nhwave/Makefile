#-----------BEGIN MAKEFILE-------------------------------------------------
         SHELL         = /bin/sh
         DEF_FLAGS     = -P -C -traditional 
         EXEC          = nhwave
#==========================================================================
#--------------------------------------------------------------------------
#        PRECISION          DEFAULT PRECISION: SINGLE                     
#                           UNCOMMENT TO SELECT DOUBLE PRECISION
#--------------------------------------------------------------------------
         FLAG_1 = -DDOUBLE_PRECISION
         FLAG_2 = -DPARALLEL
         FLAG_3 = -DLANDSLIDE
#         FLAG_4 = -DBUBBLE
#         FLAG_5 = -DSEDIMENT
#         FLAG_6 = -DINTEL
#--------------------------------------------------------------------------
#  mpi defs 
#--------------------------------------------------------------------------
         CPP      = /usr/bin/cpp 
         CPPFLAGS = $(DEF_FLAGS)
#         FC       = gfortran
         FC        = mpif90
         DEBFLGS  = 
         OPT      = -O3
         CLIB     = 
#==========================================================================
         FFLAGS = $(DEBFLGS) $(OPT) 
         MDEPFLAGS = --cpp --fext=f90 --file=-
         RANLIB = ranlib
#--------------------------------------------------------------------------
#  CAT Preprocessing Flags
#--------------------------------------------------------------------------
         CPPARGS = $(CPPFLAGS) $(DEF_FLAGS) $(FLAG_1) $(FLAG_2) $(FLAG_3) \
                   $(FLAG_4) $(FLAG_5) $(FLAG_6) 
#--------------------------------------------------------------------------
#  Libraries           
#--------------------------------------------------------------------------
#         LIBS  = -L/Users/fengyanshi15/Public/hypre-2.8.0b/src/hypre/lib -lHYPRE
#         INCS  = -I/Users/fengyanshi15/Public/hypre-2.8.0b/src/hypre/include
         LIBS  = -L/home/1047/hypre-2.8.0b/src/hypre/lib -lHYPRE
         INCS  = -I/home/1047/hypre-2.8.0b/src/hypre/include
#--------------------------------------------------------------------------
#  Preprocessing and Compilation Directives
#--------------------------------------------------------------------------
.SUFFIXES: .o .f90 .F .F90 

.F.o:
	$(CPP) $(CPPARGS) $*.F > $*.f90
	$(FC)  -c $(FFLAGS) $(INCS) $*.f90
#	\rm $*.f90
#--------------------------------------------------------------------------
#  NHWAVE Source Code.
#--------------------------------------------------------------------------

MODS  = mod_global.F	mod_util.F   slide.F

MAIN  = nhwave.F 

SRCS = $(MODS)  $(MAIN)

OBJS = $(SRCS:.F=.o) nspcg.o

#--------------------------------------------------------------------------
#  Linking Directives               
#--------------------------------------------------------------------------

$(EXEC):	$(OBJS)
		$(FC) $(FFLAGS) $(LDFLAGS) -o $(EXEC) $(OBJS) $(LIBS)
#--------------------------------------------------------------------------
#  Cleaning targets.
#--------------------------------------------------------------------------

clean:
		/bin/rm -f *.o *.mod

clobber:	clean
		/bin/rm -f *.f90 *.o nhwave







