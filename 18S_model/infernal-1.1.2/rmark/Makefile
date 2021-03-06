top_srcdir = ..
srcdir     = .


CC       = gcc
CFLAGS   = -g -O2 
SIMDFLAGS= -msse2 -msse3
CPPFLAGS = 
LDFLAGS  = 
DEFS     = -DHAVE_CONFIG_H
LIBS     = -lhmmer -leasel    -lm
MPILIBS  = 
IMPLDIR  = impl_sse

SRCDIR   = ../src
ESLDIR   = easel
HMMERDIR = hmmer

ESLINC   = -I../${ESLDIR} -I${top_srcdir}/easel
SRCINC   = -I${SRCDIR} -I${top_srcdir}/src

PROGS    = rmark-create\
	   rmark-rocplot

PROGOBJS  =\
	rmark-create.o\
	rmark-rocplot.o

# beautification magic stolen from git 
QUIET_SUBDIR0 = +${MAKE} -C #space separator after -c
QUIET_SUBDIR1 = 
ifndef V
	QUIET_CC      = @echo '    ' CC $@;
	QUIET_GEN     = @echo '    ' GEN $@;
	QUIET_AR      = @echo '    ' AR $@;
	QUIET_SUBDIR0 = +@subdir=
	QUIET_SUBDIR1 = ; echo '    ' SUBDIR $$subdir; \
		        ${MAKE} -C $$subdir
endif

.PHONY: all dev distclean clean

all:    ${PROGS}
dev:    ${PROGS}

${PROGS}: %: %.o  ../${HMMERDIR}/src/libhmmer.a ../${ESLDIR}/libeasel.a 
	${QUIET_GEN}${CC} ${CFLAGS} ${DEFS} ${LDFLAGS} -L../${HMMERDIR}/src -L../${ESLDIR} -L${SRCDIR} -o $@ $@.o ${LIBS}

${PROGOBJS}: ../${HMMERDIR}/src/hmmer.h ../${HMMERDIR}/src/p7_config.h ../src/infernal.h ../src/config.h

.c.o: 
	${QUIET_CC}${CC} ${ESLINC} ${SRCINC} ${CFLAGS} ${SIMDFLAGS} ${DEFS} -o $@ -c $<		

distclean: clean
	-rm -f Makefile

clean:
	-rm -f *.o *~  ${PROGS} TAGS
	-rm -f *.gcno

################################################################
# Infernal - inference of RNA secondary structure alignments
# Version 1.1.2; July 2016
# Copyright (C) 2016 Howard Hughes Medical Institute.
# Other copyrights also apply. See the COPYRIGHT file for a full list.
# 
# Infernal is distributed under the terms of the GNU General Public License
# (GPLv3). See the LICENSE file for details.
################################################################

