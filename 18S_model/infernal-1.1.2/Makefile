# Top level Makefile for Infernal
# 
# On most systems, to build Infernal you should only need:
#     % ./configure; make
#
# Optionally, you can run a test suite:
#     % make check
#
# And optionally, you can install everything more permanently:
#     % make install
#
# For more information, see the Installation chapter of the Infernal
# User's Guide.
#


# VPATH and shell configuration
#
top_srcdir     = .
srcdir         = .

SHELL          = /bin/sh

# location of easel
ESLDIR         = easel

# location of hmmer
HMMERDIR       = hmmer

# location of suffix array library
SADIR          = hmmer/libdivsufsort

# impl directory (currently only either 'impl_sse' or '')
IMPLDIR        = impl_sse

# Package information
#
PACKAGE         = Infernal
PACKAGE_VERSION = 1.1.2
PACKAGE_TARNAME = infernal

# Installation targets
#
prefix      = /usr/local
exec_prefix = ${prefix}
datarootdir = ${prefix}/share
bindir      = ${exec_prefix}/bin
libdir      = ${exec_prefix}/lib
includedir  = ${prefix}/include
docdir      = ${datarootdir}/doc/${PACKAGE_TARNAME}
pdfdir      = ${docdir}
mandir      = ${datarootdir}/man
man1dir     = ${mandir}/man1
man1ext     = .1

# Compiler configuration
#
CC        = gcc
CFLAGS    = -g -O2
LDFLAGS   = -static 
SIMDFLAGS = -msse2 -msse3
CPPFLAGS  = 

# Other tools
#
AR        = /usr/bin/ar 
RANLIB    = ranlib
INSTALL   = /usr/bin/install -c

# beautification magic stolen from git 
#
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

# Programs installed by Infernal
#
PROGS = cmalign\
	cmbuild\
	cmcalibrate\
	cmemit\
	cmfetch\
	cmscan\
	cmsearch\
	cmstat

SSEPROGS = sse_cmsearch\
           cmsearch_dual

.PHONY: all dev check pdf install uninstall clean distclean TAGS

# all: Compile all documented executables.
#      (Excludes test programs.)
#
all: 
	${QUIET_SUBDIR0}${ESLDIR}        ${QUIET_SUBDIR1} all
	${QUIET_SUBDIR0}${SADIR}         ${QUIET_SUBDIR1} all
	${QUIET_SUBDIR0}${HMMERDIR}      ${QUIET_SUBDIR1} all
	${QUIET_SUBDIR0}src              ${QUIET_SUBDIR1} all
	${QUIET_SUBDIR0}rmark            ${QUIET_SUBDIR1} all
ifdef IMPLDIR
	${QUIET_SUBDIR0}src/impl_sse     ${QUIET_SUBDIR1} all
endif

# dev: compile all executables, including drivers.
#
dev: 
	${QUIET_SUBDIR0}${ESLDIR}       ${QUIET_SUBDIR1} dev
	${QUIET_SUBDIR0}${SADIR}        ${QUIET_SUBDIR1} all
	${QUIET_SUBDIR0}${HMMERDIR}     ${QUIET_SUBDIR1} dev
	${QUIET_SUBDIR0}src             ${QUIET_SUBDIR1} dev
	${QUIET_SUBDIR0}rmark           ${QUIET_SUBDIR1} all

# tests: compile all test drivers for 'make check'
#
tests:
	${QUIET_SUBDIR0}${ESLDIR}       ${QUIET_SUBDIR1} tests
	${QUIET_SUBDIR0}${SADIR}        ${QUIET_SUBDIR1} all
	${QUIET_SUBDIR0}${HMMERDIR}     ${QUIET_SUBDIR1} tests
	${QUIET_SUBDIR0}src             ${QUIET_SUBDIR1} tests

# check: Run test suites.
#
# We do a 'check' in HMMERDIR/testsuite, not HMMERDIR/ because we
# don't need to run the Easel testsuite twice.
#

# Also, HMMER's 'make check' requires its own easel subdir complete
# with compiled miniapps. This may or may not exist. In working copies
# of the trunk it will exist but won't be compiled (infernal/hmmer/
# will use infernal/easel, not infernal/hmmer/easel). In a released
# version infernal/hmmer/easel won't even exist since we delete it
# prior to creating a release.  So, we take care to have the calling
# of hmmer's testsuite work in either case. If hmmer/easel/ exists we
# copy the required miniapps into it and run 'make check' from
# hmmer/testsuite then delete them. If it doesn't, we create a symlink
# to infernal/easel as infernal/hmmer/easel then we run the testsuite
# and then we delete the symlink. Also, if srcdir/hmmer/easel doesn't
# exist, we symlink it to srcdir/easel because the hmmer testsuite
# needs easel/formats/BLOSUM62 which will only be in srcdir/easel
# (yuck!).

check:
	${QUIET_SUBDIR0}${ESLDIR}              ${QUIET_SUBDIR1} tests
	${QUIET_SUBDIR0}${SADIR}               ${QUIET_SUBDIR1} all
	${QUIET_SUBDIR0}${HMMERDIR}            ${QUIET_SUBDIR1} tests
	${QUIET_SUBDIR0}src                    ${QUIET_SUBDIR1} tests
	${QUIET_SUBDIR0}${ESLDIR}              ${QUIET_SUBDIR1} check
	if   test -d ${HMMERDIR}/easel && ! test -e  ${HMMERDIR}/easel/miniapps/esl-afetch;   then cp ${ESLDIR}/miniapps/esl-afetch   ${HMMERDIR}/easel/miniapps; fi
	if   test -d ${HMMERDIR}/easel && ! test -e  ${HMMERDIR}/easel/miniapps/esl-reformat; then cp ${ESLDIR}/miniapps/esl-reformat ${HMMERDIR}/easel/miniapps; fi
	if   test -d ${HMMERDIR}/easel && ! test -e  ${HMMERDIR}/easel/miniapps/esl-sfetch;   then cp ${ESLDIR}/miniapps/esl-sfetch   ${HMMERDIR}/easel/miniapps; fi
	if   test -d ${HMMERDIR}/easel && ! test -e  ${HMMERDIR}/easel/miniapps/esl-sfetch;   then cp ${ESLDIR}/miniapps/esl-shuffle  ${HMMERDIR}/easel/miniapps; fi
	if ! test -d ${HMMERDIR}/easel; then cd ${HMMERDIR}; ln -s ../${ESLDIR} .; cd ..; fi
	if ! test -d ${srcdir}/${HMMERDIR}/easel; then ln -s ${srcdir}/${ESLDIR} ${srcdir}/${HMMERDIR}/easel; fi
	${QUIET_SUBDIR0}${HMMERDIR}/testsuite  ${QUIET_SUBDIR1} check
	if   test -L ${HMMERDIR}/easel; then rm ${HMMERDIR}/easel; fi
	if   test -L ${srcdir}/${HMMERDIR}/easel; then rm ${srcdir}/${HMMERDIR}/easel; fi
	if   test -e ${HMMERDIR}/easel/miniapps/esl-afetch;   then rm ${HMMERDIR}/easel/miniapps/esl-afetch;   fi
	if   test -e ${HMMERDIR}/easel/miniapps/esl-reformat; then rm ${HMMERDIR}/easel/miniapps/esl-reformat; fi
	if   test -e ${HMMERDIR}/easel/miniapps/esl-sfetch;   then rm ${HMMERDIR}/easel/miniapps/esl-sfetch;   fi
	if   test -e ${HMMERDIR}/easel/miniapps/esl-shuffle;  then rm ${HMMERDIR}/easel/miniapps/esl-shuffle;  fi
	${QUIET_SUBDIR0}testsuite              ${QUIET_SUBDIR1} check

# devcheck: Run Infernal developer test suite (slower, but more comprehensive than 'check').
#
devcheck:
	${QUIET_SUBDIR0}testsuite  ${QUIET_SUBDIR1} devcheck

# pdf: compile the User Guides.
#
pdf:
	${QUIET_SUBDIR0}documentation ${QUIET_SUBDIR1} pdf

# install: installs the binaries in ${bindir}/
#          When man pages are done, will install man pages in MANDIR/man1/  (e.g. if MANSUFFIX is 1)
#          Creates these directories if they don't exist.
#          Prefix those paths with ${DESTDIR} (rarely used, usually null;
#          may be set on a make command line when building contrib RPMs).
install: 
	${INSTALL} -d ${DESTDIR}${bindir}
	${INSTALL} -d ${DESTDIR}${man1dir}
	${INSTALL} -d ${DESTDIR}${pdfdir}
# src 'install' moves into impl_sse and does 'install' (if IMPLDIR is defined (system has SSE))
	${QUIET_SUBDIR0}src           ${QUIET_SUBDIR1} install
	${QUIET_SUBDIR0}documentation ${QUIET_SUBDIR1} install

# uninstall: Reverses the steps of "make install".
#
uninstall: 
	${QUIET_SUBDIR0}src           ${QUIET_SUBDIR1} uninstall
	${QUIET_SUBDIR0}documentation ${QUIET_SUBDIR1} uninstall

# "make clean" removes almost everything except configuration files.
#
clean:
	${QUIET_SUBDIR0}src           ${QUIET_SUBDIR1} clean
	${QUIET_SUBDIR0}rmark         ${QUIET_SUBDIR1} clean
	${QUIET_SUBDIR0}testsuite     ${QUIET_SUBDIR1} clean
	${QUIET_SUBDIR0}documentation ${QUIET_SUBDIR1} clean
	${QUIET_SUBDIR0}${HMMERDIR}   ${QUIET_SUBDIR1} clean
	${QUIET_SUBDIR0}${ESLDIR}     ${QUIET_SUBDIR1} clean
	-rm -f *.o *~ Makefile.bak core TAGS gmon.out

# "make distclean" leaves a pristine source distribution.
#
distclean:
	${QUIET_SUBDIR0}src           ${QUIET_SUBDIR1} distclean
	${QUIET_SUBDIR0}rmark         ${QUIET_SUBDIR1} distclean
	${QUIET_SUBDIR0}testsuite     ${QUIET_SUBDIR1} distclean
	${QUIET_SUBDIR0}documentation ${QUIET_SUBDIR1} distclean
	${QUIET_SUBDIR0}${HMMERDIR}   ${QUIET_SUBDIR1} distclean
# no need to distclean in easel/ because distclean in hmmer/ already did that (and easel/Makefile has been deleted)
	-rm config.log config.status
	-rm -rf autom4te.cache
	-rm -f *.o *~ Makefile.bak core TAGS gmon.out
	-rm -f cscope.po.out cscope.out cscope.in.out cscope.files
	-rm -f release-notes/LICENSE.sh
	-rm -f src/impl
	-rm Makefile

# "make TAGS": emacs etags database, for searching the whole source tree.
TAGS:
	./makeTAGS.sh



