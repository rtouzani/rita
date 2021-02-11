## Process this file with automake to produce Makefile.in

DDIR = $(datadir)/rita/doc

AUX_DIST = config.guess \
           config.sub \
           install-sh \
           ltconfig \
           ltmain.sh \
           mdate-sh \
           missing \
           mkinstalldirs

AUX_DIST_EXTRA = readline.m4 \
                 sys_errlist.m4 \
                 sys_siglist.m4

EXTRA_DIST = bootstrap ChangeLog
AUTOMAKE_OPTIONS = foreign
SUBDIRS = tutorial doc src
bin_SCRIPTS = test_rita.sh

MAINTAINERCLEANFILES = aclocal.m4 configure $(AUX_DIST)

DISTCLEANFILES =

ACLOCAL = aclocal

tests:
	chmod a+x test_rita.sh
	./test_rita.sh

clean-local:
	-rm -rf autom4te.cache *.log src/stamp-h1
