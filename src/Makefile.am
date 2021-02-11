AUTOMAKE_OPTIONS = no-dependencies
AM_CPPFLAGS = -std=c++1y

bin_PROGRAMS =	rita

rita_SOURCES = rita.cpp \
               rita.h \
               ritaException.h \
               approximation.cpp \
               approximation.h \
               cmd.cpp \
               cmd.h \
               configure.cpp \
               configure.h \
               data.cpp \
               data.h \
               eigen.cpp \
               eigen.h \
               equa.cpp \
               equa.h \
               help.h \
               integration.cpp \
               integration.h \
               mesh.cpp \
               mesh.h \
               optim.cpp \
               optim.h \
               runAE.cpp \
               runODE.cpp \
               runPDE.cpp \
               solve.cpp \
               solve.h \
               stationary.cpp \
               stationary.h \
               transient.cpp \
               transient.h

clean-local:
	-rm -f stamp-h1
