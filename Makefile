#
# Date:      2011/07/11 17:55
# Author:    Jan Faigl
#

CXX:=ccache $(CXX)

include Mk/libs.mk

CPPFLAGS+=$(LOCAL_CFLAGS)
LDFLAGS+=$(LOCAL_LDFLAGS)

CPPFLAGS+=$(CRL_CFLAGS) $(LOG4CXX_CPPFLAGS)  $(BOOST_CFLAGS) $(CAIRO_CFLAGS) 
LDFLAGS+=$(CRL-ALGORITHM) $(CRL-GUI_LDFLAGS) $(CRL_LDFLAGS) $(CAIRO_LDFLAGS) $(BOOST_LDFLAGS) $(LOG4CXX_LDFLAGS)

CXXFLAGS+=-std=c++11
CXXFLAGS+=-O2 -march=native
CXXFLAGS+=-g

OBJS=\
     src/tgsoa.o\
     src/gsoa.o\
     src/gsoa_ring.o\
     src/gsoa_learning.o\
     src/route_path_utils.o\
     src/coords.o\
     src/glkh.o

TARGET=tgsoa

include Mk/comrob.mk

# clean rebuild 4 threads
x:
	make clean
	make -j4