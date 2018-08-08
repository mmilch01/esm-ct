SRCRT		= ~/src/mipl_light
RELEASE 	= ~/bin
PROG	= imfeat/imfeat

include mlib3.mk

CXXSRCS	= ${PROG}.cpp
OBJS	= ${CXXSRCS:.cpp=.o}

LIBOBJS	= ${NMOBJS} ${TRXOBJS} ${MLIB3OBJS}

CXX = c++
CXXFLAGS = -O2 -fpermissive -Wno-deprecated  -I $(TRXDIR) -I $(MLIB3DIR)

$(PROG): ${OBJS} ${LIBOBJS}
	$(CXX) $(CXXFLAGS) -o $@ ${OBJS} ${LIBOBJS} 
.c.o:
	gcc -c $(CXXFLAGS) $< -o $@
.cpp.o:
	$(CXX) $(CXXFLAGS) $< -o $@
clean:
	rm -f ${OBJS}
	rm -f ${PROG}
cleanall:
	rm -f ${OBJS}
	rm -f ${LIBOBJS}
	rm -f ${PROG}
release: ${PROG}
	chmod 771 ${PROG}
	cp -f ${PROG} ${RELEASE}
