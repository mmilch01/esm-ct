#extenstion convention: .mk - library package; .mak - program package.
MLIB3DIR 	= mlib3
NEWMDIR 	= mlib3/newm
TRXDIR 		= mlib3/TRX

MLIB3SRC = ${MLIB3DIR}/mlib3analyze.cpp ${MLIB3DIR}/mlib3common.cpp ${MLIB3DIR}/mlib3math.cpp ${MLIB3DIR}/mlib3volume.cpp ${MLIB3DIR}/anyoption.cpp
MLIB3OBJS= ${MLIB3SRC:.cpp=.o}
NEWMSRCS= ${NEWMDIR}/newmatrm.cpp ${NEWMDIR}/myexcept.cpp ${NEWMDIR}/newmat1.cpp ${NEWMDIR}/newmat2.cpp ${NEWMDIR}/newmat3.cpp ${NEWMDIR}/newmat4.cpp ${NEWMDIR}/newmat5.cpp ${NEWMDIR}/newmat6.cpp ${NEWMDIR}/newmat7.cpp ${NEWMDIR}/newmat8.cpp ${NEWMDIR}/newmatex.cpp ${NEWMDIR}/bandmat.cpp ${NEWMDIR}/submat.cpp ${NEWMDIR}/cholesky.cpp ${NEWMDIR}/evalue.cpp ${NEWMDIR}/fft.cpp ${NEWMDIR}/hholder.cpp ${NEWMDIR}/jacobi.cpp ${NEWMDIR}/newfft.cpp ${NEWMDIR}/sort.cpp ${NEWMDIR}/svd.cpp ${NEWMDIR}/nm_misc.cpp
NMOBJS	= ${NEWMSRCS:.cpp=.o}
TRXSRCS	= ${TRXDIR}/Getifh.c ${TRXDIR}/endianio.c ${TRXDIR}/rec.c
TRXOBJS	= ${TRXSRCS:.c=.o}
