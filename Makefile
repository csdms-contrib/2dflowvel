
GFORTRAN=gfortran
FFLAGS=
PROG=2dflowvel
VERSION=0.1
SOURCES=2dflowvel.f
OBJS=${SOURCES:.f=.o}

all: 2dflowvel

2dflowvel: ${SOURCES}
	${GFORTRAN} ${FFLAGS} -o 2dflowvel ${SOURCES}

dist:
	@mkdir ${PROG}-${VERSION}
	@cp ${SOURCES} ${PROG}-${VERSION}
	@tar cvfz ${PROG}-${VERSION}.tar.gz ${PROG}-${VERSION} 
	@rm -rf ${PROG}-${VERSION}

clean:
	@rm -f ${PROG} ${OBJS} core

