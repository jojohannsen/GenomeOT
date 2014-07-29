COMPILER = gcc
DFLAGS = -ansi -pedantic -Wall -g
CFLAGS = -c
OFLAGS = -o
EXECNAME = suffixtree
CENTROMERE = centromere
CHRCOMPARE = chrcompare
ST_SCAN = st_scan

suffixtree:	main.o suffix_tree.o
	${COMPILER} ${DFLAGS} main.o suffix_tree.o ${OFLAGS} ${EXECNAME}

centromere:	centromere.o suffix_tree.o
	${COMPILER} ${DFLAGS} centromere.o suffix_tree.o ${OFLAGS} ${CENTROMERE}

chrcompare:	chrcompare.o suffix_tree.o
	${COMPILER} ${DFLAGS} chrcompare.o suffix_tree.o ${OFLAGS} ${CHRCOMPARE}

st_scan:	st_scan.o suffix_tree.o
	${COMPILER} ${DFLAGS} st_scan.o suffix_tree.o ${OFLAGS} ${ST_SCAN}

suffix_tree.o:	suffix_tree.c suffix_tree.h
	${COMPILER} ${DFLAGS} ${CFLAGS} suffix_tree.c

main.c:	suffix_tree.h 
	${COMPILER} ${DFLAGS} ${CFLAGS} main.c 

centromere.c: suffix_tree.h
	${COMPILER} ${DFLAGS} ${CFLAGS} centromere.c 

chrcompare.c: suffix_tree.h
	${COMPILER} ${DFLAGS} ${CFLAGS} chrcompare.c 

st_scan.c: suffix_tree.h
	${COMPILER} ${DFLAGS} ${CFLAGS} st_scan.c 

clean: 
	rm *.o 
	rm ${EXECNAME}
	rm ${CENTROMERE}
	rm ${CHRCOMPARE}

