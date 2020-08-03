SRC_EXT = .c
OUT_EXT = .out
CFLAGS  = -std=c11 -Wall -O2
DEBAG_FLAG = -g
LIBES = -lm -llapack -lblas
HEADER_DIR = ../Headers

ALL_O  = $(wildcard ${HEADER_DIR}/*.o)
ALL_H  = $(wildcard ${HEADER_DIR}/*.h)

TARGET =
EXECS = ${TARGET:.out=${OUT_EXT}}
SOURCES = ${TARGET:.out=.c}

all : ${FRC}
	# ${MAKE} "CFLAGS=${CFLAGS}" "FRC=${FRC}" ${EXECS}

Makefile : ${FRC}
	cp basic.mk $@
	@echo '# Automatically-generated dependencies list:' >> $@
	# @${CC} ${CFLAGS} -MM ${SOURCES} >> $@

NONE :
	:

.o.out:
	${CC} ${CFLAGS} -o $@ $< ${ALL_O} ${LIBES}

force_rebuild :
