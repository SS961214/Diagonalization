CC = nvcc
SRC_EXT = .cu
OUT_EXT = .out
CFLAGS  = -O2 -DADD_ -arch=compute_75
DEBAG_FLAG = -g
CFLAGS += $(DEBAG_FLAG)
LIBES = $(shell pkg-config --cflags --libs magma)
HEADER_DIR = ../Headers

ALL_O  = $(wildcard ${HEADER_DIR}/*.o)
ALL_H  = $(wildcard ${HEADER_DIR}/*.h)

TARGET  = testing_cheevd_GPU.out testing_zheevd_GPU.out
EXECS   = ${TARGET:.out=${OUT_EXT}}
SOURCES = ${TARGET:.out=${SRC_EXT}}

.SUFFIXES: .cu

all : ${FRC}
	${MAKE} "CFLAGS=${CFLAGS}" "FRC=${FRC}" ${EXECS}

Makefile : ${FRC}
	cp basic.mk $@
	@echo '# Automatically-generated dependencies list:' >> $@
	@${CC} ${CFLAGS} -MM ${SOURCES} >> $@

.cu.out:
	${CC} ${CFLAGS} -o $@ $< ${ALL_O} ${LIBES}

force_rebuild :
