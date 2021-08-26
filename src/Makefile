CC=gcc

CFLAGS=-g

CFLAGSS=-std=c99 -Wall -DUNIX -g -DDBUG

SRC_DIR:=./

BIN_DIR:=../bin

BUILD_DIR:=../build

SRCS:=$(shell find $(SRC_DIR) -name '*.c')

# Check for 32-bit vs 64-bit
PROC_TYPE = $(strip $(shell uname -m | grep 64))
 
# Check for Mac OS
OS = $(shell uname -s 2>/dev/null | tr [:lower:] [:upper:])
DARWIN = $(strip $(findstring DARWIN, $(OS)))

# MacOS System
ifneq ($(DARWIN),)
	CFLAGS += -DMAC
	LIBS=-framework OpenCL

	ifeq ($(PROC_TYPE),)
		CFLAGS+=-arch i386
	else
		CFLAGS+=-arch x86_64
	endif
else

# Linux OS
LIBS=-lOpenCL
ifeq ($(PROC_TYPE),)
	CFLAGS+=-m32
else
	CFLAGS+=-m64
endif

# Check for Linux-AMD
ifdef AMDAPPSDKROOT
   INC_DIRS=. $(AMDAPPSDKROOT)/include
	ifeq ($(PROC_TYPE),)
		LIB_DIRS=$(AMDAPPSDKROOT)/lib/x86
	else
		LIB_DIRS=$(AMDAPPSDKROOT)/lib/x86_64
	endif
else

# Check for Linux-Nvidia
ifdef CUDA
   INC_DIRS=. $(CUDA)/OpenCL/common/inc
endif

endif
endif



TAR_TOOLS:=nd_gen
TAR_OCL:=launcher test

TARGETS=$(TAR_TOOLS) $(TAR_OCL)

all: $(TARGETS)

$(TAR_TOOLS):%: %.c.o lbm.c.o
	$(CC) $^ -o $@

$(TAR_OCL):%: %.c.o lbm.c.o lbmcl.c.o
	$(CC) $(CFLAGSS) -o $@ $^ $(INC_DIRS:%=-I%) $(LIB_DIRS:%=-L%) $(LIBS)

%.c.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

.PHONY: clean install

clean:
	rm -f *.c.o $(TARGETS)

install:
	cp $(TARGETS) -t $(BIN_DIR)
	cp *.cl -t $(BIN_DIR)
