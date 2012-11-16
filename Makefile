# Makefile -- Use this to build on *NIX systems.

# Options set on command line.
debug      = not-set
MPI        = not-set
verbose    = not-set

# This proxies everything to the builddir cmake.

cputype = $(shell uname -m | sed "s/\\ /_/g")
systype = $(shell uname -s)

BUILDDIR := build/$(systype)-$(cputype)
CONFIG_FLAGS = -DUNIX=1

# Process configuration options.

# Verbose builds?
ifeq ($(verbose), 1)
  CONFIG_FLAGS += -DCMAKE_VERBOSE_MAKEFILE=1
endif

# MPI
ifeq ($(MPI), 1)
  BUILDDIR := ${BUILDDIR}-MPI
  CC = mpicc
  CONFIG_FLAGS += -DUSE_MPI=1
else
#  ifeq ($(systype), Darwin)
#    CC  = clang
#    CXX = clang++
#  else
    CC  = cc
    CXX = c++
#  endif
  CONFIG_FLAGS += -DUSE_MPI=0
endif

BUILDDIR := ${BUILDDIR}-${CC}
CONFIG_FLAGS += -DCC=${CC} -DCXX=${CXX}

# Debugging symbols
ifneq ($(debug), not-set)
  BUILDDIR := ${BUILDDIR}-Debug
  CONFIG_FLAGS += -DCMAKE_BUILD_TYPE=Debug
else
  BUILDDIR := ${BUILDDIR}-Release
  CONFIG_FLAGS += -DCMAKE_BUILD_TYPE=Release
endif

# Special considerations for specific systems.
ifeq ($(systype), Darwin)
  CONFIG_FLAGS += -DAPPLE=1
else 
  ifeq ($(systype), Linux)
    CONFIG_FLAGS += -DLINUX=1
  endif
endif

define run-config
mkdir -p $(BUILDDIR)
cd $(BUILDDIR) && cmake $(CURDIR) $(CONFIG_FLAGS)
endef

all test clean install:
	@if [ ! -f $(BUILDDIR)/Makefile ]; then \
		more INSTALL; \
	else \
		make -C $(BUILDDIR) $@ $(MAKEFLAGS); \
	fi

config: distclean
	$(run-config)

distclean:
	rm -rf $(BUILDDIR)

#dist:
#	utils/mkdist.sh $(PKGNAME)

.PHONY: config distclean all clean install uninstall 
