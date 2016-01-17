
# The C++ complier to use with options
CPP=g++ -Wall -g -O0

#
# Note: It is not necessary to compile this Makefile. This only builds the
# tests which are not necessary to use the Library. To build the tests,
# simply set GTEST_INCL and LIBGTEST to point to the right places as shown
# in their example values.
#
# Path to the gtest include source (e.g. gtest-1.6.0/include/)
# XXX TODO FIXME Set this path
 GTEST_INCL=gtest-1.6.0/include/

# Path to libgtest.a
# XXX TODO FIXME Set this path
# LIBGTEST=gtest-1.6.0/lib/.libs/libgtest.a

# Path to the hazy template library, (e.g. hazytl/include)
HTL_INCL=../hazytl/include

ifndef GTEST_INCL
$(error "Could not find google tests library -- See Makefile")
endif

# All binaries
ALL=bin/basic-test bin/run_tests

UNAME=$(shell uname)
ifneq ($(UNAME), Darwin)
	LIB_RT=-lrt
endif

all: $(ALL)

bin/basic-test:
	$(CPP) src/test/test_hogwild.cc -o bin/basic-test -I$(HTL_INCL) -Iinclude/ $(LIBGTEST) -lpthread $(LIB_RT)

bin/run_tests:
	$(CPP) src/test/test.cc -o bin/run_tests -I$(GTEST_INCL) -I$(HTL_INCL) -Iinclude/ $(LIBGTEST) -lpthread $(LIB_RT)


clean:
	rm -f $(ALL)
