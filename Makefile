# C++ Compiler and options
CPP=g++-5 -Wall -ggdb3 -O3 -march=native

# Path to hogwild (e.g. hogwildtl/include)
HOG_INCL=hogwildtl/include
# Path to Hazy Template Library (e.g. hazytl/include)
HTL_INCL=hazytl/include

LIBS=-lpthread -lnuma

# Conversion tools
TOOLS=bin/convert_matlab bin/convert bin/unconvert
UNAME=$(shell uname)
ifneq ($(UNAME), Darwin)
	LIB_RT=-lrt
endif

ALL= $(TOOLS) obj/frontend.o bin/svm bin/numasvm

all: $(ALL)

obj/frontend.o:
	$(CPP) -c src/frontend_util.cc -o obj/frontend.o

bin/svm: obj/frontend.o
	$(CPP) -o bin/svm src/svm_main.cc -I$(HOG_INCL) -I$(HTL_INCL) $(LIBS) $(LIB_RT) \
		obj/frontend.o

bin/numasvm: obj/frontend.o
	$(CPP) -o bin/numasvm src/numasvm_main.cc -I$(HOG_INCL) -I$(HTL_INCL) $(LIBS) $(LIB_RT) \
		obj/frontend.o

bin/bbsvm: obj/frontend.o
	$(CPP) -o bin/bbsvm src/bbsvm_main.cc -I$(HOG_INCL) -I$(HTL_INCL) $(LIBS) $(LIB_RT) \
		obj/frontend.o

bin/tracenorm: obj/frontend.o
	$(CPP) -o bin/tracenorm src/tracenorm.cc -I$(HOG_INCL) -I$(HTL_INCL) $(LIBS) $(LIB_RT) \
		obj/frontend.o

bin/predict: 
	$(CPP) -o bin/predict src/tracenorm/predict.cc -I$(HOG_INCL) -I$(HTL_INCL) $(LIBS) $(LIB_RT) \
		obj/frontend.o

bin/bbtracenorm: obj/frontend.o
	$(CPP) -o bin/bbtracenorm src/bbtracenorm.cc -I$(HOG_INCL) -I$(HTL_INCL) $(LIBS) $(LIB_RT) \
		obj/frontend.o

bin/multicut: obj/frontend.o
	$(CPP) -o bin/multicut src/multicut.cc -I$(HOG_INCL) -I$(HTL_INCL) $(LIBS) $(LIB_RT) \
		obj/frontend.o

bin/bbmulticut: obj/frontend.o
	$(CPP) -o bin/bbmulticut src/bbmulticut.cc -I$(HOG_INCL) -I$(HTL_INCL) $(LIBS) $(LIB_RT) \
		obj/frontend.o

bin/convert: src/tools/tobinary.cc
	$(CPP) -o bin/convert src/tools/tobinary.cc -I$(HOG_INCL) -I$(HTL_INCL) 

bin/convert_matlab: src/tools/tobinary.cc
	$(CPP) -o bin/convert_matlab src/tools/tobinary.cc -I$(HOG_INCL) -I$(HTL_INCL) -DMATLAB_CONVERT_OFFSET=1

bin/unconvert: src/tools/unconvert.cc
	$(CPP) -o bin/unconvert src/tools/unconvert.cc -I$(HOG_INCL) -I$(HTL_INCL) 

clean:
	rm -f $(ALL)
