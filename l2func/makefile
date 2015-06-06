# macbook air
#SRCDIR=${HOME}/src/git/opt_cbf/l2func
#LIBGTEST=-lpthread ${HOME}/local/gtest-1.7.0/build/libgtest_main.a
#LIBPATH=-L${HOME}/local/lib -I${HOME}/local/gtest-1.7.0/include/  ${LIBGTEST}
#CXX=clang++ -std=c++11


# ubuntu (CAUSION: SRCDIR must be modify)
#SRCDIR=${HOME}/src/git/opt_cbf/l2func
#LIBGTEST= -lgtest_main -lgtest -lpthread
#LIBPATH=-L${HOME}/local/lib

# rcclsc
#SRCDIR=${HOME}/src/project/l2func
#LIBGTEST= -lpthread ${HOME}/local/src/gtest-1.7.0/make/gtest_main.a
#LIBPATH=-L${HOME}/local/lib

include ../local.mk
CXXFLAGS=${LIBPATH}
OBJS=utest.o lcomb.o prim.o erfc.o fact.o

l2.a: lcomb.o prim.o erfc.o fact.o
	ar -cr -o $@ ${CXX_FLAGS} lcomb.o prim.o erfc.o fact.o

utest: utest.o l2.a
	${CXX} -o $@ ${CXX_FLAGS} ${LIBGTEST} utest.o l2.a
	./utest
