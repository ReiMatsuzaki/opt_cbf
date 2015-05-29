# macbook air
#SRCDIR=${HOME}/src/git/opt_cbf/l2func
#LIBGTEST=-lpthread ${HOME}/local/gtest-1.7.0/build/libgtest_main.a
#LIBPATH=-L${HOME}/local/lib -I${HOME}/local/gtest-1.7.0/include/  ${LIBGTEST}
#CXX=clang++ -std=c++11


# ubuntu (CAUSION: SRCDIR must be modify)
SRCDIR=${HOME}/src/git/opt_cbf/l2func
LIBGTEST= -lgtest_main -lgtest -lpthread
LIBPATH=-L${HOME}/local/lib

# rcclsc
#SRCDIR=${HOME}/src/project/l2func
#LIBGTEST= -lpthread ${HOME}/local/src/gtest-1.7.0/make/gtest_main.a
#LIBPATH=-L${HOME}/local/lib

CXXFLAGS=${LIBPATH}

utest: utest.o
	${CXX} -o $@ ${CXX_FLAGS} ${LIBGTEST} utest.o
	./utest
