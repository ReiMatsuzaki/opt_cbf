SRCDIR=${HOME}/src/project/l2func
LIBGTEST= -lpthread ${HOME}/local/src/gtest-1.7.0/make/gtest_main.a
LIBPATH=-L${HOME}/local/lib
CXXFLAGS=${LIBPATH}

utest: utest.o
	${CXX} -o $@ ${CXX_FLAGS} ${LIBGTEST} utest.o
	./utest
