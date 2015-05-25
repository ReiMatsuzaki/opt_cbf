SRCDIR=${HOME}/src/project/opt_cbf_h
L2_DIR=${SRCDIR}/l2func
LIBGTEST= -lpthread ${HOME}/local/src/gtest-1.7.0/make/gtest_main.a
LIBPATH=-L${HOME}/local/lib \
	-I${HOME}/local/src/eigen-eigen-10219c95fe65 \
	-I${L2_DIR}
CXXFLAGS=${LIBPATH} -Wall

utest: utest.o
	${CXX} -o $@ ${CXX_FLAGS} ${LIBGTEST} utest.o
	./utest

clean:
	rm *.o
	rm utest




