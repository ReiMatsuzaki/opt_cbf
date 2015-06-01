# ubuntu in vb lab
#SRCDIR=${HOME}/source/project/opt_cbf
#L2_DIR=${SRCDIR}/l2func
#LIBGTEST= -lgtest_main -lgtest -lpthread
#LIBPATH=-I${HOME}/local/src/eigen-eigen-10219c95fe65\
#	-I${SRCDIR} -I${L2_DIR}	${LIBGTEST}

# rcclsc
#SRCDIR=${HOME}/src/project/opt_cbf_h
#L2_DIR=${SRCDIR}/l2func
#LIBGTEST= -lpthread ${HOME}/local/src/gtest-1.7.0/make/gtest_main.a
#LIBPATH=-L${HOME}/local/lib \
#	-I${HOME}/local/src/eigen-eigen-10219c95fe65 \
#	-I${L2_DIR}

include local.mk

CXXFLAGS=${LIBPATH} -Wall

utest: utest.o
	${CXX} -o $@ ${CXXFLAGS} utest.o
	./utest

clean:
	rm *.o
	rm utest




