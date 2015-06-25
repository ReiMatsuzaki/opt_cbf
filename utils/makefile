
include ../local.mk

CXXFLAGS=${LIBPATH} -O3
OBJS=keys_values.o timer.o

utest: utest.o ${OBJS}
	${CXX} -o $@ ${CXXFLAGS} utest.o ${OBJS} ${LIBGTEST}
	./utest

clean:
	rm *.o
	rm utest
