
include ../local.mk

CXXFLAGS=${LIBPATH} -Wall

utest: utest.o keys_values.o
	${CXX} -o $@ ${CXXFLAGS} utest.o keys_values.o
	./utest

clean:
	rm *.o
	rm utest
