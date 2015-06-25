include ../local.mk
CXXFLAGS=${LIBPATH} -O3

OBJS=hatom.o lcomb.o prim.o erfc.o fact.o

l2.a: ${OBJS}
	ar -cr -o $@ ${OBJS}

utest.o: utest.cpp
	${CXX} -c -o $@ ${CXXFLAGS} utest.cpp

utest: utest.o l2.a
	${CXX} -o $@ ${CXXFLAGS} ${LIBGTEST} utest.o l2.a
	./utest

clean:
	rm -f *.o
	rm -f *.a
	rm -f utest
