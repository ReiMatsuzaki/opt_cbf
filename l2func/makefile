include ../local.mk
CXXFLAGS=${LIBPATH} -Wall

OBJS=hatom.o lcomb.o prim.o erfc.o fact.o

# this command run correctly in rcclsc
#l2.a: ${OBJS}
#	ar -cr -o $@ ${OBJS}

l2.a: ${OBJS}
	ar r $@ ${OBJS}


utest.o: utest.cpp
	${CXX} -c -o $@ ${CXXFLAGS} utest.cpp

utest: utest.o l2.a
	${CXX} -o $@ ${CXXFLAGS} ${LIBGTEST} utest.o l2.a
	./utest

clean:
	rm -f *.o
	rm -f *.a
	rm -f utest
