include ../local.mk
CXXFLAGS=${LIBPATH}

OBJS=hatom.o lcomb.o prim.o erfc.o fact.o

l2.a: ${OBJS}
	ar -cr -o $@ ${CXX_FLAGS} ${OBJS}

utest: utest.o l2.a
	${CXX} -o $@ ${CXX_FLAGS} ${LIBGTEST} utest.o l2.a
	./utest

clean:
	rm -f *.o
	rm -f *.a
	rm -f utest
