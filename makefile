include local.mk
CXXFLAGS=${LIBPATH}

utest: utest.o 
	${CXX} ${CXXFLAGS} utest.o
