include local.mk
CXXFLAGS=${LIBPATH}
OPT_CBF_OBJS=run.o controller.o ${UTILS_DIR}/keys_values.o

utest: utest.o 
	${CXX} -o utest ${CXXFLAGS} utest.o ${LIBGTEST}
	./utest

opt_cbf: ${OPT_CBF_OBJS}
	${CXX} -o opt_cbf ${CXXFLAGS} ${OPT_CBF_OBJS}

clean:
	rm *.o
	rm utest
