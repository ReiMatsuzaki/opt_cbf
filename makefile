include local.mk
CXXFLAGS=${LIBPATH} -g -Wall
OPT_CBF_OBJS= opt_cbf.o driv.o opt.o l_algebra.o
RUN_OBJS=run.o controller.o ${OPT_CBF_OBJS} ${UTILS_DIR}/keys_values.o ${UTILS_DIR}/timer.o ${L2_DIR}/l2.a

utest: utest.o ${OPT_CBF_OBJS}
	${CXX} -o utest ${CXXFLAGS} utest.o ${OPT_CBF_OBJS} ${L2_DIR}/l2.a ${LIBGTEST}
	./utest

opt_cbf: ${RUN_OBJS}
	${CXX} -o opt_cbf ${CXXFLAGS} ${RUN_OBJS}
	./opt_cbf
	cat supply/sample.out
clean:
	rm -f *.o
	rm -f *.a
	rm -f utest
