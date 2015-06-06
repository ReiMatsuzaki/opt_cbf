include local.mk
CXXFLAGS=${LIBPATH}
OPT_CBF_OBJS= opt.o l_algebra.o
RUN_OBJS=run.o controller.o ${OPT_CBF_OBJS} ${UTILS_DIR}/keys_values.o ${L2_DIR}/l2.a

utest: utest.o ${OPT_CBF_OBJS}
	${CXX} -o utest ${CXXFLAGS} utest.o ${OPT_CBF_OBJS} ${L2_DIR}/l2.a ${LIBGTEST}
	./utest

opt_cbf: ${OBJS}
	${CXX} -o opt_cbf ${CXXFLAGS} ${OBJS}

