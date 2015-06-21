include local.mk
CXXFLAGS=${LIBPATH} -g -Wall
OPT_CBF_OBJS= opt_cbf.o driv.o opt.o l_algebra.o
RUN_OBJS=run.o controller.o ${OPT_CBF_OBJS} ${UTILS_DIR}/keys_values.o ${UTILS_DIR}/timer.o ${L2_DIR}/l2.a

${UTILS_DIR}/keys_values.o :
	cd ${UTILS_DIR} \
	make keys_values.o \
	cd ..

${UTILS_DIR}/timer.o:
	cd ${UTILS_DIR} \
	make timer.o \
	cd ..

${L2_DIR}/l2.a:
	cd ${L2_DIR} \
	make l2.a \
	cd ..

utest: utest.o ${OPT_CBF_OBJS}
	${CXX} -o utest ${CXXFLAGS} utest.o ${OPT_CBF_OBJS} ${L2_DIR}/l2.a ${LIBGTEST}
	./utest

opt_cbf: ${RUN_OBJS}
	${CXX} -o opt_cbf ${CXXFLAGS} ${RUN_OBJS}
	./opt_cbf supply/sample.in supply/sample.out
	cat supply/sample.out
clean:
	rm -f *.o
	rm -f *.a
	rm -f utest

install: opt_cbf
	cp opt_cbf ${INSTALL_PATH}
