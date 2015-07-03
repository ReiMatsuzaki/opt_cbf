include local.mk
DEBUGS=-g3 -Wall -O0
OPTS=-O3 -Wall
CXXFLAGS=${LIBPATH} ${OPTS}
OPT_CBF_OBJS= factory.o opt_cbf.o driv.o opt.o restrict.o l_algebra.o
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

utest: utest.o ${OPT_CBF_OBJS} ${UTILS_DIR}/keys_values.o
	${CXX} -o utest ${CXXFLAGS} utest.o ${OPT_CBF_OBJS} ${L2_DIR}/l2.a ${UTILS_DIR}/keys_values.o ${LIBGTEST}
	./utest

factory.o: factory.cpp factory.hpp
	${CXX} -o $@ -c ${CXXFLAGS} factory.cpp
utest_factory: utest_factory.o ${OPT_CBF_OBJS} ${UTILS_DIR}/keys_values.o
	${CXX} -o $@ ${CXXFLAGS} utest_factory.o ${OPT_CBF_OBJS} ${L2_DIR}/l2.a ${UTILS_DIR}/keys_values.o ${LIBGTEST}

opt_cbf: ${RUN_OBJS}
	${CXX} -o opt_cbf ${CXXFLAGS} ${RUN_OBJS}
	./opt_cbf samples/2sto/sample.in samples/2sto/sample.out
	cat samples/2sto/sample.out
clean:
	rm -f *.o
	rm -f *.a
	rm -f utest

cleanall:
	cd utils
	make clean
	cd ..
	cd l2func
	make clean
	cd ..

install: opt_cbf
	cp opt_cbf ${INSTALL_PATH}
