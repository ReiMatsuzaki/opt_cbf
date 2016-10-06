include local.mk
DEBUGS=-g3 -Wall -O0
OPTS=-O3 -Wall
CXXFLAGS=${LIBPATH} ${OPTS}
OPT_CBF_OBJS= opt_cbf.o opt.o restrict.o l_algebra.o
RUN_OBJS=run.o controller.o factory.o ${OPT_CBF_OBJS} ${UTILS_DIR}/keys_values.o ${UTILS_DIR}/timer.o ${L2_DIR}/l2.a

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

utest: utest.o ${OPT_CBF_OBJS} ${UTILS_DIR}/keys_values.o ${UTILS_DIR}/timer.o
	${CXX} -o utest ${CXXFLAGS} utest.o ${OPT_CBF_OBJS} ${L2_DIR}/l2.a ${UTILS_DIR}/keys_values.o ${UTILS_DIR}/timer.o ${LIBGTEST}
	./utest
test: test.o  opt.o restrict.o l_algebra.o
	${CXX} -o test  test.o opt.o restrict.o l_algebra.o ${L2_DIR}/l2.a ${CXXFLAGS} -lgtest

factory.o: factory.cpp factory.hpp
	${CXX} -o $@ -c ${CXXFLAGS} factory.cpp
utest_factory: utest_factory.o factory.o ${OPT_CBF_OBJS} ${UTILS_DIR}/keys_values.o
	${CXX} -o $@ ${CXXFLAGS} utest_factory.o factory.o ${OPT_CBF_OBJS} ${L2_DIR}/l2.a ${UTILS_DIR}/keys_values.o ${LIBGTEST}

opt_cbf: ${RUN_OBJS}
	${CXX} -o opt_cbf ${CXXFLAGS} ${RUN_OBJS}
	./opt_cbf samples/2sto/sample.in samples/2sto/sample.out
	cat samples/2sto/sample.out

run_1skp_l_sto: run_1skp_l_sto.o opt.o restrict.o l_algebra.o
	${CXX} -o $@ run_1skp_l_sto.o opt.o restrict.o l_algebra.o  ${L2_DIR}/l2.a ${CXXFLAGS} -lgtest
run_1skp_l_sto_2basis: run_1skp_l_sto
	./run_1skp_l_sto --w 1.1 --zs 0.8-0.1j,0.4-0.6j

run_delta_1skpl_cutsto: run_delta_1skpl_cutsto.o opt.o restrict.o l_algebra.o
	${CXX} -o $@ run_delta_1skpl_cutsto.o opt.o restrict.o l_algebra.o  ${L2_DIR}/l2.a ${CXXFLAGS} -lgtest
run_delta_1skpl_cutsto_1basis: run_delta_1skpl_cutsto
	./run_delta_1skpl_cutsto --w 1.1 --zs 0.3-1.1j --r0 10.0

run_delta_1skpl_cutsto_2basis: run_delta_1skpl_cutsto
	./run_delta_1skpl_cutsto --w 1.1 --zs 0.3-1.1j,0.4-0.7j --r0 10.0 --maxit 1000 --file_psi psi.dat

.PHONY: check
check: test
	./test

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
