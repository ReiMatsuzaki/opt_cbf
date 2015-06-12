# rcclsc
SRCDIR=${HOME}/src/project/opt_cbf
L2_DIR=${SRCDIR}/l2func
UTILS_DIR=${SRCDIR}/utils
LIBGTEST= -lpthread ${HOME}/local/src/gtest-1.7.0/make/gtest_main.a
LIBPATH=-L${HOME}/local/lib \
	-I${HOME}/local/src/eigen-eigen-10219c95fe65 \
	-I${L2_DIR} -I${UTILS_DIR}
INSTALL_PATH=${HOME}/bin

