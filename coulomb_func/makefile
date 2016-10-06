coulomb_bind.so: coulomb.cpp 
	${CXX} -I`python -c 'from distutils.sysconfig import *; print get_python_inc()'` -DPIC -bundle -fPIC -o $@ coulomb.cpp -lgsl -lblas ${CXXFLAGS} -lboost_python  -framework Python

test: coulomb_bind.so
	python test.py
install: coulomb_bind.so coulomb.py
	cp coulomb.py ~/local/Python/
	cp coulomb_bind.so ~/local/Python/
