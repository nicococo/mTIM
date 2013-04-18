MATLAB_DIR = /usr/local

MEX = $(MATLAB_DIR)/MATLAB-R2011b/bin/mex

MATLAB_INCL = $(MATLAB_DIR)/MATLAB-R2011b/extern/include

hinge: 
	${MEX} loss_hinge_native.cpp

sqrhinge: 
	${MEX} loss_sqrhinge_native.cpp

logistic: 
	${MEX} loss_logistic_native.cpp


all: hinge logistic sqrhinge


clean:	
	rm -f *.o *.mexa64 *.mex
