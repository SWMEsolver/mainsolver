mytests.o: mytests.f90\
        finitevolume_vars.o

nonlinearsolver.o: nonlinearsolver.f90\
        finitevolume_vars.o\
        momentmodels.o\
        physicsframe.o\
        mytests.o

momentmodels.o: momentmodels.f90\
        finitevolume_vars.o\
        reconstruction.o\
        physicsframe.o\
        mytests.o

physicsframe.o: physicsframe.f90 \
    	finitevolume_vars.o \
	reconstruction.o\
        mytests.o

finitevolume_vars.o: finitevolume_vars.f90

implicitSource.o: implicitSource.f90 \
        finitevolume_vars.o

finitevolume.o: finitevolume.f90 \
	physicsframe.o \
	finitevolume_vars.o \
	reconstruction.o \
 	shocksindicator.o \
	mesh.o \
	momentmodels.o\
        mytests.o

mesh.o: mesh.f90 \
	finitevolume_vars.o\
        mytests.o

output.o: output.f90 \
	physicsframe.o \
	finitevolume_vars.o\
        mytests.o

parameters.o: parameters.f90 \
	finitevolume_vars.o\
	momentmodels.o\
        finitevolume.o\
        timediscretization.o\
        mytests.o
   
reconstruction.o: reconstruction.f90 \
	finitevolume_vars.o\
        mytests.o

shocksindicator.o: shocksindicator.f90 \
	finitevolume_vars.o\
        mytests.o

timediscretization.o: timediscretization.f90 \
	physicsframe.o \
	finitevolume.o \
	finitevolume_vars.o \
	momentmodels.o \
	output.o\
        mytests.o

main.o: main.f90 \
	finitevolume.o \
	mesh.o \
	parameters.o \
	timediscretization.o \
	implicitSource.o\
	nonlinearsolver.o\
        mytests.o
