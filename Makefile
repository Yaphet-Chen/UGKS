SHELL = /bin/bash -O extglob

#--------------------------------------------------
#variables
#--------------------------------------------------
#compilation flag
OMP = yes
# FC = ifort
FC = gfortran
FCFLAGS = -module $(BIN) -O3
OMPFLAG = -openmp -parallel -fpp

ifeq ($(FC),gfortran)
    # FCFLAGS = -J$(BIN) -O3 -ffree-line-length-none
	FCFLAGS = -J$(BIN) -O3 -ffree-line-length-none -march=native
    OMPFLAG = -fopenmp
	# gfortran -ffree-line-length-none -march=native -funroll-loops -flto -pipe -O3
endif

ifeq ($(OMP),yes)
    FCFLAGS += $(OMPFLAG)
endif

#directory for executables
BIN = bin

#--------------------------------------------------
#compiling
#--------------------------------------------------
all: checkdir 1D 2D

#mkdir
checkdir:
	mkdir -p $(BIN)

#build executables
1D : checkdir
	$(FC) $(FCFLAGS) -o $(BIN)/UGKS1D src/UGKS1D.f90

2D : checkdir
	$(FC) $(FCFLAGS) -o $(BIN)/UGKS2D src/UGKS2D.f90

#build manual
manual: 
	cd doc; latex  -shell-escape manual
	cd doc; bibtex manual
	cd doc; latex  -shell-escape manual
	cd doc; latex  -shell-escape manual
	cd doc; dvips  manual
	cd doc; ps2pdf manual.ps

#clean
clean:
	rm -f bin/*
	cd doc; rm -f !(*.tex|*.bib|*.pdf|*.eps)
