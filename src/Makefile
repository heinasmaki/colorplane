FC = gfortran
FFLAGS ?= -O3
PYTHON ?= python3
OPENMP ?= 0
OMPFLAG = $(if $(filter 1,$(OPENMP)),--openmp,)

.PHONY: all colorplane test
all: colorplane

colorplane:
	$(PYTHON) ../scripts/build.py --compiler "$(FC)" --flags="$(FFLAGS)" $(OMPFLAG) --output "$(CURDIR)/colorplane"

test:
	$(PYTHON) ../tests/test_reliability.py --compiler "$(FC)"
	$(PYTHON) ../tests/test_parallel.py --compiler "$(FC)"
