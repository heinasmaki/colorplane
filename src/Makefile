FC = gfortran
FFLAGS ?= -O3
PYTHON ?= python3

.PHONY: all colorplane test
all: colorplane

colorplane:
	$(PYTHON) ../scripts/build.py --compiler "$(FC)" --flags="$(FFLAGS)" --output "$(CURDIR)/colorplane"

test:
	$(PYTHON) ../tests/test_reliability.py --compiler "$(FC)"
