.PHONY: all lint test test-cov install dev clean distclean

PYTHON ?= python

all: ;

lint:
	q2lint
	flake8

test: all
	py.test

test-cov: all
	py.test --cov=q2_assembly

install: all
	pip install git+https://github.com/ablab/quast.git@de7e1f4891a3487f3c0df6ae27cbfba38734d686
	$(PYTHON) setup.py install

dev: all
	pip install -e .

clean: distclean

distclean: ;
