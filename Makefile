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
	pip install git+https://github.com/ablab/quast.git@bc4af762a7f53176d66bd5e6c5b7d28376d28e11
	$(PYTHON) setup.py install

dev: all
	pip install -e .

clean: distclean

distclean: ;
