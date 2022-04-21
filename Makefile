.PHONY: all lint test test-cov install dev clean distclean

PYTHON ?= python

all: ;

lint:
	q2lint
	flake8

test: all
	py.test

test-cov: all
	coverage run -m pytest
	coverage xml

install: all
	pip install git+https://github.com/ablab/quast.git@de7e1f4891a3487f3c0df6ae27cbfba38734d686
	$(PYTHON) setup.py install

dev: all
	pip install git+https://github.com/ablab/quast.git@de7e1f4891a3487f3c0df6ae27cbfba38734d686 coverage
	pip install -e .

clean: distclean

distclean: ;
