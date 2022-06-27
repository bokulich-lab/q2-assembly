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
	$(PYTHON) setup.py install

dev: all
	pip install pre-commit
	pip install -e .
	pre-commit install

clean: distclean

distclean: ;
