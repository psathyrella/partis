.PHONY : clean-pyc clean-cython clean-python clean inplace build test

PYTHON ?= python

CYTHON_COMPILED  = vdjalign/bwa.c

CYTHON_FILES = $(wildcard vdjalign/*.pyx)

LIBRARIES = $(CYTHON_FILES:pyx=so)

inplace: $(EXTS)
	$(PYTHON) setup.py build_ext --inplace

$(EXTS): $(CYTHON_FILES)

build:
	$(PYTHON) setup.py build

clean: clean-python clean-pyc clean-cython

clean-pyc:
	find vdjalign -name '*.pyc' -exec rm {} \;

clean-cython:
	rm -f $(CYTHON_COMPILED) $(LIBRARIES)
	rm -rf build dist

clean-python:
	$(PYTHON) setup.py clean
	rm -rf build dist

test:
	$(PYTHON) setup.py test
