.ONESHELL:

SHELL := /bin/bash

CONDA_ENV_NAME=T_A
# Note that the extra activate is needed to ensure that the activate floats env to the front of PATH
CONDA_ACTIVATE=source $$(micromamba info --base)/etc/profile.d/micromamba.sh ; micromamba activate ; micromamba activate

build-conda-from-req: ## Build the conda environment
	micromamba create -n $(CONDA_ENV_NAME) --copy -y python=$(PY_VERSION)
	$(CONDA_ACTIVATE) $(CONDA_ENV_NAME)
	python -s -m pip install -r requirements.txt

build-conda-from-env:
	micromamba env create -n $(CONDA_ENV_NAME) -f environment.yml

install-attnpacker:
	micromamba activate attnpacker
	git clone https://github.com/rostro36/AttnPacker.git
	pip install -r AttnPacker/requirements.txt

install-environment:
	micromamba env create -n $(CONDA_ENV_NAME) -f T_A.yml

install-faspr:
	$(CONDA_ACTIVATE) $(CONDA_ENV_NAME)
	git clone https://github.com/tommyhuangthu/FASPR.git
	cd FASPR
	g++ -O3 --fast-math -o FASPR src/*.cpp

