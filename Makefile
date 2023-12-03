.ONESHELL:

SHELL := /bin/bash

CONDA_ENV_NAME=T_A
# Note that the extra activate is needed to ensure that the activate floats env to the front of PATH
CONDA_ACTIVATE=source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate

install-attnpacker:
	$(CONDA_ACTIVATE) $(CONDA_ENV_NAME)
	git clone https://github.com/MattMcPartlon/AttnPacker.git
	pip install -r AttnPacker/requirements.txt

install-faspr:
	$(CONDA_ACTIVATE) $(CONDA_ENV_NAME)
	git clone https://github.com/tommyhuangthu/FASPR.git
	cd FASPR
	g++ -O3 --fast-math -o FASPR src/*.cpp

build-conda-from-req: ## Build the conda environment
	conda create -n $(CONDA_ENV_NAME) --copy -y python=$(PY_VERSION)
	$(CONDA_ACTIVATE) $(CONDA_ENV_NAME)
	python -s -m pip install -r requirements.txt

build-conda-from-env:
	conda env create -n $(CONDA_ENV_NAME) -f environment.yml

new-env:
	$(CONDA_ACTIVATE) $(CONDA_ENV_NAME)
	conda env export > environment.yml

new-req:
	$(CONDA_ACTIVATE) $(CONDA_ENV_NAME)
	pip freeze > requirements.txt

format:
	$(CONDA_ACTIVATE) $(CONDA_ENV_NAME)
	@echo -n "==> Checking that imports are properly sorted with isort..."
	@echo -n ""
	@isort .
	@echo -n "==> Checking that code is autoformatted with black..."
	@echo -n ""
	@black .
