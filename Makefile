SHELL := /bin/bash
.DEFAULT_GOAL := help

PYTHON_EXEC = python
ifeq (, $(shell which python))
	PYTHON_EXEC = python3
endif

help: ## Show this help
	@echo Dependencies: $(PYTHON_EXEC)
	@egrep -h '\s##\s' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'

install: ## Install project dependencies
	$(PYTHON_EXEC) -m pip install -r requirements.txt

data: ## Download datasets
	$(PYTHON_EXEC) download_data.py

check_occurrences: ## Calculate the number of occurrences of TEA015198 in a dataset
	$(PYTHON_EXEC) check_occurrences.py