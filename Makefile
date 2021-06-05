SHELL := /bin/bash
.DEFAULT_GOAL := help

.PHONY: install data analyze simulate docs

PYTHON_EXEC = python
ifeq (, $(shell which python))
	PYTHON_EXEC = python3
endif

help: ## Show this help
	@echo Dependencies: $(PYTHON_EXEC) [pdoc3]
	@egrep -h '\s##\s' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'

docs: ## Generates the documentation files (requires pdoc3)
	pdoc3 --html --output-dir docs --force .
	cp docs/theanine-synthetase/* docs/ -R
	rm -rf docs/theanine-synthetase

install: ## Install project dependencies
	$(PYTHON_EXEC) -m pip install -r requirements.txt

data: ## Download datasets
	$(PYTHON_EXEC) download_data.py

analyze: ## Run the analysis
	$(PYTHON_EXEC) analysis.py

simulate: ## Run the Monte-Carlo simulation
	$(PYTHON_EXEC) simulation_orchestrator.py