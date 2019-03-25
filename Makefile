# Name fo the project (default: infered from directory name)
PROJECT_NAME?=$(shell pwd | rev | cut -d '/' -f 1 | rev)
# Version number (using: git describe --all)
PROJECT_VERSION?=$(shell git describe --all | rev | cut -f 1 -d '/' | rev | sed 's/HEAD/master/')
# ImageId used in CWL step definition files
PROJECT_IMAGE?=$(PROJECT_NAME)

# Build docker image and tags it with required name
build:
	@echo "Building docker file"
	docker build -t $(PROJECT_IMAGE):$(PROJECT_VERSION) -f Dockerfile .
	#--no-cache
.PHONY: build

# Run pipeline with test data
test-docker:
	@echo "Testing pipeline"
	docker run -ti -w `pwd` -v `pwd`:`pwd` pipeline-structural-variation:master snakemake -p --snakefile Snakemake.smk all
.PHONY: test-docker

clean:
	@echo "Cleaning working dir"
	snakemake --snakefile Snakemake.smk -p clean_workdir
.PHONY: clean

test:
	@echo "Testing pipeline"
	snakemake --snakefile Snakemake.smk -p all
.PHONY: test

docker-env:
	@echo "Starting dev env in docker"
	docker run -ti -w `pwd` -v `pwd`:`pwd` pipeline-structural-variation:master bash
.PHONY: docker-env

conda-install:
	@echo "Installing conda env"
	conda env create -n ${PROJECT_NAME} -f env.yml
.PHONY: conda-install

conda-uninstall:
	@echo "Removing conda env"
	conda remove -n ${PROJECT_NAME} --all -y
.PHONY: conda-uninstall