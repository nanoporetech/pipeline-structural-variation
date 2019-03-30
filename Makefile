# Name fo the project (default: infered from directory name)
PROJECT_NAME?=$(shell pwd | rev | cut -d '/' -f 1 | rev)
# Version number (using: git describe --all)
PROJECT_VERSION?=$(shell grep -oE "[0-9]+.[0-9]+.[0-9]+" lib/pipeline_structural_variation/__init__.py)


# Build docker image and tags it with required name
build:
	@echo "Building docker file"
	docker build -t $(PROJECT_NAME):$(PROJECT_VERSION) -f Dockerfile .
	#--no-cache
.PHONY: build

# Run pipeline with test data
test-docker:
	@echo "Testing pipeline"
	docker run -ti -w `pwd` -v `pwd`:`pwd` $(PROJECT_NAME):$(PROJECT_VERSION) snakemake -p --snakefile Snakefile all
.PHONY: test-docker

clean:
	@echo "Cleaning working dir"
	snakemake -p clean_workdir
.PHONY: clean

test:
	@echo "Testing pipeline"
	snakemake --snakefile Snakefile -p all
.PHONY: test

docker-env:
	@echo "Starting dev env in docker"
	docker run -ti -w `pwd` -v `pwd`:`pwd` $(PROJECT_NAME):$(PROJECT_VERSION) bash
.PHONY: docker-env

conda-install:
	@echo "Installing conda env"
	conda env create -n ${PROJECT_NAME} -f env.yml
.PHONY: conda-install

conda-uninstall:
	@echo "Removing conda env"
	conda remove -n ${PROJECT_NAME} --all -y
.PHONY: conda-uninstall

version:
	@echo ${PROJECT_VERSION}

eval:
	snakemake -p eval --config bam_folder=/Users/prescheneder/Analysis/sv-benchmark-data/chr12_data/GM24385_minimap2_sv_q7_chr12.bam target=/Users/prescheneder/Analysis/sv-benchmark-data/chr12_data/target.bed
.PHONY: eval