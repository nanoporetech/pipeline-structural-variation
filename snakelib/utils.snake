import re

def generate_help(sfile):
    """Parse out target and help message from file."""
    handler = open(sfile, "r")
    for line in handler:
        match = re.match(r'^rule\s+([a-zA-Z_-]+):.*?## (.*)$$', line)
        if match:
            target, help = match.groups()
            print("%-20s %s" % (target, help))

rule help: ## print list of all targets with help
    input:
        workflow.included
    run:
        print("--------------------------------------------")
        [generate_help(sfile) for sfile in input]
        print("--------------------------------------------")

rule clean_workdir: ## delete working directory. WARNING: all data will be lost!
    input:
        wdir = {WORKDIR}
    shell: "rm -r {input.wdir}"

rule clean_resdir: ## delete results directory. WARNING: all data will be lost!
    input:
        res = RESDIR
    shell: "rm -fr {input.res}/*"

rule info: ## print pipeline information
    params:
        name = config["pipeline"],
        wdir = WORKDIR,
        repo = config["repo"],
        res  = config["resdir"],
    run:
        print("Pipeline name: ", params.name)
        print("Pipeline working directory: ", params.wdir)
        print("Pipeline results directory: ", params.res)
        print("Pipeline repository: ", params.repo)
    
