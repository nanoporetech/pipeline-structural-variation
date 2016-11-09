
import os
from os import path

configfile: "config.yml"
workdir: path.join(config["workdir_top"], config["pipeline"])

WORKDIR = path.join(config["workdir_top"], config["pipeline"])
RESDIR =  config["resdir"]
SNAKEDIR = path.dirname(workflow.snakefile)
PY2_EXEC = "cd {}/scripts; python2 ".format(SNAKEDIR)

include: "snakelib/utils.snake"
