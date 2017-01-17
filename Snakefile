
import os
from os import path

configfile: "config.yml"
workdir: path.join(config["workdir_top"], config["pipeline"])

WORKDIR = path.join(config["workdir_top"], config["pipeline"])
RESDIR =  config["resdir"]
SNAKEDIR = path.dirname(workflow.snakefile)
PY2_EXEC = "python2 {}/scripts".format(SNAKEDIR)

include: "snakelib/utils.snake"
