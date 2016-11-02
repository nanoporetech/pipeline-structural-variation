WORKDIR = path.join(config["workdir_top"], config["pipeline"])
RESDIR =  config["resdir"]

import os
from os import path

configfile: "config.yml"
workdir: path.join(config["workdir_top"], config["pipeline"])

include: "snakelib/utils.snake"
