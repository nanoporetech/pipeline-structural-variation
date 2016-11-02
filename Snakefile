
import os
from os import path

configfile: "config.yml"
workdir: config["workdir"]

include: "snakelib/utils.snake"
