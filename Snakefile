
import os
from os import path

configfile: "config.yml"
workdir: config["workdir"]
