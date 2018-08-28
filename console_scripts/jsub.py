#!/bin/env python
from sys import argv
from squid.jobs import run_nbs_cmd as rnc

cmd = "/opt/voyager/nbs/bin/jsub"
if len(argv) > 1:
    cmd += " " + " ".join(argv[1:])

print rnc(cmd).stdout.read()

