#!/bin/env python
from sys import argv
from squid.jobs import run_nbs_cmd as rnc

print rnc("jsub %s" % ' '.join(argv[1:])).stdout.read()

