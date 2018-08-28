#!/bin/env python

from squid.jobs import run_nbs_cmd as rnc

print rnc("jlist").stdout.read()

