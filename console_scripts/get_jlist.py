#!/usr/bin/env python

'''
This script is used by jAutoTab for jdel autocomplete
'''

from squid.jobs import get_all_jobs

print " ".join(get_all_jobs(detail=0))
