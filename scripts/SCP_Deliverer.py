#!/usr/bin/env python

import argparse
import sys, os, shutil
import subprocess
import glob

# Command line options
parser = argparse.ArgumentParser(description='Some dumb f*cked script I am writing')
parser.add_argument("-l","--localpath",
					type=str,
					help="local file or directory for transfer")
parser.add_argument("-c","--clusterpath",
                    type=str,
                    help="remote/cluster file or directory for transfer")
parser.add_argument("-s","--server",
					type=str,
					default="@xfer01-ext.palmetto.clemson.edu",
					help="Server name (default: @xfer01-ext.palmetto.clemson.edu)")
parser.add_argument("-u","--username",
                    type=str,
                    default="rrautsa",
                    help="Username for cluster computer (default: rrautsa)")
parser.add_argument("-d","--direction",
                    type=str,
                    help="Direction of Transfer (e.g. c2l= cluster to local; l2c = local to cluster)")
args=parser.parse_args()

if args.direction == 'c2l':
	command = "scp -r " + args.username + args.server + ":" + args.clusterpath + " " + args.localpath
	subprocess.call(command,shell=True)
elif args.direction == 'l2c':
	command = "scp -r " + args.localpath + " " + args.username + args.server + ":" + args.clusterpath
	subprocess.call(command,shell=True)
else:
	exit()	
	
	
