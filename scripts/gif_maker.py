#!/usr/bin/env python

import sys, os, shutil
import argparse
from argparse import RawTextHelpFormatter
import subprocess

parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description='This script is designed to take a .mov formatted screen recording or video and transform it into a gif. Requires installation of ffmpeg and gifsicle which can be installed via Anaconda. \n\n\
::CITE:: \n\
https://github.com/reptilerhett/RandomScripts\n\n')

###############################################

parser.add_argument("-i","--input",
					type=str,
                    default="input.mov",
                    help="Input in .mov format")
parser.add_argument("-o","--output",
					type=str,
                    default="out.gif",
					help="Output in .gif format")
parser.add_argument("-f","--fps",
					type=str,
                    default="10",
					help="Frames per second. Default = 10")
parser.add_argument("-s","--size",
					type=str,
                    default="720x480",
					help="Max width and height. Default = 720x480")
args=parser.parse_args()

###############################################

input = args.input
output = args.output
size = args.size
fps = args.fps

###############################################

command = 'ffmpeg -i ' + input + ' -s ' + size + ' -pix_fmt rgb24 -r ' + fps + ' -f gif - | gifsicle --optimize=3 --delay 3 > ' + output
subprocess.call(command,shell=True)
