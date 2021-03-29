#!/usr/bin/env python
from sys import argv
input = argv[1]
output = input[23:27]
from pyraf import iraf
iraf.imcopy(input=input, output=output, mode='hl', verbose='no')