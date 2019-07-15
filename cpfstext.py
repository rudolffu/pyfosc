#!/usr/bin/env python
from sys import argv
input = argv[1]
output = input[6:15]
from pyraf import iraf
iraf.imcopy(input=input+'[1]', output=output, mode='hl', verbose='no')
