#!/usr/bin/env python

import sys
import math

def round_coordinate_down(coord, r):
  coord = int(coord)
  coord = math.floor(coord/r) *r
  return(str(coord))

def round_coordinate_up(coord, r):
  coord = int(coord)
  coord = math.ceil(coord/r) *r
  return(str(coord))
  
if len(sys.argv) == 1:
  window = 5000
else:
  window = int(sys.argv[1])

for line in sys.stdin:
  line = line.strip().split("\t")
  line[1] = round_coordinate_down(line[1],window)
  line[2] = round_coordinate_up(line[2],window)
  line[4] = round_coordinate_down(line[4],window)
  line[5] = round_coordinate_up(line[5],window)
  print("{}".format("\t".join(line)))
