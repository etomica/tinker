#!/usr/bin/python3

import os
import sys
import getopt

if not len(sys.argv) == 4:
  print("Usage: extract.py input.prm output.prm tag1,tag2,tag3...")
  sys.exit()

if not os.path.exists(sys.argv[1]):
  print("Usage: extract.py input.prm output.prm tag1,tag2,tag3...")
  print("   {0} does not exist".format(sys.argv[1]))
  sys.exit()

intagss = sys.argv[3].split(",")
intags = [int(intagss[i]) for i in range(len(intagss))]

inpf = open(sys.argv[1], 'r')
outf = open(sys.argv[2], 'w')
tags = {}
classes = {}

classkw = {"vdw": [1],
           "bond": [1,2],
           "angle": [1,2,3],
           "strbnd": [1,2,3],
           "torsion": [1,2,3,4]}
tagkw = {"multipole": [1,2,3],
         "polarize": [1]}

reachedAtom = False

for line in inpf:
  if line[0] == ' ':
    continue
  bits = line.split()
  if len(bits) == 0:
    continue
  bit1 = bits[0].lower()
  if bit1 == "atom":
    reachedAtom = True
    thistag = int(bits[1])
    if not thistag in intags:
      continue
    newtag = len(tags)+1
    bits[1] = str(newtag)
    tags[thistag] = newtag
    thisclass = int(bits[2])
    if not thisclass in classes:
      newclass = len(classes)+1
      classes[thisclass] = newclass
    newclass = classes[thisclass]
    bits[2] = str(newclass)
    outf.write(" ".join(bits)+"\n")
  elif not reachedAtom:
    outf.write(line)
  elif bit1.lower() in tagkw:
    kw = bit1.lower()
    kw = tagkw[kw]
    skipme = False
    for i in range(len(kw)):
      itag = int(bits[kw[i]])
      if not itag in tags:
        skipme = True
        break
      bits[kw[i]] = str(tags[itag])
    if not skipme:
      outf.write(" ".join(bits)+"\n")
  elif bit1.lower() in classkw:
    kw = bit1.lower()
    kw = classkw[kw]
    skipme = False
    for i in range(len(kw)):
      iclass = int(bits[kw[i]])
      if not iclass in classes:
        skipme = True
        break
      bits[kw[i]] = str(classes[iclass])
    if not skipme:
      outf.write(" ".join(bits)+"\n")

inpf.close()
outf.close()

  
