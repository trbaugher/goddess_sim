#!/usr/bin/python
import sys

def main(argv):
  minX = 999e999
  minY = 999e999
  minZ = 999e999
  maxX = 0
  maxY = 0
  maxZ = 0
  f = open(argv[1], 'r')
  print argv[1]
  for lines in f:
    line = lines.split()
    if (line[0] == "vertex"):
      if (float(line[1]) > maxX):
        maxX = float(line[1])
      if (float(line[2]) > maxY):
        maxY = float(line[2])
      if (float(line[3]) > maxZ):
        maxZ = float(line[3])
      if (float(line[1]) < minX):
        minX = float(line[1])
      if (float(line[2]) < minY):
        minY = float(line[2])
      if (float(line[3]) < minZ):
        minZ = float(line[3])



  print "stlMinX=",minX,";stlMinY=",minY,";stlMinZ=",minZ,";"
  print "stlMaxX=",maxX,";stlMaxY=",maxY,";stlMaxZ=",maxZ,";"
  print "offset = G4ThreeVector(-((stlMaxX-stlMinX)/2.0+stlMinX), -((stlMaxY-stlMinY)/2.0+stlMinY), -((stlMaxZ-stlMinZ)/2.0+stlMinZ));"

if __name__ == "__main__":
  main(sys.argv)
