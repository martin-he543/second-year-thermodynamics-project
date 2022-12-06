import numpy as np
import readline as rd

file1 = open('simulation.txt', 'r')
Lines = file1.readlines()

for line in Lines:  print("{}".format(line.rstrip()))