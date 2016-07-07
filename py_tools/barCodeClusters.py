##Line Graph showing the coverage of barcodes across the genome
##Approximates bargraph with large data-set, loads fast
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cbook as cbook



f = open('/Users/awr/Desktop/barCodeFreq.csv', "r")
lines = f.read().split("\n") # "\r\n" if needed

x = []
y = []

for line in lines:
    if line != "": #skips headers
        cols = line.split(",")
        x.append(cols[0])
        y.append(cols[1])
        