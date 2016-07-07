###=====================================
### Creates visualations of Barcode Data
### By Will Richards
###=====================================

import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cbook as cbook
import pylab as pl
import collections
from collections import Counter
matplotlib.use('Agg')

f = open('/uufs/chpc.utah.edu/common/home/u0401321/TenX/bin/tools/barCodeFreq.csv', "r")
lines = f.read().split("\n") # "\r\n" if needed

x = []
y = []

for line in lines:
    if line != "": # add other needed checks to skip titles
        cols = line.split(",")
        x.append(cols[0])
        y.append(cols[1])



###==============================================================
### Creates a bargraph showing the count of each unique barcode
###===============================================================
import pylab as pl

x_axis = range(0,len(x))
xTicks = x
pl.xticks(range(len(x)), xTicks, rotation=45) #writes strings with 45 degree angle
pl.plot(x_axis,y,'*')
pl.show()



###========================================================================
###Creates A Histogram of the count of the frequency of each unique barcode
###========================================================================

y=map(int,y)
labels, values = zip(*Counter(y).items())
indexes = np.arange(len(labels))
width = 1

plt.title("Barcode Frequency Distribution")
plt.xlabel('Barcode Frequency')
plt.ylabel('Count of Frequencies')
plt.xticks(indexes + width * 0.5  , labels)
plt.bar(indexes, values, width )
plt.show()
