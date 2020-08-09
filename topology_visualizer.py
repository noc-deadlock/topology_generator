import sys
import subprocess
import numpy as np
# do the file-handling stuff first
input_ = sys.argv[1]
with open(input_, "r") as f:
    data_ = f.readline()
    data = data_.split(" ") # this creates a 'list' of original string contained in variable 'data_'
    rows_ = int(data[0])
    cols_ = int(data[1])

    next(f) #skip the line not necessary
    # now create a list and read all the file
    # content into that list, later refer to that list
    # to make connections in Mesh-figure
    # make an empty list

    topology = [[0 for x in range(cols_*rows_)] for y in range(rows_*cols_)]
    x, y = 0, 0
    for line in f:
        for word in line.split():
            topology[x][y] = int(word)
            y = y + 1
        x = x + 1 # increment row
        if((x == rows_ * cols_) and
                (y == rows_ * cols_)):
            break
        y = 0 # resent 'col' index on each row completion

# rows = [0, 1, 2, 3, 4, 5, 6, 7]
# rows = np.linspace(0, 8, 8, endpoint=False, dtype=int)
rows = range(0, rows_)
# cols = [0, 0, 0, 0, 0, 0, 0, 0]
cols = [0] * rows_
subprocess.call('tabs 4', shell=True) # set the spacing correctly
for r in range(len(rows)-1, -1, -1):
    for c in range(len(cols)):
        print("{0:d}\t".format(rows[r]*rows_+c)),
        if c < (len(cols) - 1):
            src = rows[r]*rows_+c   # router-id
            dst = src + 1       # router-id
            if topology[src][dst] == 1:
                assert (topology[dst][src] == 1)
                print("---"), # adding row-wise connections
            else:
                print("\t"),
    # insert a blank line after every row
    print(" ")
    for c in range(len(cols)):
        if (r > 0):
            src = rows[r] * rows_ + c  # router-id
            dst = src - rows_
            if topology[src][dst] == 1:
                assert(topology[dst][src] == 1)
                print("|\t")+str('    '),  # adding column-wise connections
            else:
                print("\t")+str("    "),
    print(" ")