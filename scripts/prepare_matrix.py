#!/usr/bin/env python                                                       
#                                                                           #
#                                                                           #
# by Mar Gonzalez-Ramirez (2022)                                            #
#############################################################################

import sys
import os

# Control for number of arguments
if len(sys.argv) != 3:
    print("Incorrect number of arguments.")
    print("Usage:")
    print("prepare_matrix.py <ChIPlevels> <expressionfile>")
    sys.exit(0)
    
# Read input
ChIPlevels = sys.argv[1]
expressionfile = sys.argv[2]

# check if ChIPlevels exists
if os.path.exists(ChIPlevels) == False:
        print(ChIPlevels+" does not exist")
        sys.exit(0)

# check if expressionfile exists
if os.path.exists(expressionfile) == False:
        print(expressionfile+" does not exist")
        sys.exit(0)
        
# save ChIP levels
f = open(ChIPlevels)
ChIPlevels_dict = {}
for l in f:
    ChIPlevelsfile,name = l.split()
    # check if ChIPlevelsfile exists
    if os.path.exists(ChIPlevelsfile) == False:
        print(ChIPlevelsfile+" does not exist")
        sys.exit(0)
    ChIPlevels_save = []
    targetgenes = []
    f2 = open(ChIPlevelsfile)
    for l2 in f2:
        col1,col2,col3,col4,col5,col6,col7 = l2.split()
        ChIPlevels_save.append(col5)
        targetgenes.append(col4)
        ChIPlevels_dict[name] = ChIPlevels_save
    f2.close()    
f.close()

# save expression
f = open(expressionfile)
genes = []
fpkms = []
for l in f:
    col1,col2 = l.split()
    genes.append(col1)
    fpkms.append(col2)
f.close()
expression = []
for tg in targetgenes:
    for i in range(len(genes)):
        if tg == genes[i]:
            expression.append(fpkms[i])

# prepare output
header = "expression"
for key in ChIPlevels_dict:
    header += "\t"
    header += key
print(header)
for i in range(len(expression)):
    line = expression[i]
    for key in ChIPlevels_dict:
        line += "\t"
        line += ChIPlevels_dict[key][i]
    print(line)
