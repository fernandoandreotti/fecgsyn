###!/usr/bin/env python
## This script extract headers from MATLAB functions for the FECGSYN toolbox
import os

startdir = "/media/andreotti/Data/git/fecgsyn/"

print "Opening the file..."
docs = open(startdir+'docsout', 'w')

# loops through files
for root, dirs, files in os.walk(startdir):
    for file in files:
        if file.endswith(".m"):
             print("Obtaining header for "+os.path.join(root, file)+"...")
             # here we have the directory for each
             docs.write(file[:-2]+'\n')
             print file
docs.close() # you can omit in most cases as the destructor will call it


