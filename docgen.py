###!/usr/bin/env python
# This script extract headers from MATLAB functions for the FECGSYN toolbox
import os, string

limdepth = 1 # only using 3 first levels on documentation
path = os.path.normpath(".")

print "Opening the file..."
docs = open(path+'docsout2', 'w')

# loops through files
for root, dirs, files in os.walk(path):
    depth = root[len(path) + len(os.path.sep):].count(os.path.sep)
    if depth > limdepth:
        continue
    docs.write('<h1 id="'+root[len(path):]+'">'+root[len(path):]+'</h1>\n\n')
    print 'Sub-directory '+root
    for file in files:
        if file.endswith(".m"):
             print("Obtaining header for "+os.path.join(root, file)+"...")
             # here we have the directory for each
             docs.write('<h3 id="'+file[:-2]+'">'+file[:-2]+'</h3>\n\n')
             
docs.close() # you can omit in most cases as the destructor will call it


