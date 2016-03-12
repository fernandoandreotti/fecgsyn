###!/usr/bin/env python
# This script extract headers from MATLAB functions for the FECGSYN toolbox
import os, string

limdepth = 1 # only using 3 first levels on documentation
path = os.path.normpath(".")
docs = open(path+'/docsout', 'w') # open output file
path = os.path.normpath("./subfunctions/") # scan subfunctions folder only

# loops through files
for root, dirs, files in os.walk(path):
    depth = root[len(path) + len(os.path.sep):].count(os.path.sep)
    if depth > limdepth: # limiting depth of search
        continue
    if depth == 0:
         docs.write('<h1 id="'+root[len(path)+1:]+'">'+root[len(path)+1:]+'</h1>\n\n')
    else:
         docs.write('<h2 id="'+root[len(path)+1:]+'">'+root[len(path)+1:]+'</h2>\n\n')
    print 'Sub-directory '+root
    for file in files:
        if file.endswith(".m"):
             print("Obtaining header for "+os.path.join(root, file)+"...")
             # here we have the directory for each
             docs.write('<h3 id="'+file[:-2]+'">'+file[:-2]+'</h3>\n\n')
             
docs.close() # you can omit in most cases as the destructor will call it


