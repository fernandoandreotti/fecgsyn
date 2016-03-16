###!/usr/bin/env python
# This script extract headers from MATLAB functions for the FECGSYN toolbox
# Escape strings are '% --' or '% fecgsyn toolbox, version'
#
import os, string

limdepth = 1 # only using 3 first levels on documentation
path = os.path.normpath(".")
docs = open(path+'/docsout', 'w') # open output file
path = os.path.normpath("./subfunctions/") # scan subfunctions folder only

# loops through files
for root, dirs, files in os.walk(path):
    depth = root[len(path) + len(os.path.sep):].count(os.path.sep)
    if "/libs" in root: # must escape libs directories in all levels
         continue 
    if depth > limdepth: # limiting depth of search
        continue
    if depth == 0:
         docs.write('<section class="bs-docs-section"> \n')
         docs.write('<h1 id="'+root[len(path)+1:]+'">'+root[len(path)+1:]+'</h1><br> \n')
    else:
         docs.write('<h2 id="'+root[len(path)+1:]+'">'+root[len(path)+1:]+'</h2><br> \n')
    print 'Sub-directory '+root
    for file in files:
        if file.endswith(".m"):
             print("Obtaining header for "+os.path.join(root, file)+"...")
             # here we have the directory for each
             docs.write('<h3 id="'+file[:-2]+'">'+file[:-2]+'</h3><br> \n')
             # now have to open each file, exclude first line and copy every line until 
             f = open(os.path.join(root, file), 'r') # open output file            
             while  f.readline(): # skip first line
                     line = f.readline()
                     if (line.rstrip()=="%"):
                         continue # skip blank lines
                     if ("% --" in line) |( "fecgsyn toolbox, version" in line) :  # where to stop reading file
                        line = "quit"
                        break
                     if ("input:" in line.lower())|("inputs:" in line.lower()):
                        line = "% <b>Input:</b> \n"
                     elif ("output:" in line.lower())|("outputs:" in line.lower()):
                        line = "% <b>Output:</b> \n"
                     elif "examples" in line.lower():
                        line = "% <b>Examples:</b> <br> \n"
                     elif "% function" in line.lower():
                        line = "% <b>Call: </b> <code>" + line[2:] + "</code> <br> \n"
                     elif "see also:" in line.lower():
                        docs.write("<b>See also:</b> <br> \n")
                        while f.readline():                            
                            line = f.readline()
                            if ("% --" in line) |(line.rstrip()=="%")| ( "fecgsyn toolbox, version" in line) :
                                line = "quit"
                                break
                            line = '<a href="{{site.github.url}}/pages/documentation.html#'+line[2:]+'"><code>'+line[2:]+'</code></a></code> \n'
                            docs.write(line+'<br>')                          
                     if line == "quit":
                        break
                     else:
                        docs.write(line[2:]+'<br> \n')
             f.close()
    docs.write('</section>\n\n')
docs.close() # you can omit in most cases as the destructor will call it


