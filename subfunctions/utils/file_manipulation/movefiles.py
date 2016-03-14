###!/usr/bin/env python
# This script is used to divide simulations into smaller subfolders for further uploading to the Physionet website
import os, string, shutil, fnmatch
 
 # subfunction to generate path names
def gen_find(pattern,filelist):
    for name in fnmatch.filter(filelist,pattern):
        yield os.path.join(path,name)



limdepth = 1# only using 3 first levels on documentation
path = os.path.normpath(".")

# loops through files
for root, subdir, files in os.walk(path):
    depth = root[len(path) + len(os.path.sep):].count(os.path.sep)    
    if depth > limdepth: # limiting depth of search
        continue    
    if len(files)<2:
        continue
    for snr in xrange (0, 13, 3):        
        pattern = "_snr%02ddB_" % snr
        filesToMove = gen_find("*"+pattern+"*",files)
        print('The current folder is ' + root)
        print pattern[1:-1]
        subd = root+"/" + pattern[1:-1]
        if not os.path.exists(subd):
            os.makedirs(subd)
        for name in filesToMove:
            print "Moving .. " + name[2:]
            shutil.move(root+name[1:], subd+name[1:])
        #shutil.move(filename, root)
