import zipfile
import os
import sys, shutil
from contextlib import closing
from zipfile import ZipFile, ZIP_DEFLATED


def zipdir(basedir, archivename):
    assert os.path.isdir(basedir)
    with closing(ZipFile(archivename, "w", ZIP_DEFLATED)) as z:
        for root, dirs, files in os.walk(basedir):
            #NOTE: ignore empty directories
            for fn in files:
                absfn = os.path.join(root, fn)
                zfn = absfn[len(basedir)+len(os.sep):] #XXX: relative path
                z.write(absfn, zfn)
                
for sub in xrange (2,11):  
    path = './sub%02d/'% sub
    folder_list = os.walk(path).next()[1]

    #Start zipping the folders
    for each_folder in folder_list:
        foldername = os.getcwd() + path[1:]+ each_folder
        zipname = path[2:-1]+"_"+each_folder
        print "Compressing .. " + foldername
#        shutil.make_archive(zipname, 'zip',foldername )
        zipdir(foldername,zipname)

