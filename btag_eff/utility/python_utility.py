# Adding two lists by adding each element in the lists
# Unless that element is string. In that case, take the element from first list
import sys
import os
import math
import glob
import numpy
from array import array

def SumColumn( M ) :
    newlist = []
    for i in range(len(M[0])) :
        icol = [ row[i] for row in M ]
        if type(icol[0]) is str : newlist.append(icol[0])
        else : newlist.append(sum(icol)) 
    return(newlist) 

def SumRow( M ) :
    return [ sum(row) for row in M]

# Find all non-identical possibilities in the list
def GetListChoices(l):
    tmp = []
    for item in l :
        if item not in tmp : tmp.append(item)
    return tmp

def ListCompare( a,b) :
    diff_a = [ i for i in a if i not in b ]
    diff_b = [ i for i in b if i not in a ] 
    if len(diff_a) >= len(diff_b) : return diff_a
    else : return diff_b
 
def GetSomeFiles(allfiles,startfile,filesperjob):
    # Only keep certain number of input files for fexibility
    if filesperjob <= 0 : somefiles = allfiles   # maxfiles<= 0 indicates run all files
    else :
        headfile = startfile
        endfile = headfile+filesperjob
        endfile = min(endfile,len(allfiles)) # end file cannot be larger than total number of files
        headfile = min(headfile,endfile)   # headfile must be smaller or equal to endfile
        print 'Will process file',headfile,'to',endfile
        somefiles = [allfiles[i] for i in range(headfile,endfile)]  
    return somefiles

# Make directory if not exist already
def MakeDirectory(dirname):
    prefix = dirname
    # check if plotting dir is made. If not , make it now
    if not os.path.exists(prefix):
        os.mkdir(prefix)
        print 'Making '+prefix
    else :
        print 'dir',dirname,'exists already'
    os.system('cp utility/index.php '+dirname)

def Upload(dirname,webdirname=''):
    prefix = dirname
    webdir = 'pha:/home/lfeng/public_html/research/Dump/'+webdirname
    # check if plotting dir is made. If not , make it now
    if not os.path.exists(prefix):
        print 'Dir ',dirname,'does not exist. Cannot upload to webpage!'
    else :
        os.system('scp -r '+dirname+' '+webdir)
        print 'Uploading',dirname,'to',webdir    

