"""
count total weight for mc@NLO file, or , any weighted MC
"""

from Legacy_Afb.Tools.fwlite_boilerplate import *
from Legacy_Afb.Tools.root_utility import *
import sys
import glob

argv = sys.argv[1:]

files = argv.pop(0)
files = glob.glob(files)
nfiles = len(files)
print 'Getting %i files'%nfiles

if argv!=[]:
    maxevt = argv.pop(0)
else:
    maxevt = 10000

outfile = open('total_weight.txt','w')
towrite = ''
# MC@NLO only
GenEventHandle = Handle("GenEventInfoProduct"); 
GenEventLabel  = ("generator","")

total_w = 0

events = Events(files)
evt_size = events.size() 
print 'Getting',evt_size,'events'

for i,evt in enumerate(events):
   if i==maxevt : 
       print 'Reach %i evt will break'%i
       break
   if i%10000 == 1: print 'Loop over %i evt'%i

   evt.getByLabel(GenEventLabel,GenEventHandle)
   if GenEventHandle.isValid():
       GenEvent = GenEventHandle.product()
       weight_gen = GenEvent.weight()
       total_w += weight_gen 

print '\n'
towrite += 'Total num of files = %i\n'%nfiles
towrite += 'Total num of evts = %i\n'%evt_size
if maxevt > 0 :
    evt_size = min(evt_size,maxevt)
towrite += 'Total weight = %.2f\n'%total_w
towrite += 'Average weight = %.2f'%(total_w/evt_size)

outfile.write(towrite)
outfile.close()

print towrite
