import ROOT
import sys
import glob

argv = sys.argv[1:]

infile = argv.pop(0)
if len(argv)>0 :
    cut = argv.pop(0)
else:
    cut = ''

chain = ROOT.TChain("selected")
files= glob.glob(infile)
for item in files:
    print item
    chain.Add(item)

deno = ttree.GetEntries(cut)
if cut is not '':
    cut2 = '%s && %s'%(cut,'lep_iso<0.1')
    cut3 = '%s && %s'%(cut,'lep_iso>0.2&&lep_iso<1.2')
else:
    cut2 = 'lep_iso<0.1'
    cut3 = 'lep_iso>0.2&&lep_iso<1.2'
num_B = ttree.GetEntries(cut2)
num_A = ttree.GetEntries(cut3)
iso_eff = num_B*1.0/deno
rej_fact = num_B*1.0/num_A
print '\ncuts: %s ,%i evts passed'%(cut,deno) 
print 'isolation efficiency of iso<0.1 = %.4f'%iso_eff
print 'Rejection factor = %.4f'%rej_fact

tfile.Close()
