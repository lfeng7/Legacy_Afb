from utility import *

f_ = loadroot('selected_files/ejets/all/DY1JetsToLL_M_selected.root')
ttree = f_.Get('selected')
n_pdgid0 = 0
n_jets = 0
for i in range(ttree.GetEntries()):
    ttree.GetEntry(i)
    for item in ttree.jets_flavor:
        if item == 0 : n_pdgid0 += 1
        n_jets+= 1
print 'num of entries',ttree.GetEntries()
print 'num of jets',n_jets
print 'num jets with pdgid 0',n_pdgid0
