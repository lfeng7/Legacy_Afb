from utility import *

f_yields = open('yields.txt','w')

sample_yields = [['a',3,2,1],['b',3,2,1],['a',3,2,1]]

MC_yields = []

all_types = GetListChoices[ sample[0] for sample in sample_yields]

for itype in all_types :
    ilist = [ sample for sample in sample_yields if sample[0] == itype]
    MC_yields.append(SumColumn(ilist))

mc_total_yields = SumColumn(MC_yields)
mc_total_yields.pop(0)
MC_yields.append(mc_total_yields)
# Yields for data
data_yields = ['data',9,6,3]
# Write into yields output files
f_yields.write('sample,nocut, el selection, loose mu veto, dilep veto, jets selection, b-tagging, MET \n')
for row in MC_yields :
    for item in row :
        f_yields.write(str(item)+',')
    f_yields.write('\n')
for item in data_yields :
    f_yields.write(str(item)+',')

f_yields.close()

