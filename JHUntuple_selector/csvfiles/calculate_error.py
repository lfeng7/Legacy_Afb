# Calculate systematics of final yields

import math

with open('yields_scaled.csv','r') as f0:
    nominal = f0.readlines()
with open('yields_scaledup.csv','r') as f1:
    up = f1.readlines()
with open('yields_scaleddown.csv','r') as f2:
    down = f2.readlines()
fout = open('yields_with_error.csv','w')
# sample_type nevts error

n_singletop,n_wjets,n_zjets,n_ttbar,n_mc,n_data = 0,0,0,0,0,0
u_singletop,u_wjets,u_zjets,u_ttbar,u_mc        = 0,0,0,0,0
d_singletop,d_wjets,d_zjets,d_ttbar,d_mc        = 0,0,0,0,0


for i in range(len(nominal)):	
    if nominal[i].split(' ')[0]=='singletop':
        n_singletop += float(nominal[i].split(' ')[1])
        u_singletop += float(up[i].split(' ')[1])
        d_singletop += float(down[i].split(' ')[1])
    if nominal[i].split(' ')[0]=='wjets':
        n_wjets += float(nominal[i].split(' ')[1])
        u_wjets += float(up[i].split(' ')[1])
        d_wjets += float(down[i].split(' ')[1])
    if nominal[i].split(' ')[0]=='zjets':
        n_zjets += float(nominal[i].split(' ')[1])
        u_zjets += float(up[i].split(' ')[1])
        d_zjets += float(down[i].split(' ')[1])
    if nominal[i].split(' ')[0]=='ttbar':
        n_ttbar += float(nominal[i].split(' ')[1])
        u_ttbar += float(up[i].split(' ')[1])
        d_ttbar += float(down[i].split(' ')[1])
    if nominal[i].split(' ')[0]!='data':
        n_mc += float(nominal[i].split(' ')[1])
        u_mc += float(up[i].split(' ')[1])
        d_mc += float(down[i].split(' ')[1])
    if nominal[i].split(' ')[0]=='data':
        n_data += float(nominal[i].split(' ')[1])

err_singletop = int(max(abs(n_singletop-u_singletop),abs(n_singletop-d_singletop)))
err_wjets = int(max(abs(n_wjets-u_wjets),abs(n_wjets-d_wjets)))
err_zjets = int(max(abs(n_zjets-u_zjets),abs(n_zjets-d_zjets)))
err_ttbar = int(max(abs(n_ttbar-u_ttbar),abs(n_ttbar-d_ttbar)))
err_mc = int(max(abs(n_mc-u_mc),abs(n_mc-d_mc)))   
err_data = int(math.sqrt(n_data))     

fout.write('singletop '+str(int(n_singletop))+' '+str(err_singletop)+' '+str(int(u_singletop))+' '+str(int(d_singletop))+'\n')
fout.write('wjets '+str(int(n_wjets))+' '+str(err_wjets)+' '+str(int(u_wjets))+' '+str(int(d_wjets))+'\n')
fout.write('zjets '+str(int(n_zjets))+' '+str(err_zjets)+' '+str(int(u_zjets))+' '+str(int(d_zjets))+'\n')
fout.write('ttbar '+str(int(n_ttbar))+' '+str(err_ttbar)+' '+str(int(u_ttbar))+' '+str(int(d_ttbar))+'\n')
fout.write('mc '+str(int(n_mc))+' '+str(err_mc)+' '+str(int(u_mc))+' '+str(int(d_mc))+'\n')
fout.write('data '+str(int(n_data))+' '+str(err_data))

fout.close()

