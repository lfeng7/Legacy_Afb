fout = open('ana.listOfJobs','w')
num_files = 282 
filesperjob =  10
startfile = 0
txtfile = 'TT_8TeV-mcatnlo.txt'
sampletype = 'TT_8TeV-mcatnlo'
mcordata = 'mc'
while startfile < num_files :
	toprint = 'python ./tardir/selection.py --txtfiles tardir/inputfiles/'+txtfile+' --makeplots no --mcordata '+mcordata
	toprint += ' --signal yes --maxevts -1 --type '+sampletype+' --grid yes --maxfiles '+str(filesperjob)+' --startfile '+str(startfile)+'\n'
	fout.write(toprint)
	startfile += filesperjob

fout.close()
