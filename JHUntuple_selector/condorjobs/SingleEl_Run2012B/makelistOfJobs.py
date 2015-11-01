fout = open('ana.listOfJobs','w')
num_files = 636 
filesperjob =  10
startfile = 0
txtfile = 'SingleEl_Run2012B.txt'
sampletype = 'SingleEl_Run2012B'
mcordata = 'data'
while startfile < num_files :
	toprint = 'python ./tardir/selection.py --txtfiles tardir/inputfiles/'+txtfile+' --makeplots no --mcordata '+mcordata
	toprint += ' --applyHLT yes --maxevts -1 --type '+sampletype+' --grid yes --maxfiles '+str(filesperjob)+' --startfile '+str(startfile)+'\n'
	fout.write(toprint)
	startfile += filesperjob

fout.close()
