fout = open('ana.listOfJobs','w')
num_files = 5 
filesperjob =  10
startfile = 0
txtfile = 'Tbar_tW.txt'
sampletype = 'Tbar_tW'
mcordata = 'mc'
while startfile < num_files :
	toprint = 'python ./tardir/selection.py --txtfiles tardir/inputfiles/'+txtfile+' --makeplots no --mcordata '+mcordata
	toprint += ' --mctype singletopbar  --maxevts -1 --type '+sampletype+' --grid yes --maxfiles '+str(filesperjob)+' --startfile '+str(startfile)+'\n'
	fout.write(toprint)
	startfile += filesperjob

fout.close()
