fout = open('ana.listOfJobs','w')
num_files = 5
filesperjob =  10
startfile = 0
txtfile = 'SingleEl_Run2012C_part2.txt'
sampletype = 'SingleEl_Run2012C_part2'
mcordata = 'data'
while startfile < num_files :
	toprint = 'python ./tardir/selection.py --txtfiles tardir/inputfiles/'+txtfile+' --makeplots no --mcordata '+mcordata
	toprint += ' --maxevts -1 --type '+sampletype+' --grid yes --maxfiles '+str(filesperjob)+' --startfile '+str(startfile)+'\n'
	fout.write(toprint)
	startfile += filesperjob

fout.close()
