#! usr/local/python
# look up root files that is already created, and find those jobs that contains not generated output root file, write into a new 
# ana.listOfJobs_resubmit

import glob

# find all output root files
finished_files = glob.glob('*.root')
finished_index = [ item.split('output_')[1].split('.root')[0] for item in finished_files]

# read original ana file
fin = open('ana.listOfJobs')
jobs = fin.readlines()
fin.close()
# create a new ana file
fout = open('ana.listOfJobs','w')

# find the input file txt
inputfile_txt = jobs[0].split('--txtfiles tardir/')[1].split('--')[0].strip()

# # find the number of files per job
# maxfiles = jobs[0].split('--maxfiles')[1].split('--')[0].strip()

# read input file dirs
with open(inputfile_txt,'r') as tmp:
	all_input = tmp.readlines()
# Get a list of input root file index 
all_input.pop(0)
all_input_index = [item.split('numEvent')[1].split('_')[1].split('.root')[0] for item in all_input]

# select job that has no output seen
for i,line in enumerate(jobs):
	# find first root file index for current job 
	start_file = line.split('--startfile')[1].split('--')[0].strip()
	start_file = int(start_file)
	start_file_index = all_input_index[start_file]

	if start_file_index not in finished_index:
		fout.write(line)
		print 'input file %s not in outputs, start file %s, %sth job'%(start_file_index,start_file,i)
		print 

# finish
fout.close()



