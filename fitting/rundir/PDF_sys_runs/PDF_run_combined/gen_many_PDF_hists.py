import os
num_PDFs = 53 
for i in range(num_PDFs):
	cmd = 'python mc_run.py --n hist_PDF_'+str(i)+' --pdf '+str(i)
	os.system(cmd)

