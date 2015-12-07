import os
# Decide how many PDF members there are in the sample root files
PDF_members = 53 

f_run_MC = open("run_MC_jobs.sh","w")
f_run_data_separated = open("run_data_separated.sh","w")
f_run_data_combined = open("run_data_combined.sh","w")
f_clean_MC = open("clean_MC.sh","w")
f_clean_data = open("clean_data.sh","w")
run_MC_cmd = ''
run_data_separated_cmd = ''
run_data_combined_cmd = ''
clean_MC_cmd = ''
clean_data_cmd = ''
for i in range(PDF_members) :
	#Make ana.listOfJobs for each PDF
        f_job = open("ana.listOfJobs","w")
        job_cmd =  'root -l -b -q tardir/main_MC.C(0,1,"test1","tardir/MC_input_with_bkg.txt",'+str(i)+')'
        f_job.write(job_cmd)
        f_job.close()
	#Make the run dir and copy the ana.listOfJobs files just made
	cmd = 'cp standard_grid_run PDF_'+str(i)+' -r \n'
	cmd += 'cp ana.listOfJobs PDF_'+str(i) 
	#print cmd 
        os.system(cmd)
	#Commmand to run all jobs
	run_MC_cmd += 'cd PDF_'+str(i)+'\nsource grid_sub_MC.csh\ncd ..\n'
        run_data_separated_cmd += 'cd PDF_'+str(i)+'\nsource grid_sub_separated.csh\ncd ..\n'
        run_data_combined_cmd += 'cd PDF_'+str(i)+'\nsource grid_sub_combined.csh\ncd ..\n'
	#Command to clean MC or data grid files
	clean_MC_cmd += 'cd PDF_'+str(i)+'\nsource cleanup_MC.sh\ncd ..\n'
        clean_data_cmd += 'cd PDF_'+str(i)+'\nsource cleanup_grid.sh\ncd ..\n'
#write commands to run all jobs
f_run_MC.write(run_MC_cmd)
f_run_data_separated.write(run_data_separated_cmd)
f_run_data_combined.write(run_data_combined_cmd)
f_clean_data.write(clean_data_cmd)
f_clean_MC.write(clean_MC_cmd)
