import csv

input = open('summary_combined.csv')
results = list(input)
sys_type_list = ['nominal_1M','btag_eff_reweight','lepID_reweight','lepIso_reweight','tracking_reweight','trigger_reweight']
fit_result_lists = (len(sys_type_list)*2-1)*['']
uncertainty_list = []
fit_parameter = ['Rqqbar','Rbck','Delta','Afb']
uncertainty_list.append(['corrections']+fit_parameter)

for line in results :
	#line = [word.strip() for word in line] 
	if sys_type_list[0] in line :
		fit_result_lists[0] = line
	i = 0
	for types in sys_type_list :
		type_hi = types+'_hi'
		type_low = types+'_low'
		if type_hi in line :
			fit_result_lists[2*i-1] = line
		if type_low in line :
			fit_result_lists[2*i] = line
		i += 1
for i in range( len(sys_type_list) ) :
	print i
	if i == 0 :
		continue
	if i>0 :
		tmp_nom = fit_result_lists[0].split(',')
		tmp_high = fit_result_lists[2*i-1].split(',')
		tmp_low = fit_result_lists[2*i].split(',')
		tmp_uncertainty = []
		tmp_uncertainty.append( sys_type_list[i] )
		tmp_uncertainty.append( (abs(float(tmp_high[1])-float(tmp_nom[1]))+abs(float(tmp_low[1])-float(tmp_nom[1])))/2)
                tmp_uncertainty.append( (abs(float(tmp_high[3])-float(tmp_nom[3]))+abs(float(tmp_low[3])-float(tmp_nom[3])))/2)
                tmp_uncertainty.append( (abs(float(tmp_high[7])-float(tmp_nom[7]))+abs(float(tmp_low[7])-float(tmp_nom[7])))/2)
                tmp_uncertainty.append( (abs(float(tmp_high[9])-float(tmp_nom[9]))+abs(float(tmp_low[9])-float(tmp_nom[9])))/2)
		uncertainty_list.append(tmp_uncertainty)

output = open('detector_sys_combined.csv','w') 
csvwriter = csv.writer(output)
for row in range(len(uncertainty_list)):
	csvwriter.writerow(uncertainty_list[row])
input.close()
output.close()

