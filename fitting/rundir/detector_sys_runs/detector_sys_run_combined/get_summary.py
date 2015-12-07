f = open('summary_list.txt')
result = list(f)
result = [word.strip() for word in result]
cmd = ''
for run in result :
	cmd += 'cat '+run+' >> summary_combined.csv \n'
print cmd
