for i in range(0,53):
	output = 'hist_PDF_'+str(i)+'_output'
	print 'cp -r std_grid_fit_combined/*.* '+output
	print 'cd '+output
	print 'source grid_sub.csh'
	print 'cd ..'
