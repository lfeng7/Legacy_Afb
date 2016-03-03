mkdir data_histo_files
mkdir final_fit_stuff
mkdir fit_results
mv *_data_histos.root data_histo_files
mv f_*_distribution.pdf final_fit_stuff
mv all_data_histos*.* final_fit_stuff
mv all_data_angles.* final_fit_stuff
mv fit_comparison.pdf final_fit_stuff
mv sideband.pdf final_fit_stuff
mv data_run_combined_output/final_fit_stuff/fit_results.txt fit_results/
mv fit_comparison.png final_stack.root  fit_results/
mv fit_results data_run_combined_output
