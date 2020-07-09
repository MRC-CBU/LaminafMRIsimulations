# LaminafMRIsimulations
Based on PhD work of Huang Pei (HP)

## List of files
### Simulation files
attention_simulation.m - main function used to generate the simulation results. Calls the following subfunctions: attention_simulation_iteration.m, covdiag.m, create_glm.m, deming.m and detrend_timecourse.m.
attention_simulation_plots.m - Searches the list of outputs from attention_simulation.m for the corresponding examples and generates the graphs.
attention_simulation_contri_plots.m - Similar to attention_simulation_plots.m, but plots the relative contributions of laminar bias and attentional modulation on the data.

### Real data files
attention_plots.m - Plots the results from real_fMRI_results.mat in a similar fashion to attention_simulation_plots.m.
tstat_comparison.m - Outputs the numerical comparison and graphs for the tstat of real and simulated data for matching purposes

### Sample results directory
att_sim_results.mat - Results of running attention_simulation.m with n_neurons_face_sigma = 1.1 and n_neurons_house_sigma = 0.5
Standard/ - Graphical outputs from att_sim_results.mat
sim_scatter/ - Scatterplots of TaskD+/TaskD- for the standard condition

att_sim_results_samepref.mat - Results of running attention_simulation.m with n_neurons_face_sigma = 1.1 and n_neurons_house_sigma = 1.1
SamePref/ - Graphical outputs from att_sim_results_samepref.mat

real_fMRI_results.mat - Results from a real fMRI study for comparison. In addition to summary metrics for comparison, raw tstat values and contrast estimates for all voxels is available under results_exclude.full
real_data/ - Graphical outputs from the real fMRI data
