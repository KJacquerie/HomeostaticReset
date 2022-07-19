# Description of the Matlab codes used to generate the images of the Article


* * Datasheet_plot_STDP.m *
Function that allows to obtain the STDP curves plotted in the datasheet of all models
It used the function * function_fig_STDP.m * to generate the graph.
> Resulting figures are available in the folder Fig_datasheet_STDP

* *Datasheet_plot_wf.m*
Function that allows to obtain the frequency-dependency in the pairing protocol [Sjostrom,2001]. It shows the ration between the final weight and the initial weight as a function of the pairing frequency (wf)
It used the function * function_fig_wf.m * to generate the graph and SJO.mat the experimental data.
> Resulting figures are available in the folder Fig_datasheet_wf.


* *Datasheet_reset.m*
It shows the homeostatic reset for all models coded in soft-bounds weight dependency and saturation in models implementing the hard bounds dependency. 
It uses the function *function_plot_wtime_Supp.m* to generate the graph.
> Resulting figures are available in the folder Fig_datasheet_reset.

* *DEMO_SAT_simu.m*
This code permits to compare and plot the slope of saturation computed on the simulation and the analytical solution in phenomenological models.

* *DEMO_wSAT_simu.mat*
It contains the variable and the parameters to compute the slope of saturation in phenomenological models.

* *Fig1_plot_scenario.m*
Matlab code to obtain Figure 1 providing the evolution of the synaptic weights in 6 circuits. 

* *Fig2_plot_wf_superpose.m*
Compute the mean and the standard deviation of all the frequency-dependency in the pairing protocol [Sjostrom,2001] computed in each model (implementing soft bounds).
It uses the function *function_fig_wf.m*
—> Resulting figures are available in the fold Fig2.

* *Fig2_plot_STDP_superpose.m*
Similar as the previous code but it compares the STDP results in each model (implementing hard bounds). 
It uses the function *function_fig_STDP.m*

* *Fig3_trace.m (trace V en A, B)*
It provides the membrane potential evolution (A) for different applied currents (data are available in the folder Calcium/model2/trace/burst/...) and (B) for intrinsic variabilities (data are available in the folder Calcium/model2/results_Variab/...)
—> Resulting figures are available in the folder Fig3/Fig3_trace

* *Fig3AB_wtime_burst.m*
It provides the homeostatic reset shown in the article (corresponding to model 2) and saturation (corresponding to model 9). The other results are shown in Fig_datasheet_reset
It uses the function *function_plot_wtime_Article.m* to display the graph.
—> Resulting figures are available in the folder Fig3/Fig3_wtime

* *Fig3B_varibility.m*
It computes the mean of the synaptic weight evolution computed in 10 circuits with variabilities. 
> moyenne les données des 10 circuits pour un courant 
> avec l’évolution de la moyenne au cours du temps pour les diff w0
—> save dans Fig3_wtime

* *Fig3C_Network.m (plot w(t) et les traces)*
It uses the data from the Julia .jl files/Network_heterogeneous/traceX (X indicates the number of the considered network).
--> Resulting figures are available in Fig3/Fig3_ntk.

* *Folder Fig 4*
display_correlation.m: plot the correlation matrix
Fig4_scatter.m: display the scatter plot. 
scatter_simuvspredict.m: function to retrieve the simulated value of the reset or the predicted value.

* *Fig4_wtime_burst *
[appendix figure]
Function to obtain the homeostatic reset and the saturation curve shown in Figure 4.
It uses the function *function_plot_wtime_Oneline.m* to display only one curve (from the Figure 3). 
--> Resulting figures are available in Fig4

* * Fig4_temporel_calcium.m*
[appendix figure]
Function to obtain the zoom of the calcium traces as well as a zoom in the synaptic weight evolution to explain intuitively the homeostatic reset mechanism in calcium-based models.
--> It saves the resulting figures in Fig4

* *Fig4_temporel_phenom.m*
[appendix figure]
Function to obtain the zoom of the pre and postsynaptic traces as well as a zoom in the synaptic weight evolution to explain intuitively the homeostatic reset mechanism in phenomenological models.
--> It saves the resulting figures in Fig4

* *SI_plot_Vburst.m*
Function to provide the membrane potential evolution for 8 different levels of NMOD.
--> It saves the resulting figures in Fig_SI

* *SI_wf.m*
Function that compares the frequency-rate in a pairing protocol for 75 pulses (as done in [Graupner,2016]) and 7x15 with a period of silence (as done in [Sjostrom,2001]).
--> It saves the resulting figures in Fig_SI/75vsSJO.eps

* *plot_V.m*
Function that plots the evolution of the membrane potential. For example, it was used to show the frequency-dependency pairing protocol in Fig2
—> Resulting figures are saved in the related Fig folders. 