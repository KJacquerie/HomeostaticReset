clear all
close all
clc

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');


%%
function_plot_wtime_Supp('model1','Gamma','model1', 'weak') %#1

function_plot_wtime_Supp('model2','Control','model2', 'weak') % #2

function_plot_wtime_Supp('model3','kernel_DP','model3', 'weak') % #3
%%
function_plot_wtime_Supp('model4','cortex','model4', 'weak') % #4
function_plot_wtime_Supp('model5','DP','model5', 'weak') % #5

function_plot_wtime_Supp('model6','CC','model6', 'weak') % #6


function_plot_wtime_Supp('model7','Dep','model7', 'weak') % #7

function_plot_wtime_Supp('model8','STD','model8', 'weak') % #8
%%

function_plot_wtime_Supp('model9','DD_rel_half','model9', 'weak') % #9
function_plot_wtime_Supp('model10','CC_rel','model10', 'weak') % #10
%%
function_plot_wtime_Supp('model11','Dep_DD_rel','model11', 'weak') % #11

function_plot_wtime_Supp('model12','Dep_CC_rel','model12', 'weak') % #12
%%

%function_plot_wtime_Supp('STDxShouval','DD','model8bis', 'weak') % =model 8

function_plot_wtime_Supp('model13','STD_DD_rel','model13', 'weak') %  # 13
%%
function_plot_wtime_Supp('model14','STD_CC_rel','model14', 'weak') %  # 14

%%


% 
% %% %%%% COSYNE %%%%%
% function_plot_wtime_Oneline_Cosyne('GB_DD','AllWdep_Depression_thetaHIGH', 'weak', 2)
% function_plot_wtime_Cosyne('GB_DD','OmegaWdep_Born', 'weak', 3)
% function_plot_wtime_Oneline_Cosyne('Caro_Gonz','T10', 'cosyne', 6)
% function_plot_wtime_Oneline_Cosyne('GB_DD','Control', 'weak', 6)
% function_plot_wtime_Oneline('GB_DD','Control', 'weak', 6)
% function_plot_wtime_Oneline('GBxShouval','DD_rel_half', 'weak', 7)
% function_plot_wtime_Oneline('GB_DD','AllWdep_Depression_thetaHIGH', 'weak', 2)
% function_plot_wtime_Oneline('Caro_Gonz','T10', 'cosyne', 6)
% function_plot_wtime_Article('GB_DD','OmegaWdep_Born', 'weak', 3)
% 
% 
% 
% 
% 
% 
% 
% %% %%%%%%%%% NMOD %%%%%%%%
% function_plot_wtime('GB_DD','OmegafctW', 'weak')
% 
% function_plot_wtime('GB_kernel','kernel_DP', 'strong')
% 
% function_plot_wtime('GB_kernel','kernel_D', 'weak')
% function_plot_wtime('GB_kernel','kernel_D', 'strong')
% 
% function_plot_wtime('GB_kernel','kernel_Ddown', 'weak')
% function_plot_wtime('GB_kernel','kernel_Ddown', 'strong')
% 
% %%
% function_plot_wtime('GB_gonz','Gonz_DU', 'weak')
% function_plot_wtime('GB_gonz','Gonz_DU', 'strong')
% 
% function_plot_wtime('GB_gonz','Gonz_DP_thetaLOW', 'weak')
% function_plot_wtime('GB_gonz','Gonz_DP_thetaLOW', 'strong')
% 
% 
% function_plot_wtime('GB_gonz','Gonz_DP', 'weak')
% function_plot_wtime('GB_gonz','Gonz_DP', 'strong')
% 
% function_plot_wtime('GB_gonz','Gonz_DPP', 'weak')
% function_plot_wtime('GB_gonz','Gonz_DPP', 'strong')
% 
% 
% function_plot_wtime('GB_DD','thetaHIGH_OmegaLOW', 'weak')
% function_plot_wtime('GB_DD','thetaHIGH_OmegaLOW', 'strong')
% 
% function_plot_wtime('GB_DD','OmegaWdep', 'weak')
% function_plot_wtime('GB_DD','OmegaWdep', 'strong')
% 
% function_plot_wtime('GB_DD','OmegaWdep_BORN', 'weak')
% function_plot_wtime('GB_DD','OmegaWdep_BORN', 'strong')
% 
% function_plot_wtime('GB_DD','OmegaWdep_BORN_thetaHIGH', 'weak')
% function_plot_wtime('GB_DD','OmegaWdep_BORN_thetaHIGH', 'strong')
% 
% function_plot_wtime('GB_DD','AllWdep_Born_thetaLOW', 'weak')
% function_plot_wtime('GB_DD','AllWdep_Born_thetaLOW', 'strong')
% 
% 
% function_plot_wtime('GB_DD','AllWdep_Born_thetaHIGH', 'weak')
% function_plot_wtime('GB_DD','AllWdep_Born_thetaHIGH', 'strong')
% 
% function_plot_wtime('GB_DD','AllWdep_Depression_thetaHIGH', 'weak')
% function_plot_wtime('GB_DD','AllWdep_Depression_thetaHIGH', 'strong')
% 
% 
% function_plot_wtime('GB_DD','4levels_LOW', 'weak')
% function_plot_wtime('GB_DD','4levels_LOW', 'strong')
% 
% function_plot_wtime('GB_DD','4levels_HIGH', 'weak')
% function_plot_wtime('GB_DD','4levels_HIGH', 'strong')
% 
% %%
% 
% function_plot_wtime('GB_DD','TAG', 'weak')
% function_plot_wtime('GB_DD','TAG', 'strong')
% 
% %%
% 
% 
% 
% function_plot_wtime('GB2012','cortex', 'strong')
% %%
% function_plot_wtime('GB2012','DP', 'weak')
% function_plot_wtime('GB2012','DP', 'strong')
% 
% function_plot_wtime('GB2012','DP_pink', 'weak')
% function_plot_wtime('GB2012','DP_pink', 'strong')
% 
% function_plot_wtime('GB2012','Donly', 'weak')
% function_plot_wtime('GB2012','Donly', 'strong')
% 
% function_plot_wtime('GB2012','Dprime', 'weak')
% function_plot_wtime('GB2012','Dprime', 'strong')
% 
% %%
% function_plot_wtime('GBxShouval','DD', 'weak')
% function_plot_wtime('GBxShouval','DD', 'strong')
% %%
% function_plot_wtime('GBxShouval','CC', 'weak')
% function_plot_wtime('GBxShouval','CC', 'strong')
% %%
% function_plot_wtime('GBxShouval','CC_rel', 'weak')
% function_plot_wtime('GBxShouval','CC_rel', 'strong')
% %%
% function_plot_wtime('GBxShouval','DD_rel', 'weak')
% function_plot_wtime('GBxShouval','DD_rel', 'strong')
% 
% %%
% function_plot_wtime('GBxShouval','DD_rel_half', 'weak')
% function_plot_wtime('GBxShouval','DD_rel_half', 'strong')

%%

% 
% function_plot_wtime_Oneline('GB_DD','Control','_RESET_Oneline', 'weak', 6) % #1
% close all 
% function_plot_wtime_Oneline('GBxShouval','DD_rel_half','_SATURATION_Oneline', 'weak', 7) % #9
% close all
% %%
% function_plot_wtime_Oneline('GB_DD','Var','_VARIABILITY_Oneline', 'weak', 6) % #1 with variability in intrinsic parameters
% 
% 
% %%
% close all
% function_plot_wtime_Article('GB_DD','Control','_RESET', 'weak', 6) % #1
% %%
% close all
% function_plot_wtime_Article('GBxShouval','DD_rel_half','_SATURATION', 'weak', 3) % #9
% close all
% function_plot_wtime_Article('GB_DD','Var','_VARIABILITY', 'weak', 6) % #1 with variability in intrinsic parameters
% 
% %%
% 
% function_plot_wtime_variability('_VARIABILITY', 1, 4)
% %%
% function_plot_wtime_Network('_network1', 1)
% 
% %%
% 
% 
% %%
% function_plot_wtime_Network('Network_heterogeneous', 1)
% function_plot_wtime_Network('Network_heterogeneous', 2)
% function_plot_wtime_Network('Network_heterogeneous', 3)
% function_plot_wtime_Network('Network_heterogeneous', 4)



