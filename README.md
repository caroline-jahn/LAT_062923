# LAT_062923
Codes for "Learning attentional templates for value-based decision-making"

Caroline Jahn - 06/29/23

For now, the raw data is not available. We will deposit the data and the codes before publication.

Set the path!
You need to download:
- VBMC toolbox (https://github.com/lacerbi/vbmc)
- NPMK toolbox (https://github.com/BlackrockNeurotech/NPMK)
- readPLXFileC and Matlab Offline Files SDK (Plexon software)
- CircStat2012a toolbox (https://github.com/anne-urai/Tools/tree/master/CircStat2012a)
- Violinplot-Matlab toolbox (https://github.com/bastibe/Violinplot-Matlab)

Behavioral data (Fig1)

In Codes_Figure_1

doExecuteMasterScriptResetRWVBMC(fsroot,monkey,N_channels,analysis) (winning model: N_channels=6, analysis=‘Reset’)
save_Model_VBMC

plot_figure1 (behavioral results per monkey)
plot_models_VBMC (model comparison)
plot_models_VBMC_2L (2 learning rate model comparison)

Neural data (Fig2-6)

Pre-processing (get the FR for each neuron at each event)

1) Get the channels for each session (based on preprocessing on Plexon offline sorter)
In Split_channels
Run chanListCompile.m
Saved in /Volumes/buschman/Projects/Learning_Attentional_Templates/Data/General/NS6Directory.mat 

2) Set what the start of the session was (semaphore) for the delayed saccade task
ppParseTrials_delsaccade
Run chanListCompile.m
Saved in /Volumes/buschman/Projects/Learning_Attentional_Templates/Data/General/NS6Directory_sem.mat
NOTE: this contains the semaphores for delayed saccade task and will be overwritten when we run it for the main task

3) Align channel mua activity to trial event for spikes in delayed saccade task task for each SD 
In ppParseTrials_delsaccade
Run cutTrialsSave.m
Saved in /Volumes/buschman/Projects/Learning_Attentional_Templates/Data/Monkey/date/MUA/xxxsd

4) Get the FR for each trial in time window for the delayed saccade task
In spike_count_delsaccade
Run stimTypeFiringRate.m
Saved in /Volumes/buschman/Projects/Learning_Attentional_Templates/Data/Monkey/date/MUA/stimTypeFR

5) Set what the start of the session was (semaphore) for the main task
In ppParseTrials_exploreexploit
Run chanListCompile.m
Saved in /Volumes/buschman/Projects/Learning_Attentional_Templates/Data/General/NS6Directory_sem.mat

6) Align channel activity to trial event for spikes in main task (exploreexploit)
In /pParseTrials_exploreexploit
Run cutTrialsSave.m
Saved in /Volumes/buschman/Projects/Learning_Attentional_Templates/Data/Monkey/date

7) Set which channels belong to which ROI and were active
In Grid
Run ChannelsActivityToGrid.m
Saved in /Volumes/buschman/Projects/Learning_Attentional_Templates/Data/General/NS6Directory_AC.mat

8) Find out which channels were responding to direction of prepared saccade (delayed saccade task, MUA) for each SD threshold
In /stats_delsaccade
Run compute_stats.m
Saved in /Volumes/buschman/Projects/Learning_Attentional_Templates/Analysis/Electrophy_analysis/delsacc/xxxsd

9) Map the sensitive channels (set the SD for the MUA, for paper)
In Vstats_delsaccade
Run Save_LIP_sensitive_channels.m 
Saved in /Volumes/buschman/Projects/Learning_Attentional_Templates/Data/General/NS6Directory_SC.mat

10) Get the FR for each trial in time window
In spike_count_exploreexploit
Run extract_data_for_belief_analysis_Reset_RW.m or extract_data_for_belief_analysis_reward_end_RW.m
Saved in /Volumes/buschman/Projects/Learning_Attentional_Templates/Analysis/Electrophy_analysis/exploreexploit/Reset_RW_model

Figure 2 and S2

In spike_count_exploreexploit

Create the window [-600 300] around target, size = 900ms (step 10 in pre-processing)
Create the sliding window around target, -600 to 350ms, size = 300ms
Create the sliding window around reward, 0 to 300ms, size = 300ms 

In Codes_Figure_2

models_FR_figure2_value_function.m (for window [-600 300])
plot_FigS2_prop_sig_models.m

prepare_pseudopop_peak_belief_classifier_prog_restricted_FEF.m 
prepare_bootstrap_peak_belief_classifier_across_time_prog.m

Run_belief_classifier_combined_time_prog_restricted_FEF.m (for window [-600 300])
Run_belief_classifier_combined_time_prog_NN_restricted_FEF.m (for window [-600 300])
get_projection_classifier_combined_time_prog_restricted_FEF.m (for window [-600 300])
plot_Fig2_ET_pseudo_pop.m

Run_belief_classifier_combined_time_prog_restricted_FEF.m (for sliding windows)
Run_belief_classifier_combined_time_prog_NN_restricted_FEF.m (for sliding windows)
get_projection_classifier_combined_time_prog_restricted_FEF.m (for sliding windows)

plot_Fig2_cross_temporal_decoding.m

plot_Fig2_PCA.m

models_FR_figure2_true_template.m (for window [-600 300])
plot_FigS2_true_template_prop_sig.m

models_FR_figure2_value_function.m  (for sliding windows)
plot_FigS2_prop_sig_models_across_time.m

Figure 3, S3 a and S4

In Codes_Figures_3_and_4

Run belief_peak_decoder_single_session_update_NN_restricted_FEF.m
Plot_fig3_and_4.m

Figure 5 and S5

In spike_count_exploreexploit
Create the sliding window around target, -500 to 800ms, size = 200ms
Create the sliding window around response, -500 to 600ms, size = 200ms

In Codes_Figure_5

prepare_pseudopop_choice_prog_classifier_across_time.m
prepare_bootstrap_peak_belief_choice_across_time.m
Run_classifier_peak_belief_choice_restricted_FEF.m
Accuracy_classifier_peak_belief_choice.m

prepare_pseudopop_peak_belief_value_classifier_across_time.m
prepare_bootstrap_peak_belief_value_across_time.m
Run_classifier_peak_belief_value_restricted_FEF.m
Accuracy_classifier_peak_belief_value.m

prepare_pseudopop_peak_belief_cc_classifier_across_time.m
prepare_bootstrap_peak_belief_cc_across_time.m
Accuracy_classifier_peak_belief_chosen_color_2bins.m

prepare_pseudopop_stim_color_across_time.m
prepare_bootstrap_stim_color_across_time.m
Run_stim_color_classifier_across_time_with_NN_restricted_FEF.m
Accuracy_classifier_stim_color.m

Figure 6 and S6

In Codes_Figure_6

get_explained_variance_split_restricted_FEF.m 
get_boot_glm_all_values_split_correct.m
get_mean_EV_across_split_and_locs.m

plot_correlation_chosen_unchosen.m
plot_correlation_local_global.m
plot_correlation_across_locations.m
plot_EV_across_loc.m








