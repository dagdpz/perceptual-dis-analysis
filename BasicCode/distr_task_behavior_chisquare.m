function stat_struct_chi2 = distr_task_behavior_chisquare(data_struct_per_session, stim_con, Settings)

trial_cons = {'targ_targ','distr_distr','targ_distr','single_targ','single_distr'};
avs = {'contra_choice','fixation_choice','ipsi_choice'};
avnames = {'num_choice_contra','num_choice_fixation','num_choice_ipsi'};
comparisons = [1 2; 1 3; 1 4];

for tcon = 1:length(trial_cons) % loop through trial conditions
    
    stimulus_sides = fieldnames(data_struct_per_session.(trial_cons{tcon}));
    
    for stimside = 1:length(stimulus_sides) % loop through side on which stimulus/stimuli was/were presented 
        
        for av = 1:length(avs) % loop through dependent variables
            
            ncolors = size(data_struct_per_session.(trial_cons{tcon}).(stimulus_sides{stimside}).(stim_con{1}),1);
            
            for col = 1:ncolors % loop through stimulus colors
                
                for ses = 1:size(data_struct_per_session.(trial_cons{tcon}).(stimulus_sides{stimside}).(stim_con{1}),2) % loop through sessions
                    
                    for c = 1:size(comparisons, 1)
                        
                        % create vector with binary variable for each trial defining whether the respective choice was made
                        stat_struct_chi2.(trial_cons{tcon}).(stimulus_sides{stimside}).(avs{av}){col,ses}{c,1} = ...
                            [ones(1,data_struct_per_session.(trial_cons{tcon}).(stimulus_sides{stimside}).(stim_con{comparisons(c,1)}){col,ses}.(avnames{av})),...
                            zeros(1,data_struct_per_session.(trial_cons{tcon}).(stimulus_sides{stimside}).(stim_con{comparisons(c,1)}){col,ses}.num_total_trials - ...
                            data_struct_per_session.(trial_cons{tcon}).(stimulus_sides{stimside}).(stim_con{comparisons(c,1)}){col,ses}.(avnames{av})),...
                            ones(1,data_struct_per_session.(trial_cons{tcon}).(stimulus_sides{stimside}).(stim_con{comparisons(c,2)}){col,ses}.(avnames{av})),...
                            zeros(1,data_struct_per_session.(trial_cons{tcon}).(stimulus_sides{stimside}).(stim_con{comparisons(c,2)}){col,ses}.num_total_trials - ...
                            data_struct_per_session.(trial_cons{tcon}).(stimulus_sides{stimside}).(stim_con{comparisons(c,2)}){col,ses}.(avnames{av}))];
                        
                        % create vector with binary variable for each trial defining which stimulation condition the trial belongs to
                        stat_struct_chi2.(trial_cons{tcon}).(stimulus_sides{stimside}).(avs{av}){col,ses}{c,2} = ...
                            [repmat(comparisons(c,1),1,data_struct_per_session.(trial_cons{tcon}).(stimulus_sides{stimside}).(stim_con{comparisons(c,1)}){col,ses}.num_total_trials),...
                             repmat(comparisons(c,2),1,data_struct_per_session.(trial_cons{tcon}).(stimulus_sides{stimside}).(stim_con{comparisons(c,2)}){col,ses}.num_total_trials)];

                        % add stimulation conditions that are compared to each other to output structure
                        stat_struct_chi2.(trial_cons{tcon}).(stimulus_sides{stimside}).(avs{av}){col,ses}{c,3} = stim_con{comparisons(c,1)}; 
                        stat_struct_chi2.(trial_cons{tcon}).(stimulus_sides{stimside}).(avs{av}){col,ses}{c,4} = stim_con{comparisons(c,2)}; 
                           
                        [tbl,chi2,p,labels] = crosstab(stat_struct_chi2.(trial_cons{tcon}).(stimulus_sides{stimside}).(avs{av}){col,ses}{c,1},stat_struct_chi2.(trial_cons{tcon}).(stimulus_sides{stimside}).(avs{av}){col,ses}{c,2});
                        
                        if p <= 0.05
                            stat_struct_chi2.(trial_cons{tcon}).(stimulus_sides{stimside}).(avs{av}){col,ses}{c,5} = 1;
                        else stat_struct_chi2.(trial_cons{tcon}).(stimulus_sides{stimside}).(avs{av}){col,ses}{c,5} = 0;
                        end
                        
                        stat_struct_chi2.(trial_cons{tcon}).(stimulus_sides{stimside}).(avs{av}){col,ses}{c,6} = p;
                        stat_struct_chi2.(trial_cons{tcon}).(stimulus_sides{stimside}).(avs{av}){col,ses}{c,7} = chi2;
                        stat_struct_chi2.(trial_cons{tcon}).(stimulus_sides{stimside}).(avs{av}){col,ses}{c,8} = (numel(unique(stat_struct_chi2.(trial_cons{tcon}).(stimulus_sides{stimside}).(avs{av}){col,ses}{c,1}))-1)*...
                            (numel(unique(stat_struct_chi2.(trial_cons{tcon}).(stimulus_sides{stimside}).(avs{av}){col,ses}{c,2}))-1);
                    end                                    
                end               
            end
        end    
    end   
end

stat_struct_chi2.sessions = Settings.folders;

if Settings.save_stat_output
    save([Settings.save_dir.STIM_gh4_5 filesep 'stat_struct_chi2.mat'], 'stat_struct_chi2');
    disp(['Saved ' Settings.save_dir.STIM_gh4_5 filesep 'stat_struct_chi2.mat']);
end