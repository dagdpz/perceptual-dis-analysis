function [ stat_results_anova, stat_struct_ttest, stat_struct_ttest_anova ] = distr_task_behavior_statistics( data_struct_per_session, Settings, stim_con )
% Statistical analysis of the distractor task.
% Input: - Data structures from distr_task_behavior_analysis
%        - Settings from distr_task_behavior_callscript
%        - stim_con variable. A cell array containing strings to identify 
%             the fieldnames representing the stimulation conditions in the
%             data structure
%
% Statistics across sessions
%
% Although the statistics are calculated across sessions, the input for 
% this test are the respective means per session, so 
% Settings.analyze_across_sessions does not need to be turned on.
%
% CHOICE:
% - T-Tests of all plots comparing choices among stimulation conditions 
%   (4 conditions -> 6 comparisons per choice & trial condition).
%   Correction for multiple testing using the False Discovery Rate approach
%   not functional yet.
% - 3 way repeated measures ANOVA with stimulus type (contra-targ
%   ipsi-distr, contra-distr ipsi-targ, contra-distr ipsi-distr),
%   stimulation onset (off, -80, Go +80) and difficulty (easy and
%   difficult distr color) as factors. 
%   - Post hoc: ttest of stimulation conditions
%     Paired ttest testing the null hypothesis that the pairwise difference 
%     between stimulation conditions has a mean equal to zero.
% - 3 way repeated measure ANOVA with single stimulus type (single targ, 
%   difficult single distr, easy single distr), side (contraversive or 
%   ipsiversive stimulus display) and stimulation (off, -80ms, Go +80ms) as
%   factors. 
%
% LATENCY:
% - 2 way repeated measures ANOVA of single difficult distractor trials / 
%   single target trials with side (contraversive & ipsiversive distractor)
%   and stimulation (off, -80ms, Go +80ms) as factors. 
% - 3 way repeated measures ANOVA of correct and error saccades in double 
%   stimulus trials (correct saccades: target-distractor & target-target;
%   error saccades: target-distractor & distractor-distractor) with choice
%   (contraversive saccade & ipsiversive saccade) and stimulation condition
%   (off, -80ms, Go +80ms) as 2nd and 3rd factors. This design is 
%   calculated separately for early and difficult distractor trials and 
%   separately for correct and error saccades. 

% stim_con = {'stim_off','minus_80ms','go_signal','plus_80ms'};
stat_results_anova = [];
stat_struct_ttest = [];
stat_struct_ttest_anova = [];
num_stim_con = numel(stim_con);
%%

%-----------------------------------------------------------------
    %%% CHOICES, LATENCY AND ERROR RATES FOR PLOTS (TTESTS) %%%
%-----------------------------------------------------------------

%%%  T-TESTS FOR SIGNIFICANCE BARS IN BAR PLOTS %%%
% Test the null hypothesis that the pairwise difference between data vectors x and y has a mean equal to zero.
if Settings.ttests
    perform_ttests = 1;
    tr_con2stat = {'targ_targ','distr_distr','targ_distr','targ_distr','single_targ','single_targ','single_distr','single_distr'};
    tr_con2stat = {'targ_targ_1HF','targ_targ_2HF','distr_distr_1HF','distr_distr_2HF','targ_distr_1HF','targ_distr_1HF','targ_distr_2HF','targ_distr_2HF','single_targ','single_targ','single_distr','single_distr'};
    num_tr_con2stat = numel(tr_con2stat);
    position = {'combined','combined','contra','ipsi','contra','ipsi','contra','ipsi'}; % refers to tr_con2stat
    %    position = {'combined','combined','left','right','left','right','left','right'}; % refers to tr_con2stat
    cho_lat_err = [Settings.analyze_choice,Settings.analyze_latency,Settings.analyze_error_rates];
    if Settings.FDR
        mode = 'compare_to_control';
    else
        mode = 'all_comparisons';
    end
    
    for tr_con = 1:num_tr_con2stat
        tr_con
        for pos = tr_con; % for the respective position in this trial condition
            for c_l_e = 1:numel(cho_lat_err) % CHOICE OR LATENCY OR ERRORS
                if cho_lat_err(c_l_e)
                    if c_l_e == 1 %%% CHOICES %%%
                        if tr_con > 4  && strcmp(position{pos},'contra') % 'single_targ' 'single_distr' conditions with contra stimulus
                            dep_var = {'contra','fixation'}; % don't create empty ipsi-choices output when only one contra stimulus is shown
                        elseif tr_con > 4  && strcmp(position{pos},'ipsi') % 'single_targ' 'single_distr' conditions with ipsi stimulus
                            dep_var = {'fixation','ipsi'}; % don't create empty contra-choices when only one ipsi stimulus is shown
                        else
                            dep_var = {'contra','fixation','ipsi'};
                        end
                    elseif c_l_e == 2 %%% LATENCY %%%
                        dep_var = {'contra','ipsi'};  %{'left','right'}; %KK changed
                    elseif c_l_e == 3 %%% ERRORS %%%
                        dep_var = {'fixation_break','target_acquis','target_hold'};
                    end
                    for dv = 1:numel(dep_var) % number choice options
                        if c_l_e == 1 %%% CHOICES %%%
                            fname_stat = [dep_var{dv},'_choice']; % fieldname for stat_data structure
                            fname_tmp = 'true_prop_choice_per_session_';
                            fname = [fname_tmp,dep_var{dv}]; % fieldname for tmp_data structure
                        elseif c_l_e == 2 %%% LATENCY %%%
                            fname_stat = [dep_var{dv},'_latency']; % fieldname for stat_data structure
                            fname = ['mean_latency_',dep_var{dv},'_choice_per_session'];
                        elseif c_l_e == 3 %%% ERRORS %%%
                            fname_stat = [dep_var{dv},'_errors']; % fieldname for stat_data structure
                            fname = ['error_rate_',dep_var{dv},'_per_session']; % fieldname for tmp_data structure
                        end
                        for stim = 1:num_stim_con
                            for col = 1:size(data_struct_per_session.(tr_con2stat{tr_con}).(position{pos}).(stim_con{stim}),1); % for all colors (rows in data_struct) in this trial condition
                                tmp_data = [data_struct_per_session.(tr_con2stat{tr_con}).(position{pos}).(stim_con{stim}){col,:}]; % tmp_data with all sessions concatenated
                                if isequal(unique([tmp_data.G_R_ratios])',[0.0859 ,1]) || isequal(unique([tmp_data.G_R_ratios]),[0])%KK changed 0.1796875 check if distractor colors are the easy and difficult distractor in all sessions, col = 1 is always difficult distr, col = 2 always easy distr
                                    % add the data to the stat_struct_ttest
                                    stat_struct_ttest = get_stat_struct_ttest (mode, stat_struct_ttest, stim, tr_con2stat{tr_con}, position{pos}, fname_stat, col, [tmp_data.(fname)], stim_con);
                                    perform_ttests = 1;

                                else
                                    warning(sprintf('The distractor colors do not correspond to the easy and the difficult distractor. \nIn 3 way repeated measures ANOVA: double stimulus trials - difficulty - stimulation \nANOVA is not performed!'));
                                    perform_ttests = 0;
                                end
                            end 
                        end
                    end
                    %% t-test & FDR correction ... added to the structure as 5 to 9 colums
                    if perform_ttests % if all distractor colors are the correct ones
                        for dv = 1:numel(dep_var) % number of choices
                            if c_l_e == 1 %%% CHOICES %%%
                                fname_stat = [dep_var{dv},'_choice']; % fieldname for stat_data structure
                                fname_tmp = 'true_prop_choice_per_session_';
                                fname = [fname_tmp,dep_var{dv}]; % fieldname for tmp_data structure
                            elseif c_l_e == 2 %%% LATENCY %%%
                                fname_stat = [dep_var{dv},'_latency']; % fieldname for stat_data structure
                                fname = ['mean_latency_',dep_var{dv},'_choice_per_session'];
                            elseif c_l_e == 3 %%% ERRORS %%%
                                fname_stat = [dep_var{dv},'_errors']; % fieldname for stat_data structure
                                fname = ['error_rate_',dep_var{dv},'_per_session']; % fieldname for tmp_data structure
                            end
                            for col = 1:size(stat_struct_ttest.(tr_con2stat{tr_con}).(position{pos}).(fname_stat),1); % for all colors (rows in data_struct) in this trial condition
                                for comp = 1:size(stat_struct_ttest.(tr_con2stat{tr_con}).(position{pos}).(fname_stat){col,1},1)
                                    x = stat_struct_ttest.(tr_con2stat{tr_con}).(position{pos}).(fname_stat){col,1}{comp,1};
                                    y = stat_struct_ttest.(tr_con2stat{tr_con}).(position{pos}).(fname_stat){col,1}{comp,2};
                                    
                                    [H,P,~,stats] = ttest(x,y,'tail','both','alpha',Settings.alpha); %ttest KK
                                    H(isnan(H)) = 0; % if x & y are vectors of zeros or if there are not enough data H will be NaN which cannot be converted into logicals in the plot_bargraph function, so this will be set to 0, i.e. not significant; P is still NaN
                                    
                                    stat_struct_ttest.(tr_con2stat{tr_con}).(position{pos}).(fname_stat){col,1}{comp,5} = H; % is this comparison significantly different?
                                    stat_struct_ttest.(tr_con2stat{tr_con}).(position{pos}).(fname_stat){col,1}{comp,6} = P; % p-value of ttest
                                    stat_struct_ttest.(tr_con2stat{tr_con}).(position{pos}).(fname_stat){col,1}{comp,7} = stats.tstat; % t-value of ttest
                                    stat_struct_ttest.(tr_con2stat{tr_con}).(position{pos}).(fname_stat){col,1}{comp,8} = stats.df; % degrees of freedom of ttest
                                end
                                if Settings.FDR
                                    PValues = [stat_struct_ttest.(tr_con2stat{tr_con}).(position{pos}).(fname_stat){col,1}{:,6}];
                                    
                                    [H_adj, crit_p, adj_ci_cvrg, P_adj]=fdr_bh(PValues,Settings.alpha,'pdep','no');
                                    H_adj(isnan(H_adj)) = 0; % if x & y are vectors of zeros or if there are not enough data H will be NaN which cannot be converted into logicals in the plot_bargraph function, so this will be set to 0, i.e. not significant; P is still NaN
                                    
                                    stat_struct_ttest.(tr_con2stat{tr_con}).(position{pos}).(fname_stat){col,1}(:,9) = num2cell(H_adj'); % corrected indicators of significance
                                    stat_struct_ttest.(tr_con2stat{tr_con}).(position{pos}).(fname_stat){col,1}(:,10) = num2cell(P_adj'); % corrected p-values of ttest
                                    
                                    %%% Check if FDR corrected results are different from uncorrected %%%
                                    if ~isequal(num2cell([stat_struct_ttest.(tr_con2stat{tr_con}).(position{pos}).(fname_stat){col,1}{:,5}]'),...
                                            stat_struct_ttest.(tr_con2stat{tr_con}).(position{pos}).(fname_stat){col,1}(:,9))
                                        disp(['stat_struct_ttest',', ',tr_con2stat{tr_con},', ',position{pos},', ',fname_stat,', ','color: ',num2str(col)]);
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    if Settings.FDR
        stat_struct_ttest.legend = {fname,['compared to ',fname],'Stimulation condition','compared to stimulation condition','Significant?','P-value','t','df','FDR Significant?','FDR adjusted p'};
    else
        stat_struct_ttest.legend = {fname,['compared to ',fname],'Stimulation condition','compared to stimulation condition','Significant?','P-value','t','df'};
    end
end
%%
%%%% THIS IS THE STUPID ANOVA DESIGN %%%%
%%%  3-WAY REAPEATED MEASURES ANOVA OF DOUBLE STIMULUS TRIALS (except targ-targ) WITH DIFFICULTY AND STIMULATION CONDITIONS AS 2ND AND 3RD FACTORS %%%
if Settings.choice_rmanova3way_DoubleStimDisplay_Difficulty_Stimulation
    stat_data = []; stat_struct_ttest_anova = [];
    assignment = [];
    perform_ANOVA = 1; anova_con = 0;
    choices_stat = {'contra','fixation','ipsi'};
    tr_con2stat = {'targ_distr','targ_distr','distr_distr'};
    num_tr_con2stat = numel(tr_con2stat);
    distr_G_R_ratio = [0.1796875;1];
    if Settings.FDR
        mode = 'compare_to_control';
    else
        mode = 'all_comparisons';
    end
    
    for tr_con = 1:num_tr_con2stat
        position = {'contra','ipsi','combined'}; % refers to tr_con2stat
        for pos = tr_con; % for the respective position in this trial condition
            for stim = 1:num_stim_con
                for col = 1:size(data_struct_per_session.(tr_con2stat{tr_con}).(position{pos}).(stim_con{stim}),1); % for all colors (rows in data_struct) in this trial condition
                    anova_con = anova_con + 1;
                    tmp_data = [data_struct_per_session.(tr_con2stat{tr_con}).(position{pos}).(stim_con{stim}){col,:}]; % tmp_data with all sessions concatenated
                    if isequal(unique([tmp_data.G_R_ratios]),distr_G_R_ratio) % check if distractor colors are the easy and difficult distractor in all sessions, col = 1 is always difficult distr, col = 2 always easy distr
                        for ch = 1:numel(choices_stat) % number choice options
                            fname_stat = [choices_stat{ch},'_choices']; % fieldname for stat_data structure
                            fname = ['true_prop_choice_per_session_',choices_stat{ch}]; % fieldname for tmp_data structure
                            stat_data.(fname_stat){1,anova_con} = [tmp_data.(fname)];
                            assignment(anova_con,:) = [tr_con,stim,col];
                            
                            %%% POST HOC T-TEST %%%
                            % Test the null hypothesis that the pairwise difference between data vectors x and y has a mean equal to zero.
                            if Settings.posthoc_ttest
%                                 stat_struct_ttest = get_stat_struct_ttest (mode, stat_struct_ttest, stim, tr_con2stat{tr_con}, position{pos}, fname_stat, col, [tmp_data.(fname)], stim_con);
                                stat_struct_ttest_anova = get_stat_struct_ttest (mode, stat_struct_ttest_anova, stim, tr_con2stat{tr_con}, position{pos}, fname_stat, col, [tmp_data.(fname)], stim_con);
                                
                            end
                        end
                    else
                        warning(sprintf('The distractor colors do not correspond to the easy and the difficult distractor. \nIn 3 way repeated measures ANOVA: double stimulus trials - difficulty - stimulation \nANOVA is not performed!'));
                        perform_ANOVA = 0;
                    end
                end
            end
            if Settings.posthoc_ttest && perform_ANOVA % if all distractor colors are the correct ones
                for ch = 1:3 % number of choices
                    fname_stat = [choices_stat{ch},'_choices']; % fieldname for stat_data structure
                    for col = 1:size(stat_struct_ttest_anova.choice_rmanova3way_DoubleStimDisplay_Difficulty_Stimulation.(tr_con2stat{tr_con}).(position{pos}).(fname_stat),1); % for all colors (rows in data_struct) in this trial condition
                        for comp = 1:size(stat_struct_ttest_anova.choice_rmanova3way_DoubleStimDisplay_Difficulty_Stimulation.(tr_con2stat{tr_con}).(position{pos}).(fname_stat){col,1},1)
                            x = stat_struct_ttest_anova.choice_rmanova3way_DoubleStimDisplay_Difficulty_Stimulation.(tr_con2stat{tr_con}).(position{pos}).(fname_stat){col,1}{comp,1};
                            y = stat_struct_ttest_anova.choice_rmanova3way_DoubleStimDisplay_Difficulty_Stimulation.(tr_con2stat{tr_con}).(position{pos}).(fname_stat){col,1}{comp,2};
                            [H,P,~,stats] = ttest(x,y,'tail','both','alpha',Settings.alpha);
                            stat_struct_ttest_anova.choice_rmanova3way_DoubleStimDisplay_Difficulty_Stimulation.(tr_con2stat{tr_con}).(position{pos}).(fname_stat){col,1}{comp,5} = H; % is this comparison significantly different?
                            stat_struct_ttest_anova.choice_rmanova3way_DoubleStimDisplay_Difficulty_Stimulation.(tr_con2stat{tr_con}).(position{pos}).(fname_stat){col,1}{comp,6} = P; % p-value of ttest
                            stat_struct_ttest_anova.choice_rmanova3way_DoubleStimDisplay_Difficulty_Stimulation.(tr_con2stat{tr_con}).(position{pos}).(fname_stat){col,1}{comp,7} = stats.tstat; % p-value of ttest
                            stat_struct_ttest_anova.choice_rmanova3way_DoubleStimDisplay_Difficulty_Stimulation.(tr_con2stat{tr_con}).(position{pos}).(fname_stat){col,1}{comp,8} = stats.df; % p-value of ttest
                        end
                    end
                end
            end
        end
    end
    if perform_ANOVA
        varnames = {'Stimulus' 'Stimulation' 'Difficulty'};
        levnames = {{'Contra-Targ Ipsi-Distr' 'Contra-Distr Ipsi-Targ' 'Contra-Distr-Ipsi-Distr'} {'None' num2str(-80) 'Go' num2str(+80)} {'Difficult' 'Easy'}};
        %         assignment = [1 1 1; 1 1 2; 1 2 1; 1 2 2; 1 3 1; 1 3 2; 1 4 1; 1 4 2;...
        %             2 1 1; 2 1 2; 2 2 1; 2 2 2; 2 3 1; 2 3 2; 2 4 1; 2 4 2;...
        %             3 1 1; 3 1 2; 3 2 1; 3 2 2; 3 3 1; 3 3 2; 3 4 1; 3 4 2];
        dep_var = fieldnames(stat_data);
        for dv = 1:numel(fieldnames(stat_data))
            stat_results_anova.choice_double_stimuli.(dep_var{dv}) = threeway_rmanova(stat_data.(dep_var{dv}), assignment, varnames, levnames, Settings.alpha);
        end
    end
end
%%%% THIS WAS THE STUPID ANOVA DESIGN %%%%
%%

%-----------------------------------------------------------------
    %%% CHOICES AND ERROR RATES (ANOVA AND POSTHOC TTESTS) %%%
%-----------------------------------------------------------------

%%%  2-WAY REAPEATED MEASURES ANOVA OF DISTR-DISTR (EASY AND DIFFICULT) & TARG-TARG TRIALS WITH STIMULATION CONDITIONS AS 2ND FACTOR %%%
%%% CHOICES OR ERROR RATES %%%
choice_or_errors = [Settings.choice_rmanova2way_DistrDistrTargTarg_Stimulation, Settings.error_rate_rmanova2way_DistrDistrTargTarg_Stimulation];
stat_struct_ttest_anova = [];
for rep = 1:numel(choice_or_errors)
    if choice_or_errors(rep) % if the current ANOVA of choices / error rates is to be calculated     
        stat_data = []; 
        assignment = [];
        perform_ANOVA = 1; anova_con = 0;
        if rep == 1 % choices
            dep_var = {'left','fixation','right'};
        elseif rep == 2 % error rates
            dep_var = {'fixation_break','target_acquis','target_hold'};
        end
        tr_con2stat = {'targ_targ','distr_distr','distr_distr'}; %%%%%%%%%%%%%%%% easy and difficult distractor combined
        num_tr_con2stat = numel(tr_con2stat);
        difficulty = [1,1,2];
        
        distr_G_R_ratio = [0.1796875;1];
        distr_G_R_ratio = [0.0859;1]; %Cornelius

        
        for tr_con = 1:num_tr_con2stat
            position = {'combined','combined','combined'}; % refers to tr_con2stat
            for pos = tr_con; % for the respective position in this trial condition
                for stim = 1:num_stim_con
                    for col = tr_con; % for all colors (rows in data_struct) in this trial condition
                        anova_con = anova_con + 1;
                        tmp_data = [data_struct_per_session.(tr_con2stat{tr_con}).(position{pos}).(stim_con{stim}){difficulty(col),:}]; % tmp_data with all sessions concatenated
                        if isequal(unique([tmp_data.G_R_ratios]),distr_G_R_ratio) || isequal(unique([tmp_data.G_R_ratios]),0) % check if distractor colors are the easy and difficult distractor in all sessions, col = 1 is always difficult distr, col = 2 always easy distr
                            for dv = 1:numel(dep_var) % number choice options
                                
                                if rep == 1 % choices
                                    fname_stat = [dep_var{dv},'_choices']; % fieldname for stat_data structure
                                    fname = ['true_prop_choice_per_session_',dep_var{dv}]; % fieldname for tmp_data structure
                                elseif rep == 2 % error rates
                                    fname_stat = [dep_var{dv},'_errors']; % fieldname for stat_data structure
                                    fname = ['error_rate_',dep_var{dv},'_per_session']; % fieldname for tmp_data structure
                                end
                                
                                stat_data.(fname_stat).values(:,anova_con) = [tmp_data.(fname)]; % choice proportion
                                stat_data.(fname_stat).sessions(:,anova_con) = [1:1:numel(tmp_data)]'; % subjects, i.e. sessions
                                stat_data.(fname_stat).tr_con(:,anova_con) = tr_con*ones(numel(tmp_data),1); % factor 1, i.e. trial stimulus
                                stat_data.(fname_stat).stimulation(:,anova_con) = stim*ones(numel(tmp_data),1); % factor 2, i.e. stimulation
                                
                            end
                        else
                            warning(sprintf('The distractor colors do not correspond to the easy and the difficult distractor. \nIn 3 way repeated measures ANOVA: double stimulus trials - difficulty - stimulation \nANOVA is not performed!'));
                            perform_ANOVA = 0;
                        end
                    end
                end
            end
        end
        %%% POST HOC T-TEST %%%
        if perform_ANOVA && Settings.posthoc_ttest && rep == 1 % choices
            
            %%% GENERAL CHOICE BIAS - T-TESTS ON SINGLE STIMULUS DISPLAYS %%%
            difficulty = {'none','difficult','easy'};
            for n = 1:numel(difficulty) % easy & difficult
                x = stat_data.contra_choices.values(:,4*(n-1)+1); % targ-targ, contra choices, easy for n = 1; difficult distr-distr, contra choices for n = 2; easy distr-distr, contra choices for n = 3;
                y = stat_data.ipsi_choices.values(:,4*(n-1)+1); % ipsi target, ipsi choices, easy for ne = 1; ipsi target, ipsi choices, difficult for n = 2;
                [H,P,~,stats] = ttest(x,y,'tail','both','alpha',Settings.alpha);
                
                stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.general_bias{n,1} = difficulty{n};
                stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.general_bias{n,2} = tr_con2stat{n};
                stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.general_bias{n,3} = mean(x);
                stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.general_bias{n,4} = nanstd(x)/sqrt(size(x,1)); % SEM
                stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.general_bias{n,5} = 'contra_choice';
                stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.general_bias{n,6} = mean(y);
                stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.general_bias{n,7} = nanstd(y)/sqrt(size(y,1)); % SEM
                stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.general_bias{n,8} = 'ipsi_choice';
                stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.general_bias{n,9} = H; % is this comparison significantly different?
                stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.general_bias{n,10} = P; % p-value of ttest
                stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.general_bias{n,11} = stats.tstat; % t-value of ttest
                stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.general_bias{n,12} = stats.df; % degrees of freedom
                
            end
            %%% FDR CORRECTION %%%
            if Settings.FDR
                PValues = [stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.general_bias{:,10}];
                
                [H_adj, ~, ~, P_adj]=fdr_bh(PValues,Settings.alpha,'pdep','no');
                H_adj(isnan(H_adj)) = 0; % if x & y are vectors of zeros or if there are not enough data H will be NaN which cannot be converted into logicals in the plot_bargraph function, so this will be set to 0, i.e. not significant; P is still NaN
                
                stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.general_bias(:,13) = num2cell(H_adj'); % corrected indicators of significance
                stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.general_bias(:,14) = num2cell(P_adj'); % corrected p-values of ttest
                
                stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.general_bias = ... % header if FDR
                    [{'difficulty','tr_con','mean','SEM','choice','mean','SEM','compared to choice','significant','p','t','df','signif_FDR','p_FDR'};...
                    stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.general_bias];
                
            else
                stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.general_bias = ... % header if no FDR
                    [{'difficulty','tr_con','mean','SEM','choice','mean','SEM','compared to choice','significant','p','t','df'};...
                    stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.general_bias];
            end
        end
        
        if perform_ANOVA
            dep_var = fieldnames(stat_data);
            for dv = 1:numel(dep_var)
                
                Y.(dep_var{dv}) = reshape(stat_data.(dep_var{dv}).values,pos*stim*numel(tmp_data),1); % concatenate to one column vector
                S = reshape(stat_data.(dep_var{dv}).sessions,pos*stim*numel(tmp_data),1); % subjects, i.e. sessions
                F1 = reshape(stat_data.(dep_var{dv}).tr_con,pos*stim*numel(tmp_data),1); % factor 1, i.e. tr_con
                F2 = reshape(stat_data.(dep_var{dv}).stimulation,pos*stim*numel(tmp_data),1); % factor 2, i.e. stimulation
                varnames = {'Stimulus' 'Stimulation'};
                if rep == 1 % choices
                    stat_results_anova.choice_DistrDistr_Targ_Targ.(dep_var{dv}) = rm_anova2(Y.(dep_var{dv}),S,F1,F2,varnames);
                    stat_results_anova.choice_DistrDistr_Targ_Targ.([dep_var{dv},'_stat_data']) = stat_data.(dep_var{dv});
                elseif rep == 2 % error rates
                    stat_results_anova.errors_DistrDistr_Targ_Targ.(dep_var{dv}) = rm_anova2(Y.(dep_var{dv}),S,F1,F2,varnames);
                    stat_results_anova.errors_DistrDistr_Targ_Targ.([dep_var{dv},'_stat_data']) = stat_data.(dep_var{dv});
                end
                
            end
        end
    end
end

%%%  3-WAY REAPEATED MEASURES ANOVA OF TARG-DISTR TRIALS (Left target vs. right target) WITH DIFFICULTY AND STIMULATION CONDITIONS AS 2ND AND 3RD FACTORS %%%
%%% CHOICES OR ERROR RATES %%%
choice_or_errors = [Settings.choice_rmanova3way_TargDistr_Difficulty_Stimulation,Settings.error_rate_rmanova3way_TargDistr_Difficulty_Stimulation];
stat_struct_ttest_anova = []; 
for rep = 1:numel(choice_or_errors)
    if choice_or_errors(rep) % if the current ANOVA of choices / error rates is to be calculated        
        stat_data = []; ttest_data = [];
        assignment = [];
        perform_ANOVA = 1; anova_con = 0;
        if rep == 1 % choices
            dep_var = {'contra','fixation','ipsi'};
        elseif rep == 2 % error rates
            dep_var = {'fixation_break','target_acquis','target_hold'};
        end
        tr_con2stat = {'targ_distr','targ_distr'};
        num_tr_con2stat = numel(tr_con2stat);
        distr_G_R_ratio = [0.1796875;1];
        for tr_con = 1:num_tr_con2stat
            position = {'contra','ipsi'}; % refers to tr_con2stat
            for pos = tr_con; % for the respective position in this trial condition
                for stim = 1:num_stim_con
                    for col = 1:size(data_struct_per_session.(tr_con2stat{tr_con}).(position{pos}).(stim_con{stim}),1); % for all colors (rows in data_struct) in this trial condition
                        anova_con = anova_con + 1;
                        tmp_data = [data_struct_per_session.(tr_con2stat{tr_con}).(position{pos}).(stim_con{stim}){col,:}]; % tmp_data with all sessions concatenated
                        if isequal(unique([tmp_data.G_R_ratios]),distr_G_R_ratio) % check if distractor colors are the easy and difficult distractor in all sessions, col = 1 is always difficult distr, col = 2 always easy distr
                            for dv = 1:numel(dep_var) % number choice / error rate options
                                
                                if rep == 1 % choices
                                    fname_stat = [dep_var{dv},'_choices']; % fieldname for stat_data structure
                                    fname = ['true_prop_choice_per_session_',dep_var{dv}]; % fieldname for tmp_data structure
                                elseif rep == 2 % error rates
                                    fname_stat = [dep_var{dv},'_errors']; % fieldname for stat_data structure
                                    fname = ['error_rate_',dep_var{dv},'_per_session']; % fieldname for tmp_data structure
                                end
                                
                                stat_data.(fname_stat){1,anova_con} = [tmp_data.(fname)];
                                assignment(anova_con,:) = [pos,stim,col];
                                % mode = 'compare_to_control';
                                % stat_struct_ttest_anova = get_stat_struct_ttest (mode, stat_struct_ttest_anova, stim, tr_con2stat{tr_con}, position{pos}, fname_stat, col, [tmp_data.(fname)], stim_con);
                            end
                        else
                            warning(sprintf('The distractor colors do not correspond to the easy and the difficult distractor. \nIn 3 way repeated measures ANOVA: double stimulus trials - difficulty - stimulation \nANOVA is not performed!'));
                            perform_ANOVA = 0;
                        end
                    end
                end
            end
        end
        
        %%% POST HOC T-TEST %%%
        if perform_ANOVA && Settings.posthoc_ttest && rep == 1 % choices
            
            %%% GENERAL CHOICE BIAS - T-TESTS ON SINGLE STIMULUS DISPLAYS %%%
            difficulty = {'difficult','easy'};
            for n = 1:numel(difficulty) % easy & difficult
                x = stat_data.contra_choices{1,0+n}; % contra target, contra choices, easy for n = 1; contra target, contra choices, difficult for n = 2; indcs: 1 and 2 of assignment
                y = stat_data.ipsi_choices{1,8+n}; % ipsi target, ipsi choices, easy for ne = 1; ipsi target, ipsi choices, difficult for n = 2; indcs: 9 and 10
                [H,P,~,stats] = ttest(x,y,'tail','both','alpha',Settings.alpha);
                
                stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.general_bias{n,1} = difficulty{n};
                stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.general_bias{n,2} = mean(x);
                stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.general_bias{n,3} = nanstd(x)/sqrt(size(x,1)); % SEM
                stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.general_bias{n,4} = position{1};
                stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.general_bias{n,5} = mean(y);
                stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.general_bias{n,6} = nanstd(y)/sqrt(size(y,1)); % SEM
                stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.general_bias{n,7} = position{2};
                stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.general_bias{n,8} = H; % is this comparison significantly different?
                stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.general_bias{n,9} = P; % p-value of ttest
                stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.general_bias{n,10} = stats.tstat; % t-value of ttest
                stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.general_bias{n,11} = stats.df; % degrees of freedom
                
            end
            %%% FDR CORRECTION %%%
            if Settings.FDR
                PValues = [stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.general_bias{:,9}];
                
                [H_adj, ~, ~, P_adj]=fdr_bh(PValues,Settings.alpha,'pdep','no');
                H_adj(isnan(H_adj)) = 0; % if x & y are vectors of zeros or if there are not enough data H will be NaN which cannot be converted into logicals in the plot_bargraph function, so this will be set to 0, i.e. not significant; P is still NaN
                
                stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.general_bias(:,12) = num2cell(H_adj'); % corrected indicators of significance
                stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.general_bias(:,13) = num2cell(P_adj'); % corrected p-values of ttest
                
                stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.general_bias = ... % header if FDR
                    [{'difficulty','mean','SEM','targ position','mean','SEM','compared to targ pos','significant','p','t','df','signif_FDR','p_FDR'};...
                    stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.general_bias];
                
            else
                stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.general_bias = ... % header if no FDR
                    [{'difficulty','mean','SEM','targ position','mean','SEM','compared to targ pos','significant','p','t','df'};...
                    stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.general_bias];
            end
            
            %%% STIMULUS X STIMULATION INTERACTION -> TAKE MEAN ACROSS DIFFICULTY %%%
            mean_idcs = {{1,2},{3,4},{5,6},{7,8},{9,10},{11,12},{13,14},{15,16}}; % which rows should be taken together and averaged across
            comparisons = {{1,2},{1,3},{1,4},{5,6},{5,7},{5,8}};
            stim_con4ttest = [stim_con;stim_con];
            position4ttest = {'contra','contra','contra','ipsi','ipsi','ipsi'};
            choices_stat = {'contra','fixation','ipsi'};
            for ch = 1:numel(choices_stat)
                fname_stat = [choices_stat{ch},'_choices']; % fieldname for stat_data structure
                for i = 1:numel(mean_idcs)
                    for ses = 1:length(stat_data.(fname_stat){1,mean_idcs{i}{1}})
                        ttest_data.(fname_stat){i}(1,ses) = mean([stat_data.(fname_stat){mean_idcs{i}{1}}(ses),stat_data.(fname_stat){mean_idcs{i}{2}}(ses)]);
                        
                    end
                end
                for comp = 1:numel(comparisons) % which elements should be compared in a t test?
                    x = ttest_data.(fname_stat){comparisons{comp}{1}};
                    y = ttest_data.(fname_stat){comparisons{comp}{2}};
                    [H,P,~,stats] = ttest(x,y,'tail','both','alpha',Settings.alpha);
                    stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.stimulus_x_stimulation.(fname_stat){comp,1} = position4ttest{comp};
                    stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.stimulus_x_stimulation.(fname_stat){comp,2} = stim_con4ttest{comparisons{comp}{1}};
                    stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.stimulus_x_stimulation.(fname_stat){comp,3} = stim_con4ttest{comparisons{comp}{2}};
                    stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.stimulus_x_stimulation.(fname_stat){comp,4} = H; % is this comparison significantly different?
                    stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.stimulus_x_stimulation.(fname_stat){comp,5} = P; % p-value of ttest
                    stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.stimulus_x_stimulation.(fname_stat){comp,6} = stats.tstat; % t-value of ttest
                    stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.stimulus_x_stimulation.(fname_stat){comp,7} = stats.df; % degrees of freedom
                end
            end
            %%% FDR CORRECTION %%%
            if Settings.FDR
                PValues = [stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.stimulus_x_stimulation.(fname_stat){:,5}];
                
                [H_adj, ~, ~, P_adj]=fdr_bh(PValues,Settings.alpha,'pdep','no');
                H_adj(isnan(H_adj)) = 0; % if x & y are vectors of zeros or if there are not enough data H will be NaN which cannot be converted into logicals in the plot_bargraph function, so this will be set to 0, i.e. not significant; P is still NaN
                
                stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.stimulus_x_stimulation.(fname_stat)(:,8) = num2cell(H_adj'); % corrected indicators of significance
                stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.stimulus_x_stimulation.(fname_stat)(:,9) = num2cell(P_adj'); % corrected p-values of ttest
                
                stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.stimulus_x_stimulation.(fname_stat) = ... % header if FDR
                    [{'target position','stimulation','compared to stimulation','significant','p','t','df','signif_FDR','p_FDR'};...
                    stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.stimulus_x_stimulation.(fname_stat)];
                
            else
                stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.stimulus_x_stimulation.(fname_stat) = ... % header if no FDR
                    [{'target position','stimulation','compared to stimulation','significant','p','t','df'};...
                    stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.stimulus_x_stimulation.(fname_stat)];
            end
            
            %%% DIFFICULTY X STIMULATION INTERACTION -> TAKE MEAN ACROSS STIMULUS %%% 
            assignment_averaged = assignment(1:8,:);
            mean_idcs = {{1,9},{2,10},{3,11},{4,12},{5,13},{6,14},{7,15},{8,16}}; % which rows in assignment should be taken together and averaged across
            comparisons = {{1,3},{1,5},{1,7},{2,4},{2,6},{2,8}}; % which of the mean_idcs are to be compared against each other
            stim_con4ttest = [stim_con';stim_con'];
            difficulty4ttest = {'difficult','difficult','difficult','easy','easy','easy'};
            for ch = 1:numel(choices_stat)
                fname_stat = [choices_stat{ch},'_choices']; % fieldname for stat_data structure
                for i = 1:numel(mean_idcs)
                    for ses = 1:length(stat_data.(fname_stat){1,mean_idcs{i}{1}})
                        ttest_data.(fname_stat){i}(1,ses) = mean([stat_data.(fname_stat){mean_idcs{i}{1}}(ses),stat_data.(fname_stat){mean_idcs{i}{2}}(ses)]);
                        
                    end
                end
                for comp = 1:numel(comparisons) % which elements should be compared in a t test?
                    x = ttest_data.(fname_stat){comparisons{comp}{1}};
                    y = ttest_data.(fname_stat){comparisons{comp}{2}};
                    [H,P,~,stats] = ttest(x,y,'tail','both','alpha',Settings.alpha);
                    stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.difficulty_x_stimulation.(fname_stat){comp,1} = difficulty4ttest{comp};
                    stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.difficulty_x_stimulation.(fname_stat){comp,2} = stim_con4ttest{comparisons{comp}{1}};
                    stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.difficulty_x_stimulation.(fname_stat){comp,3} = stim_con4ttest{comparisons{comp}{2}};
                    stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.difficulty_x_stimulation.(fname_stat){comp,4} = H; % is this comparison significantly different?
                    stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.difficulty_x_stimulation.(fname_stat){comp,5} = P; % p-value of ttest
                    stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.difficulty_x_stimulation.(fname_stat){comp,6} = stats.tstat; % t-value of ttest
                    stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.difficulty_x_stimulation.(fname_stat){comp,7} = stats.df; % degrees of freedom
                end
            end
            %%% FDR CORRECTION %%%
            if Settings.FDR
                PValues = [stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.difficulty_x_stimulation.(fname_stat){:,5}];
                
                [H_adj, ~, ~, P_adj]=fdr_bh(PValues,Settings.alpha,'pdep','no');
                H_adj(isnan(H_adj)) = 0; % if x & y are vectors of zeros or if there are not enough data H will be NaN which cannot be converted into logicals in the plot_bargraph function, so this will be set to 0, i.e. not significant; P is still NaN
                
                stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.difficulty_x_stimulation.(fname_stat)(:,8) = num2cell(H_adj'); % corrected indicators of significance
                stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.difficulty_x_stimulation.(fname_stat)(:,9) = num2cell(P_adj'); % corrected p-values of ttest
                
                stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.difficulty_x_stimulation.(fname_stat) = ... % header if FDR
                    [{'difficulty','stimulation','compared to stimulation','significant','p','t','df','signif_FDR','p_FDR'};...
                    stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.difficulty_x_stimulation.(fname_stat)];
                
            else
                stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.difficulty_x_stimulation.(fname_stat) = ... % header if no FDR
                    [{'difficulty','stimulation','compared to stimulation','significant','p','t','df'};...
                    stat_struct_ttest_anova.choice_rmanova3way_TargDistr_Difficulty_Stimulation.difficulty_x_stimulation.(fname_stat)];
            end
        end
        
        %%% ANOVA %%%
        if perform_ANOVA
            varnames = {'Stimulus' 'Stimulation' 'Difficulty'};
            levnames = {{'Contra-Targ Ipsi-Distr' 'Contra-Distr Ipsi-Targ' 'Contra-Distr-Ipsi-Distr'} {'None' num2str(-80) 'Go' num2str(+80)} {'Difficult' 'Easy'}};
            %         assignment = [1 1 1; 1 1 2; 1 2 1; 1 2 2; 1 3 1; 1 3 2; 1 4 1; 1 4 2;...
            %             2 1 1; 2 1 2; 2 2 1; 2 2 2; 2 3 1; 2 3 2; 2 4 1; 2 4 2;...
            %             3 1 1; 3 1 2; 3 2 1; 3 2 2; 3 3 1; 3 3 2; 3 4 1; 3 4 2];
            dep_var = fieldnames(stat_data);
            for dv = 1:numel(fieldnames(stat_data))
                if rep == 1 % choices
                    stat_results_anova.choice_TargDistr.(dep_var{dv}) = threeway_rmanova(stat_data.(dep_var{dv}), assignment, varnames, levnames, Settings.alpha);
                    stat_results_anova.choice_TargDistr.([dep_var{dv},'_stat_data']) = stat_data.(dep_var{dv});
                    stat_results_anova.choice_TargDistr.([dep_var{dv},'_stat_data_assignment']) = assignment;
                elseif rep == 2 % error rates
                    stat_results_anova.errors_TargDistr.(dep_var{dv}) = threeway_rmanova(stat_data.(dep_var{dv}), assignment, varnames, levnames, Settings.alpha);
                    stat_results_anova.errors_TargDistr.([dep_var{dv},'_stat_data']) = stat_data.(dep_var{dv});
                    stat_results_anova.errors_TargDistr.([dep_var{dv},'_stat_data_assignment']) = assignment;
                end
            end
        end
    end
end

%%%  3-WAY REAPEATED MEASURES ANOVA OF SINGLE STIMULUS TRIALS WITH SIDE AND STIMULATION CONDITIONS AS 2ND AND 3RD FACTORS %%%
%%% CHOICES OR ERROR RATES %%%
choice_or_errors = [Settings.choice_rmanova3way_SingleStimDisplay_Side_Stimulation,Settings.error_rate_rmanova3way_SingleStimDisplay_Side_Stimulation];
for rep = 1:numel(choice_or_errors)
    if choice_or_errors(rep) % if the current ANOVA of choices / error rates is to be calculated
        stat_data = []; assignment = []; ttest_data = [];
        perform_ANOVA = 1;  anova_con = 0;
        tr_con2stat = {'single_targ','single_distr','single_distr'};
        num_tr_con2stat = numel(tr_con2stat);
        difficulty = [1,1,2]; % first row in single_targ trials (no distr color), 1st row in single_distr (diff color), 2nd row in single_distr trials (easy color)
        position = {'contra','ipsi'}; % refers to tr_con2stat
        side_levels = [];
        dep_var = {'fixation_break','target_acquis','target_hold'};
        for tr_con = 1:num_tr_con2stat
            for pos = 1:numel(position); % for all positions in this trial condition
                for stim = 1:num_stim_con
                    for col = difficulty(tr_con) % for the respective color (rows in data_struct) per trial condition
                        anova_con = anova_con + 1; % data for different conditions are concatenated vertically
                        tmp_data = [data_struct_per_session.(tr_con2stat{tr_con}).(position{pos}).(stim_con{stim}){col,:}]; % tmp_data with all sessions concatenated
                        if rep == 1 % choices
                            if isequal(tr_con2stat{tr_con},'single_targ') && isequal(unique([tmp_data.G_R_ratios]), 0) % if it's really single_targ trials without distractors
                                stat_data.fixation{1,anova_con} = [tmp_data.true_prop_choice_per_session_fixation]; % fixation choice
                                stat_data.saccade{1,anova_con} = [tmp_data.hitrate_per_session]; % saccade choice
                            elseif isequal(tr_con2stat{tr_con},'single_distr') && isequal(unique([tmp_data.G_R_ratios]), [0.1796875;1]) % if it's really single_distr trials with the correct distractor colors
                                stat_data.fixation{1,anova_con} = [tmp_data.hitrate_per_session]; % fixation choice
                                stat_data.saccade{1,anova_con} = [tmp_data.(['true_prop_choice_per_session_',position{pos}])]; % saccade choice
                            else
                                warning(sprintf('The distractor colors do not correspond to the easy and the difficult distractor. \nIn 3 way repeated measures ANOVA: single stimulus trials - side - stimulation \nANOVA is not performed!'));
                                perform_ANOVA = 0;
                            end
                        elseif rep == 2 % error rates
                            for dv = 1:numel(dep_var) % number choice / error rate options
                                fname_stat = [dep_var{dv},'_errors']; % fieldname for stat_data structure
                                fname = ['error_rate_',dep_var{dv},'_per_session']; % fieldname for tmp_data structure
                                stat_data.(fname_stat){1,anova_con} = [tmp_data.(fname)];
                            end
                        end
                        assignment(anova_con,:) = [tr_con,pos,stim];
                    end
                end
                side_levels = [side_levels,position(pos)];
            end
        end
        
        %%% POST HOC T-TEST %%%
        if perform_ANOVA && Settings.posthoc_ttest && rep == 1 % choices
            %%% STIMULUS X STIMULATION INTERACTION -> TAKE MEAN ACROSS SIDE %%%
            assignment_averaged = assignment([1:4,9:12,17:20],:);
            mean_idcs = {{1,5},{2,6},{3,7},{4,8},{9,13},{10,14},{11,15},{12,16},{17,21},{18,22},{19,23},{20,24}}; % which rows should be taken together and averaged across
            comparisons = {{1,2},{1,3},{1,4},{5,6},{5,7},{5,8},{9,10},{9,11},{9,12}};
            stim_con4ttest = [stim_con;stim_con;stim_con];
            tr_con4ttest = {'single_targ','single_targ','single_targ','single_diff_distr','single_diff_distr','single_diff_distr','single_easy_distr','single_easy_distr','single_easy_distr'};
            choices_stat = {'fixation','saccade'};
            for ch = 1:numel(choices_stat)
                fname_stat = choices_stat{ch}; % fieldname for stat_data structure
                for i = 1:numel(mean_idcs)
                    for ses = 1:length(stat_data.(fname_stat){1,mean_idcs{i}{1}})
                        ttest_data.(fname_stat){i}(1,ses) = mean([stat_data.(fname_stat){mean_idcs{i}{1}}(ses),stat_data.(fname_stat){mean_idcs{i}{2}}(ses)]);
                        
                    end
                end
                for comp = 1:numel(comparisons) % which elements should be compared in a t test?
                    x = ttest_data.(fname_stat){comparisons{comp}{1}};
                    y = ttest_data.(fname_stat){comparisons{comp}{2}};
                    [H,P,~,stats] = ttest(x,y,'tail','both','alpha',Settings.alpha);
                    stat_struct_ttest_anova.choice_rmanova3way_SingleStimDisplay_Side_Stimulation.stimulus_x_stimulation.(fname_stat){comp,1} = tr_con4ttest{comp};
                    stat_struct_ttest_anova.choice_rmanova3way_SingleStimDisplay_Side_Stimulation.stimulus_x_stimulation.(fname_stat){comp,2} = stim_con4ttest{comparisons{comp}{1}};
                    stat_struct_ttest_anova.choice_rmanova3way_SingleStimDisplay_Side_Stimulation.stimulus_x_stimulation.(fname_stat){comp,3} = stim_con4ttest{comparisons{comp}{2}};
                    stat_struct_ttest_anova.choice_rmanova3way_SingleStimDisplay_Side_Stimulation.stimulus_x_stimulation.(fname_stat){comp,4} = H; % is this comparison significantly different?
                    stat_struct_ttest_anova.choice_rmanova3way_SingleStimDisplay_Side_Stimulation.stimulus_x_stimulation.(fname_stat){comp,5} = P; % p-value of ttest
                    stat_struct_ttest_anova.choice_rmanova3way_SingleStimDisplay_Side_Stimulation.stimulus_x_stimulation.(fname_stat){comp,6} = stats.tstat; % t-value of ttest
                    stat_struct_ttest_anova.choice_rmanova3way_SingleStimDisplay_Side_Stimulation.stimulus_x_stimulation.(fname_stat){comp,7} = stats.df; % degrees of freedom
                end
            end
            
            %%% FDR CORRECTION %%%
            if Settings.FDR
                PValues = [stat_struct_ttest_anova.choice_rmanova3way_SingleStimDisplay_Side_Stimulation.stimulus_x_stimulation.(fname_stat){:,5}];
                
                [H_adj, ~, ~, P_adj]=fdr_bh(PValues,Settings.alpha,'pdep','no');
                H_adj(isnan(H_adj)) = 0; % if x & y are vectors of zeros or if there are not enough data H will be NaN which cannot be converted into logicals in the plot_bargraph function, so this will be set to 0, i.e. not significant; P is still NaN
                
                stat_struct_ttest_anova.choice_rmanova3way_SingleStimDisplay_Side_Stimulation.stimulus_x_stimulation.(fname_stat)(:,8) = num2cell(H_adj'); % corrected indicators of significance
                stat_struct_ttest_anova.choice_rmanova3way_SingleStimDisplay_Side_Stimulation.stimulus_x_stimulation.(fname_stat)(:,9) = num2cell(P_adj'); % corrected p-values of ttest
                
                stat_struct_ttest_anova.choice_rmanova3way_SingleStimDisplay_Side_Stimulation.stimulus_x_stimulation.(fname_stat) = ... % header if FDR
                    [{'target position','stimulation','compared to stimulation','significant','p','t','df','signif_FDR','p_FDR'};...
                    stat_struct_ttest_anova.choice_rmanova3way_SingleStimDisplay_Side_Stimulation.stimulus_x_stimulation.(fname_stat)];
                
            else
                stat_struct_ttest_anova.choice_rmanova3way_SingleStimDisplay_Side_Stimulation.stimulus_x_stimulation.(fname_stat) = ...% header if not FDR
                    [{'target position','stimulation','compared to stimulation','significant','p','t','df'};...
                    stat_struct_ttest_anova.choice_rmanova3way_SingleStimDisplay_Side_Stimulation.stimulus_x_stimulation.(fname_stat)];
            end
            
            %%%% SIDE X STIMULATION INTERACTION -> TAKE MEAN ACROSS STIMULUS %%%
            assignment_averaged = assignment(1:8,:);
            mean_idcs = {{1,9,17},{2,10,18},{3,11,19},{4,12,20},{5,13,21},{6,14,22},{7,15,23},{8,16,24}}; % which rows should be taken together and averaged across
            comparisons = {{1,2},{1,3},{1,4},{5,6},{5,7},{5,8}};
            stim_con4ttest = [stim_con;stim_con];
            side4ttest = {'contra','contra','contra','ipsi','ipsi','ipsi'};
            for ch = 1:numel(choices_stat)
                fname_stat = choices_stat{ch}; % fieldname for stat_data structure
                for i = 1:numel(mean_idcs)
                    for ses = 1:length(stat_data.(fname_stat){1,mean_idcs{i}{1}})
                        ttest_data.(fname_stat){i}(1,ses) = mean([stat_data.(fname_stat){mean_idcs{i}{1}}(ses),stat_data.(fname_stat){mean_idcs{i}{2}}(ses)]);
                        
                    end
                end
                for comp = 1:numel(comparisons) % which elements should be compared in a t test?
                    x = ttest_data.(fname_stat){comparisons{comp}{1}};
                    y = ttest_data.(fname_stat){comparisons{comp}{2}};
                    [H,P,~,stats] = ttest(x,y,'tail','both','alpha',Settings.alpha);
                    stat_struct_ttest_anova.choice_rmanova3way_SingleStimDisplay_Side_Stimulation.side_x_stimulation.(fname_stat){comp,1} = side4ttest{comp};
                    stat_struct_ttest_anova.choice_rmanova3way_SingleStimDisplay_Side_Stimulation.side_x_stimulation.(fname_stat){comp,2} = stim_con4ttest{comparisons{comp}{1}};
                    stat_struct_ttest_anova.choice_rmanova3way_SingleStimDisplay_Side_Stimulation.side_x_stimulation.(fname_stat){comp,3} = stim_con4ttest{comparisons{comp}{2}};
                    stat_struct_ttest_anova.choice_rmanova3way_SingleStimDisplay_Side_Stimulation.side_x_stimulation.(fname_stat){comp,4} = H; % is this comparison significantly different?
                    stat_struct_ttest_anova.choice_rmanova3way_SingleStimDisplay_Side_Stimulation.side_x_stimulation.(fname_stat){comp,5} = P; % p-value of ttest
                    stat_struct_ttest_anova.choice_rmanova3way_SingleStimDisplay_Side_Stimulation.side_x_stimulation.(fname_stat){comp,6} = stats.tstat; % t-value of ttest
                    stat_struct_ttest_anova.choice_rmanova3way_SingleStimDisplay_Side_Stimulation.side_x_stimulation.(fname_stat){comp,7} = stats.df; % degrees of freedom
                end
            end
            
            %%% FDR CORRECTION %%%
            if Settings.FDR
                PValues = [stat_struct_ttest_anova.choice_rmanova3way_SingleStimDisplay_Side_Stimulation.side_x_stimulation.(fname_stat){:,5}];
                
                [H_adj, ~, ~, P_adj]=fdr_bh(PValues,Settings.alpha,'pdep','no');
                H_adj(isnan(H_adj)) = 0; % if x & y are vectors of zeros or if there are not enough data H will be NaN which cannot be converted into logicals in the plot_bargraph function, so this will be set to 0, i.e. not significant; P is still NaN
                
                stat_struct_ttest_anova.choice_rmanova3way_SingleStimDisplay_Side_Stimulation.side_x_stimulation.(fname_stat)(:,8) = num2cell(H_adj'); % corrected indicators of significance
                stat_struct_ttest_anova.choice_rmanova3way_SingleStimDisplay_Side_Stimulation.side_x_stimulation.(fname_stat)(:,9) = num2cell(P_adj'); % corrected p-values of ttest
                
                stat_struct_ttest_anova.choice_rmanova3way_SingleStimDisplay_Side_Stimulation.side_x_stimulation.(fname_stat) = ... % header if FDR
                    [{'target position','stimulation','compared to stimulation','significant','p','t','df','signif_FDR','p_FDR'};...
                    stat_struct_ttest_anova.choice_rmanova3way_SingleStimDisplay_Side_Stimulation.side_x_stimulation.(fname_stat)];
                
            else
                stat_struct_ttest_anova.choice_rmanova3way_SingleStimDisplay_Side_Stimulation.side_x_stimulation.(fname_stat) = ...% header if not FDR
                    [{'target position','stimulation','compared to stimulation','significant','p','t','df'};...
                    stat_struct_ttest_anova.choice_rmanova3way_SingleStimDisplay_Side_Stimulation.side_x_stimulation.(fname_stat)];
            end
        end
        
        %%% ANOVA %%%
        if perform_ANOVA
            varnames = {'Stimulus' 'Side' 'Stimulation'};
            levnames = {{'Single Targ' 'Difficult single distr' 'Easy single distr'}, side_levels, {'None' num2str(-80) 'Go' num2str(+80)}};
            %         assignment = [1 1 1; 1 1 2; 1 1 3; 1 1 4; 1 2 1; 1 2 2; 1 2 3; 1 2 4;...
            %               2 1 1; 2 1 2; 2 1 3; 2 1 4; 2 2 1; 2 2 2; 2 2 3; 2 2 4;...
            %               3 1 1; 3 1 2; 3 1 3; 3 1 4; 3 2 1; 3 2 2; 3 2 3; 3 2 4];
            dep_var = fieldnames(stat_data);
            for dv = 1:numel(fieldnames(stat_data))
                if rep == 1 % choices
                    stat_results_anova.choice_single_stimuli.(dep_var{dv}) = threeway_rmanova(stat_data.(dep_var{dv}), assignment, varnames, levnames, Settings.alpha);
                    stat_results_anova.choice_single_stimuli.([dep_var{dv},'_stat_data']) = stat_data.(dep_var{dv});
                    stat_results_anova.choice_single_stimuli.([dep_var{dv},'_stat_data_assignment']) = assignment;
                elseif rep == 2 % error rates
                    stat_results_anova.error_rate_single_stimuli.(dep_var{dv}) = threeway_rmanova(stat_data.(dep_var{dv}), assignment, varnames, levnames, Settings.alpha);
                    stat_results_anova.error_rate_single_stimuli.([dep_var{dv},'_stat_data_assignment']) = assignment;
                end
            end
        end
    end
end
%%

%-----------------------------------------------------------------
                  %%% CHOICES (BASELINE ANOVA) %%%
%-----------------------------------------------------------------

%%% BASELINE ANOVA %%%
%%%  2-WAY ANOVA OF SINGLE AND DOUBLE STIMULUS TRIALS WITH BASELINE AND CONTROL AS LEVELS OF 2ND FACTOR %%%
%%% CHOICES %%%
if Settings.choice_anova2way_Baseline_Control
    control_data_struct_per_session = data_struct_per_session; clear('data_struct_per_session');
    load('Y:\Projects\Pulv_distractor_spatial_choice\Baseline_pre_stim\baseline_pre_stim_dataset.mat','data_struct_per_session');
    baseline_data_struct_per_session = data_struct_per_session;
    data_struct_per_session = control_data_struct_per_session;
    
    stat_data = []; assignment = []; ttest_data = [];
    perform_ANOVA = 1;  
    
    % all stimulus displays
%     tr_con2stat = {'targ_targ','distr_distr','distr_distr','targ_distr','targ_distr','targ_distr','targ_distr','single_targ','single_targ','single_distr','single_distr','single_distr','single_distr'};
%     difficulty = [1,1,2,1,1,2,2,1,1,1,1,2,2]; % first row in single_targ trials (no distr color), 1st row in single_distr (diff color), 2nd row in single_distr trials (easy color)
%     position = {'combined','combined','combined','contra','ipsi','contra','ipsi','contra','ipsi','contra','ipsi','contra','ipsi'}; % refers to tr_con2stat
    
% single stimulus displays
%     tr_con2stat = {'single_targ','single_targ','single_distr','single_distr','single_distr','single_distr'};
%     difficulty = [1,1,1,1,2,2]; % first row in single_targ trials (no distr color), 1st row in single_distr (diff color), 2nd row in single_distr trials (easy color)
%     position = {'contra','ipsi','contra','ipsi','contra','ipsi'}; % refers to tr_con2stat
    
    % paired-stimulus displays
    tr_con2stat = {'targ_targ','distr_distr','distr_distr','targ_distr','targ_distr','targ_distr','targ_distr'};
    difficulty = [1,1,2,1,1,2,2]; % first row in single_targ trials (no distr color), 1st row in single_distr (diff color), 2nd row in single_distr trials (easy color)
    position = {'combined','combined','combined','contra','ipsi','contra','ipsi'}; % refers to tr_con2stat
    
    num_tr_con2stat = numel(tr_con2stat);
    dep_var = {'fixation'}; % {'contra','fixation','ipsi'}
    group = {'baseline','control'};
    for dv = 1:numel(dep_var)
        anova_con = 0;
        for gr = 1: numel(group)
            for tr_con = 1:num_tr_con2stat
                pos = tr_con;
                col = tr_con;
                if gr == 1
                    tmp_data = [control_data_struct_per_session.(tr_con2stat{tr_con}).(position{pos}).stim_off{difficulty(col),:}]; % tmp_data with all sessions concatenated
                elseif gr == 2
                    tmp_data = [baseline_data_struct_per_session.(tr_con2stat{tr_con}).(position{pos}).stim_off{difficulty(col),:}];
                end
                for ses = 1: size(tmp_data,2)
                    anova_con = anova_con + 1; % data for different conditions are concatenated vertically
                    tmp_data2 = [tmp_data.(['true_prop_choice_per_session_',dep_var{dv}])];
                    stat_data.(dep_var{dv}).values(1,anova_con) = tmp_data2(ses); % fixation choice
                    stat_data.(dep_var{dv}).group(1,anova_con) = gr;
                    stat_data.(dep_var{dv}).tr_con(1,anova_con) = tr_con;
                end
            end
        end
        %%% POST HOC T-TEST %%%
        if perform_ANOVA && Settings.posthoc_ttest
            %%% STIMULUS X STIMULATION INTERACTION -> TAKE MEAN ACROSS SIDE %%%
            i = 0; P = []; H = [];
            for tr_con = 1:num_tr_con2stat
                for gr = 1: numel(group)
                    anova_con = 0;
                    pos = tr_con;
                    col = tr_con;
                    if gr == 1
                        tmp_data = [control_data_struct_per_session.(tr_con2stat{tr_con}).(position{pos}).stim_off{difficulty(col),:}]; % tmp_data with all sessions concatenated
                    elseif gr == 2
                        tmp_data = [baseline_data_struct_per_session.(tr_con2stat{tr_con}).(position{pos}).stim_off{difficulty(col),:}];
                    end
                    for ses = 1: size(tmp_data,2)
                        anova_con = anova_con + 1; % data for different conditions are concatenated vertically
                        tmp_data2 = [tmp_data.(['true_prop_choice_per_session_',dep_var{dv}])];
                        posthoc_data.(dep_var{dv}){gr}(anova_con) = tmp_data2(ses); % fixation choice
                    end
                    
                end
                i = i + 1;
                [H,P,~,stats] = ttest2(posthoc_data.(dep_var{dv}){1},posthoc_data.(dep_var{dv}){2});
                
                stat_struct_ttest_anova.baseline_control{i,1} = tr_con2stat{tr_con};
                stat_struct_ttest_anova.baseline_control{i,2} = difficulty(col);
                stat_struct_ttest_anova.baseline_control{i,3} = position{pos};
                stat_struct_ttest_anova.baseline_control{i,4} = dep_var{dv};
                stat_struct_ttest_anova.baseline_control{i,5} = nanmean(posthoc_data.(dep_var{dv}){1}); % group 1, i.e. baseline
                stat_struct_ttest_anova.baseline_control{i,6} = nanstd(posthoc_data.(dep_var{dv}){1}/sqrt(size(posthoc_data.(dep_var{dv}){2},1))); % group 1, i.e. baseline
                stat_struct_ttest_anova.baseline_control{i,7} = 'baseline';
                stat_struct_ttest_anova.baseline_control{i,8} = mean(posthoc_data.(dep_var{dv}){2}); % group 2, i.e. control
                stat_struct_ttest_anova.baseline_control{i,9} = nanstd(posthoc_data.(dep_var{dv}){2}/sqrt(size(posthoc_data.(dep_var{dv}){2},1))); % group 2, i.e. control
                stat_struct_ttest_anova.baseline_control{i,10} = 'control';
                stat_struct_ttest_anova.baseline_control{i,11} = H;
                stat_struct_ttest_anova.baseline_control{i,12} = P;
                stat_struct_ttest_anova.baseline_control{i,13} = stats.tstat;
                stat_struct_ttest_anova.baseline_control{i,14} = stats.df;
            end
            
            %%% FDR CORRECTION %%%
            if Settings.FDR
                PValues = [stat_struct_ttest_anova.baseline_control{:,12}];
                
                [H_adj, ~, ~, P_adj]=fdr_bh(PValues,Settings.alpha,'pdep','no');
                H_adj(isnan(H_adj)) = 0; % if x & y are vectors of zeros or if there are not enough data H will be NaN which cannot be converted into logicals in the plot_bargraph function, so this will be set to 0, i.e. not significant; P is still NaN
                
                stat_struct_ttest_anova.baseline_control(:,15) = num2cell(H_adj'); % corrected indicators of significance
                stat_struct_ttest_anova.baseline_control(:,16) = num2cell(P_adj'); % corrected p-values of ttest
                
                stat_struct_ttest_anova.baseline_control = ... % header if FDR
                    [{'tr_con','difficulty','position','choice','mean','SEM','group','mean','SEM','compared to group','significant','p','t','df','signif_FDR','p_FDR'};...
                    stat_struct_ttest_anova.baseline_control];
            else
                stat_struct_ttest_anova.baseline_control = ...% header if not FDR
                    [{'tr_con','difficulty','position','choice','mean','SEM','group','mean','SEM','compared to group','significant','p','t','df'};...
                    stat_struct_ttest_anova.baseline_control];
            end
        end
        %%% ANOVA %%%
        if perform_ANOVA
            clear y;
            % Single target design - only hit latencies %
            y = stat_data.(dep_var{dv}).values;
            F1 = stat_data.(dep_var{dv}).group; % factor 1, i.e. baseline or control
            F2 = stat_data.(dep_var{dv}).tr_con; % factor 2, i.e. stimulus displays, i.e. trial condition
            varnames = {'Group' 'Display'};
            [p,tbl,stats] = anovan(y,{F1,F2},'varnames',varnames,'model','full','display','off');
            stat_results_anova.baseline_control.(dep_var{dv}).table = tbl;
            stat_results_anova.baseline_control.(dep_var{dv}).p = p;
        end
        [results,means] = multcompare(stats,'Dimension',[1 2],'CType','bonferroni');
        set(gcf,'color','w');
    end
end


%%

%-----------------------------------------------
      %%% LATENCY & LATENCY DIFFERENCES %%%
%-----------------------------------------------

%%%  2-WAY REAPEATED MEASURES ANOVA OF SINGLE DIFFICULT DISTRACTOR TRIALS WITH SIDE AND STIMULATION CONDITIONS AS FACTORS %%%
%%%  2-WAY REAPEATED MEASURES ANOVA OF SINGLE        TARGET        TRIALS WITH SIDE AND STIMULATION CONDITIONS AS FACTORS %%%
latency_or_lat_diff = [Settings.latency_rmanova2way_SingleStimDisplay_Side_Stimulation,Settings.latency_diff_rmanova2way_SingleStimDisplay_Side_Stimulation];
for rep = 1:numel(latency_or_lat_diff)
    if latency_or_lat_diff(rep)
        stat_data.hits_latency = [];
        stat_data.error_latency = [];
        perform_ANOVA = 1; anova_con_targ = 0; anova_con_distr = 0;
        tr_con2stat = {'single_targ','single_distr'};
        num_tr_con2stat = numel(tr_con2stat);
        position = {'contra','ipsi'}; % refers to tr_con2stat
        difficulty = [1,1]; % first row in single_targ trials (no distr color), 1st row in single_distr (diff color), no easy distractor included because there are not trials in all conditions
        side_levels = [];
        if rep == 1 % latency
            i = 1;
        elseif rep == 2 % latency differences
            i = 2;
        end
        for tr_con = 1:num_tr_con2stat
            for pos = 1:numel(position); % for all positions in this trial condition
                for stim = i:num_stim_con
                    for col = difficulty(tr_con) % for the respective color (rows in data_struct) per trial condition
                        tmp_data = [data_struct_per_session.(tr_con2stat{tr_con}).(position{pos}).(stim_con{stim}){col,:}]; % tmp_data with all sessions concatenated
                        if rep == 2 % latency differences
                            tmp_data_control = [data_struct_per_session.(tr_con2stat{tr_con}).(position{pos}).stim_off{col,:}]; % tmp_data_control (stim_off) with all sessions concatenated
                        end
                        if isequal(tr_con2stat{tr_con},'single_targ') && isequal(unique([tmp_data.G_R_ratios]), 0) % if it's really single_targ trials without distractors
                            anova_con_targ = anova_con_targ + 1; % data for different conditions are concatenated vertically
                            fname = 'mean_latency_hits_per_session';
                            if rep == 1 % latency
                                stat_data.hits_latency.values(:,anova_con_targ) = [tmp_data.(fname)]; % latency of hits
                            elseif rep == 2 % latency differences
                                stat_data.hits_latency.values(:,anova_con_targ) = [tmp_data.(fname)] - [tmp_data_control.(fname)]; % latency of hits in stim condition minus latency of hits in control cond
                            end
                            stat_data.hits_latency.sessions(:,anova_con_targ) = [1:1:numel(tmp_data)]'; % subjects, i.e. sessions
                            stat_data.hits_latency.side(:,anova_con_targ) = pos*ones(numel(tmp_data),1); % factor 1, i.e. side
                            stat_data.hits_latency.stimulation(:,anova_con_targ) = stim*ones(numel(tmp_data),1); % factor 2, i.e. stimulation
                            
                        elseif isequal(tr_con2stat{tr_con},'single_distr') && isequal(unique([tmp_data.G_R_ratios]), [0.1796875;1]) % if it's really single_distr trials with the correct distractor colors
                            anova_con_distr = anova_con_distr + 1; % data for different conditions are concatenated vertically
                            fname = ['mean_latency_incorrect_',position{pos},'_choice_per_session'];
                            if rep == 1 % latency
                                stat_data.error_latency.values(:,anova_con_distr) = [tmp_data.(fname)]; % latency of errors
                            elseif rep == 2 % latency differences
                                stat_data.error_latency.values(:,anova_con_distr) = [tmp_data.(fname)] - [tmp_data_control.(fname)]; % latency of errors in stim condition minus latency of hits in control cond
                            end
                            stat_data.error_latency.sessions(:,anova_con_distr) = [1:1:numel(tmp_data)]'; % subjects, i.e. sessions
                            stat_data.error_latency.side(:,anova_con_distr) = pos*ones(numel(tmp_data),1); % factor 1, i.e. side
                            stat_data.error_latency.stimulation(:,anova_con_distr) = stim*ones(numel(tmp_data),1); % factor 2, i.e. stimulation
                            
                            % stat_data.error_latency{1,anova_con} = [tmp_data.(fname)]; % latency of incorrect choice of contra/ipsi, respectively
                        else
                            warning(sprintf('The distractor colors do not correspond to the easy and the difficult distractor. \nIn 3 way repeated measures ANOVA: single stimulus trials - side - stimulation \nANOVA is not performed!'));
                            perform_ANOVA = 0;
                        end
                    end
                end
                side_levels = [side_levels,position(pos)];
            end
        end
        if perform_ANOVA
            if rep == 1 % latency
                num_stim = 4;
            elseif rep == 2 % latency differences
                num_stim = 3;
            end
            % Single target design - only hit latencies %
            Y.hits_latency = reshape(stat_data.hits_latency.values,pos*num_stim*numel(tmp_data),1); % concatenate to one column vector
            S = reshape(stat_data.hits_latency.sessions,pos*num_stim*numel(tmp_data),1); % subjects, i.e. sessions
            F1 = reshape(stat_data.hits_latency.side,pos*num_stim*numel(tmp_data),1); % factor 1, i.e. side
            F2 = reshape(stat_data.hits_latency.stimulation,pos*num_stim*numel(tmp_data),1); % factor 2, i.e. stimulation
            varnames = {'Side' 'Stimulation'};
            if rep == 1 % latency
                stat_results_anova.lat_1targ.hits_latency = rm_anova2(Y.hits_latency,S,F1,F2,varnames);
                stat_results_anova.lat_1targ.hits_latency_stat_data = stat_data.hits_latency;
            elseif rep == 2 % latency differences
                stat_results_anova.lat_diff_1targ.hits_latency = rm_anova2(Y.hits_latency,S,F1,F2,varnames);
                stat_results_anova.lat_diff_1targ.hits_latency_stat_data = stat_data.hits_latency;
            end
            
            % Single distractor design - only error latencies %
            Y.error_latency = reshape(stat_data.error_latency.values,pos*num_stim*numel(tmp_data),1); % concatenate to one column vector
            S = reshape(stat_data.error_latency.sessions,pos*num_stim*numel(tmp_data),1); % subjects, i.e. sessions
            F1 = reshape(stat_data.error_latency.side,pos*num_stim*numel(tmp_data),1); % factor 1, i.e. side
            F2 = reshape(stat_data.error_latency.stimulation,pos*num_stim*numel(tmp_data),1); % factor 2, i.e. stimulation
            varnames = {'Side' 'Stimulation'};
            if rep == 1 % latency
                stat_results_anova.lat_1hard_distr.error_latency = rm_anova2(Y.error_latency,S,F1,F2,varnames);
                stat_results_anova.lat_1hard_distr.stat_data = stat_data.error_latency;
            elseif rep == 2 % latency differences
                stat_results_anova.lat_diff_1hard_distr.error_latency = rm_anova2(Y.error_latency,S,F1,F2,varnames);
                stat_results_anova.lat_diff_1hard_distr.stat_data = stat_data.hits_latency;
            end
        end
        %%% POST HOC T-TEST %%%
        if perform_ANOVA && Settings.posthoc_ttest
            % Compare RTs of contra to ipsi saccades; across hemifield (choice)
            fname_stat = {'absolute_lat','lat_differences'};
            if rep == 1 % latency
                comparisons = {{1,5},{2,6},{3,7},{4,8}}; %
                stim_con4ttest = stim_con;
                side4ttest = {'contra','contra','contra','contra','ipsi','ipsi','ipsi','ipsi'}; % corresponding to the values in assignment(:,2) (-> choices), but only three times, because only three comparisons
            elseif rep == 2 % latency differences
                comparisons = {{1,4},{2,5},{3,6}}; %
                stim_con4ttest = [stim_con(2:4)];
                side4ttest = {'contra','contra','contra','ipsi','ipsi','ipsi'}; % corresponding to the values in assignment(:,2) (-> choices), but only three times, because only three comparisons
            end
            for typ = 1:numel(tr_con2stat)
                mean_across_stimulation = 0; % t-test on RTs (differences) averaged across stimulation conditions to test for main effect of hemifield
                if mean_across_stimulation
                    mean_idcs = {{1,2,3},{4,5,6}};
                    comparisons = {{1,2}};
                    for i = 1:numel(mean_idcs)
                        for ses = 1:length(stat_data.hits_latency.values(1,mean_idcs{i}{1}))
                            if typ == 1
                                stat_data.hits_latency.values(i,ses) = mean([stat_data.hits_latency.values(mean_idcs{i}{1},ses),...
                                    stat_data.hits_latency.values(mean_idcs{i}{1},ses),...
                                    stat_data.hits_latency.values(mean_idcs{i}{1},ses)]);
                            elseif typ == 2
                                stat_data.error_latency.values(i,ses) = mean([stat_data.hits_latency.values(mean_idcs{i}{1},ses),...
                                    stat_data.hits_latency.values(mean_idcs{i}{1},ses),...
                                    stat_data.hits_latency.values(mean_idcs{i}{1},ses)]);
                            end
                        end
                    end
                end
                for comp = 1:numel(comparisons) % which elements should be compared in a t test?
                    if typ == 1
                        x = stat_data.hits_latency.values(:,comparisons{comp}{1});
                        y = stat_data.hits_latency.values(:,comparisons{comp}{2});
                    elseif typ == 2
                        x = stat_data.error_latency.values(:,comparisons{comp}{1});
                        y = stat_data.error_latency.values(:,comparisons{comp}{2}); 
                    end
                    [H,P,~,stats] = ttest(x,y,'tail','both','alpha',Settings.alpha);
                    stat_struct_ttest_anova.(tr_con2stat{typ}).(fname_stat{rep}){comp,1} = stim_con4ttest{comp};
                    stat_struct_ttest_anova.(tr_con2stat{typ}).(fname_stat{rep}){comp,2} = mean(x);
                    stat_struct_ttest_anova.(tr_con2stat{typ}).(fname_stat{rep}){comp,3} = nanstd(x)/sqrt(size(x,2)); % SEM
                    stat_struct_ttest_anova.(tr_con2stat{typ}).(fname_stat{rep}){comp,4} = side4ttest{comparisons{comp}{1}};
                    stat_struct_ttest_anova.(tr_con2stat{typ}).(fname_stat{rep}){comp,5} = mean(y);
                    stat_struct_ttest_anova.(tr_con2stat{typ}).(fname_stat{rep}){comp,6} = nanstd(y)/sqrt(size(y,2)); % SEM
                    stat_struct_ttest_anova.(tr_con2stat{typ}).(fname_stat{rep}){comp,7} = side4ttest{comparisons{comp}{2}};
                    stat_struct_ttest_anova.(tr_con2stat{typ}).(fname_stat{rep}){comp,8} = H; % is this comparison significantly different?
                    stat_struct_ttest_anova.(tr_con2stat{typ}).(fname_stat{rep}){comp,9} = P; % p-value of ttest
                    stat_struct_ttest_anova.(tr_con2stat{typ}).(fname_stat{rep}){comp,10} = stats.tstat; % t-value of ttest
                    stat_struct_ttest_anova.(tr_con2stat{typ}).(fname_stat{rep}){comp,11} = stats.df; % degrees of freedom
                end
                
                %%% FDR CORRECTION %%%
                if Settings.FDR
                    PValues = [stat_struct_ttest_anova.(tr_con2stat{typ}).(fname_stat{rep}){:,9}];
                    
                    [H_adj, ~, ~, P_adj]=fdr_bh(PValues,Settings.alpha,'pdep','no');
                    H_adj(isnan(H_adj)) = 0; % if x & y are vectors of zeros or if there are not enough data H will be NaN which cannot be converted into logicals in the plot_bargraph function, so this will be set to 0, i.e. not significant; P is still NaN
                    
                    stat_struct_ttest_anova.(tr_con2stat{typ}).(fname_stat{rep})(:,12) = num2cell(H_adj'); % corrected indicators of significance
                    stat_struct_ttest_anova.(tr_con2stat{typ}).(fname_stat{rep})(:,13) = num2cell(P_adj'); % corrected p-values of ttest
                    
                    stat_struct_ttest_anova.(tr_con2stat{typ}).(fname_stat{rep}) = ... % header if FDR
                    [{'target position','mean','SEM','stimulation','mean','SEM','compared to stimulation','significant','p','t','df','signif_FDR','p_FDR'};...
                    stat_struct_ttest_anova.(tr_con2stat{typ}).(fname_stat{rep})];
                    
                else
                    stat_struct_ttest_anova.(tr_con2stat{typ}).(fname_stat{rep}) = ... % header if no FDR
                    [{'target position','mean','SEM','stimulation','mean','SEM','compared to stimulation','significant','p','t','df'};...
                    stat_struct_ttest_anova.(tr_con2stat{typ}).(fname_stat{rep})];
                end
            end
        end
    end
end

%%%  3-WAY REAPEATED MEASURES ANOVA OF DOUBLE STIMULUS TRIALS WITH CHOICE AND STIMULATION CONDITIONS AS 2ND AND 3RD FACTORS %%%
latency_or_lat_diff = [Settings.latency_rmanova3way_DoubleStimDisplay_Choice_Stimulation,Settings.latency_diff_rmanova3way_DoubleStimDisplay_Choice_Stimulation];
ttest_fnames = {'latency_rmanova3way_DoubleStimDisplay_Choice_Stimulation','latency_diff_rmanova3way_DoubleStimDisplay_Choice_Stimulation'};
for rep = 1:numel(latency_or_lat_diff)
    if latency_or_lat_diff(rep)
        stat_data = []; ttest_data = [];
        perform_ANOVA = 1;
        display = 'targ_distr and targ_targ'; % 'targ_distr' % 'targ_distr and targ_targ'
        switch display
            
            case 'targ_distr and targ_targ' % both targ_distr and targ_targ trials
                tr_con2stat = {{'targ_distr','targ_targ'},{'targ_distr','targ_targ'};... % hits: difficult; easy distractor
                    {'targ_distr','distr_distr'},{'targ_distr','distr_distr'}};      % errors: diff;easy
                position = {{{'contra','ipsi'},{'combined','combined'}},{{'contra','ipsi'},{'combined','combined'}};... % refers to tr_con2stat, i.e. contra target - ipsi difficult distractor, ipsi target - contra diff distr, ...
                    {{'ipsi','contra'},{'combined','combined'}},{{'ipsi','contra'},{'combined','combined'}}}; % ... error saccades because position (of target) and choice are inverted, contra position to ipsi choice and vice versa
                difficulty = {{{1,1},{1,1}},{{2,2},{1,1}};...
                    {{1,1},{1,1}},{{2,2},{2,2}}}; % difficult: first row in targ-distr trials, easy: 2nd row in targ-distr, 1st row in targ-targ and distr-distr trials (no distr color)
                distr_G_R_ratio = [0.1796875;1];
                choices_stat = {'contra','ipsi'};
                fname_stat = {'correctSaccades_difficultDistractor','correctSaccades_easyDistractor';...
                    'errorSaccades_difficultDistractor','errorSaccades_easyDistractor'};
                
            case 'targ_distr' % targ_distr trials only
                tr_con2stat = {{'targ_distr'},{'targ_distr',}}; % hits: difficult; easy distractor
                position = {{{'contra','ipsi'}},{{'contra','ipsi'}}};% refers to tr_con2stat, i.e. contra target - ipsi difficult distractor, ipsi target - contra diff distr, ...
                difficulty = {{{1,1}},{{2,2}}};
                distr_G_R_ratio = [0.1796875;1];
                choices_stat = {'contra','ipsi'};
                fname_stat = {'correctSaccades_difficultDistractor','correctSaccades_easyDistractor'};
                
            case 'targ_targ' % targ_targ trials only
                tr_con2stat = {{'targ_targ'},{'targ_targ'}}; %;... % hits: difficult; easy distractor
                position = {{{'combined','combined'}},{{'combined','combined'}}};%... % refers to tr_con2stat, i.e. contra target - ipsi difficult distractor, ipsi target - contra diff distr, ...
                difficulty = {{{1,1}},{{1,1}}};%...
                distr_G_R_ratio = [0.1796875;1];
                choices_stat = {'contra','ipsi'};
                fname_stat = {'correctSaccades','correctSaccades'};
        end
        
        if rep == 1 % latency
            i = 1;
        elseif rep == 2 % latency differences
            i = 2;
        end
        for h_e = 1:size(tr_con2stat,1) % targ_targ or distr_distr as 2nd trial condition
            for col = 1:size(difficulty,2) % difficulty - easy or difficult distractor color
                num_tr_con2stat = numel(tr_con2stat{h_e,col});
                anova_con = 0; assignment = [];
                for tr_con = 1:num_tr_con2stat
                    for pos = 1:numel(position{h_e,col}{tr_con}); % for the respective position in this trial condition
                        for stim = i:num_stim_con
                            for ch = pos % number choice options
                                anova_con = anova_con + 1; % data for different conditions are concatenated vertically
                                tmp_data = [data_struct_per_session.(tr_con2stat{h_e,col}{tr_con}).(position{h_e,col}{tr_con}{pos}).(stim_con{stim}){difficulty{h_e,col}{tr_con}{ch},:}]; % tmp_data with all sessions concatenated
                                if rep == 2 % latency differences
                                    tmp_data_control = [data_struct_per_session.(tr_con2stat{h_e,col}{tr_con}).(position{h_e,col}{tr_con}{pos}).stim_off{difficulty{h_e,col}{tr_con}{ch},:}]; % tmp_data_control (stim_off) with all sessions concatenated
                                end
                                if isequal(unique([tmp_data.G_R_ratios]),distr_G_R_ratio) || isequal(unique([tmp_data.G_R_ratios]),0)% if distr colors are the easy and difficult distr OR if no distr is present, i.e. targ-targ
                                    fname = ['mean_latency_',choices_stat{ch},'_choice_per_session']; % fieldname for tmp_data structure
                                    if rep == 1 % latency
                                        stat_data.(fname_stat{h_e,col}){1,anova_con} = [tmp_data.(fname)]; % latency of hits
                                    elseif rep == 2 % latency differences
                                        stat_data.(fname_stat{h_e,col}){1,anova_con} = [tmp_data.(fname)] - [tmp_data_control.(fname)]; % latency of hits in stim condition minus latency of hits in control cond
                                    end
                                    if rep == 1 % latency
                                        assignment(anova_con,:) = [tr_con,ch,stim];
                                    elseif rep == 2 % latency differences
                                        assignment(anova_con,:) = [tr_con,ch,stim-1];
                                    end
                                    
                                    
                                else
                                    warning(sprintf('The distractor colors do not correspond to the easy and the difficult distractor. \nIn 3 way repeated measures ANOVA: single stimulus trials - side - stimulation \nANOVA is not performed!'));
                                    perform_ANOVA = 0;
                                end
                            end
                        end
                    end
                end
            end
        end
        %%% POST HOC T-TEST %%%
        if perform_ANOVA && Settings.posthoc_ttest 
            % Choice x Stimulation interaction -> take mean across stimulus display
            if rep == 1 % latency
                assignment_averaged = assignment(1:8,:);
                mean_idcs = {{1,9},{2,10},{3,11},{4,12},{5,13},{6,14},{7,15},{8,16}}; % which rows should be taken together and averaged across
                comparisons = {{1,2},{1,3},{1,4},{5,6},{5,7},{5,8}};
                stim_con4ttest = [stim_con;stim_con];
            elseif rep == 2 % latency differences
                assignment_averaged = assignment(1:6,:);
                mean_idcs = {{1,7},{2,8},{3,9},{4,10},{5,11},{6,12}}; % which rows should be taken together and averaged across
                comparisons = {{1,2},{1,3},{2,3},{4,5},{4,6},{5,6}}; % 
                stim_con4ttest = [stim_con(2:4);stim_con(2:4)];
            end
            % Compare stimulation conditions to each other; within hemifield (choice)
            % !! test the commented code before using it!!
            %             side4ttest = {'contra','contra','contra','ipsi','ipsi','ipsi'}; % corresponding to the values in assignment(:,2) (-> choices), but only three times, because only three comparisons
            %             fname_stat = {'correctSaccades_difficultDistractor','correctSaccades_easyDistractor','errorSaccades_difficultDistractor'}; % fieldname for stat_data structure
            %             for typ = 1:numel(fname_stat)
            %                 for i = 1:numel(mean_idcs)
            %                     for ses = 1:length(stat_data.(fname_stat{typ}){1,mean_idcs{i}{1}})
            %                         ttest_data.(fname_stat{typ}){i}(1,ses) = mean([stat_data.(fname_stat{typ}){mean_idcs{i}{1}}(ses),stat_data.(fname_stat{typ}){mean_idcs{i}{2}}(ses)]);
            %                     end
            %                 end
            %                 for comp = 1:numel(comparisons) % which elements should be compared in a t test?
            %                     x = ttest_data.(fname_stat{typ}){comparisons{comp}{1}};
            %                     y = ttest_data.(fname_stat{typ}){comparisons{comp}{2}};
            %                     [H,P,~,stats] = ttest(x,y,'tail','both','alpha',Settings.alpha);
            %                     stat_struct_ttest_anova.(ttest_fnames{rep}).choice_x_stimulation.(fname_stat{typ}){comp,1} = side4ttest{comp};
            %                     stat_struct_ttest_anova.(ttest_fnames{rep}).choice_x_stimulation.(fname_stat{typ}){comp,2} = mean(x);
            %                     stat_struct_ttest_anova.(ttest_fnames{rep}).choice_x_stimulation.(fname_stat{typ}){comp,3} = nanstd(x)/sqrt(size(x,2)); % SEM
            %                     stat_struct_ttest_anova.(ttest_fnames{rep}).choice_x_stimulation.(fname_stat{typ}){comp,4} = stim_con4ttest{comparisons{comp}{1}};
            %                     stat_struct_ttest_anova.(ttest_fnames{rep}).choice_x_stimulation.(fname_stat{typ}){comp,5} = mean(y);
            %                     stat_struct_ttest_anova.(ttest_fnames{rep}).choice_x_stimulation.(fname_stat{typ}){comp,6} = nanstd(y)/sqrt(size(y,2)); % SEM
            %                     stat_struct_ttest_anova.(ttest_fnames{rep}).choice_x_stimulation.(fname_stat{typ}){comp,7} = stim_con4ttest{comparisons{comp}{2}};
            %                     stat_struct_ttest_anova.(ttest_fnames{rep}).choice_x_stimulation.(fname_stat{typ}){comp,8} = H; % is this comparison significantly different?
            %                     stat_struct_ttest_anova.(ttest_fnames{rep}).choice_x_stimulation.(fname_stat{typ}){comp,9} = P; % p-value of ttest
            %                     stat_struct_ttest_anova.(ttest_fnames{rep}).choice_x_stimulation.(fname_stat{typ}){comp,10} = stats.tstat; % t-value of ttest
            %                     stat_struct_ttest_anova.(ttest_fnames{rep}).choice_x_stimulation.(fname_stat{typ}){comp,11} = stats.df; % degrees of freedom
            %                 end
            
            % Compare RTs of contra to ipsi saccades; across hemifield (choice)
            if rep == 1 % latency
                assignment_averaged = assignment(1:8,:);
                mean_idcs = {{1,9},{2,10},{3,11},{4,12},{5,13},{6,14},{7,15},{8,16}}; % which rows should be taken together and averaged across
                comparisons = {{1,5},{2,6},{3,7},{4,8}};
                stim_con4ttest = stim_con;
                side4ttest = {'contra','contra','contra','contra','ipsi','ipsi','ipsi','ipsi'}; % corresponding to the values in assignment(:,2) (-> choices)
            elseif rep == 2 % latency differences
                assignment_averaged = assignment(1:6,:);
                mean_idcs = {{1,7},{2,8},{3,9},{4,10},{5,11},{6,12}}; % which rows should be taken together and averaged across
                comparisons = {{1,4},{2,5},{3,6}}; %
                stim_con4ttest = [stim_con(2:4)];
                side4ttest = {'contra','contra','contra','ipsi','ipsi','ipsi'}; % corresponding to the values in assignment(:,2) (-> choices), but only three times, because only three comparisons
            end

            for typ = 1:numel(fname_stat)
                for i = 1:numel(mean_idcs)
                    if strcmp(display,'targ_targ') || strcmp(display,'targ_distr')
                        ttest_data = stat_data;
                    else
                        for ses = 1:length(stat_data.(fname_stat{typ}){1,mean_idcs{i}{1}})
                            ttest_data.(fname_stat{typ}){i}(1,ses) = mean([stat_data.(fname_stat{typ}){mean_idcs{i}{1}}(ses),stat_data.(fname_stat{typ}){mean_idcs{i}{2}}(ses)]);
                            % for separate analysis of targ-distr and targ-targ trials set mean_idcs{i}{2} to mean_idcs{i}{1}
                        end
                    end
                end
                for comp = 1:numel(comparisons) % which elements should be compared in a t test?
                    x = ttest_data.(fname_stat{typ}){comparisons{comp}{1}};
                    y = ttest_data.(fname_stat{typ}){comparisons{comp}{2}};
                    [H,P,~,stats] = ttest(x,y,'tail','both','alpha',Settings.alpha);
                    stat_struct_ttest_anova.(ttest_fnames{rep}).choice_x_stimulation.(fname_stat{typ}){comp,1} = stim_con4ttest{comp};
                    stat_struct_ttest_anova.(ttest_fnames{rep}).choice_x_stimulation.(fname_stat{typ}){comp,2} = mean(x);
                    stat_struct_ttest_anova.(ttest_fnames{rep}).choice_x_stimulation.(fname_stat{typ}){comp,3} = nanstd(x)/sqrt(size(x,2)); % SEM
                    stat_struct_ttest_anova.(ttest_fnames{rep}).choice_x_stimulation.(fname_stat{typ}){comp,4} = side4ttest{comparisons{comp}{1}};
                    stat_struct_ttest_anova.(ttest_fnames{rep}).choice_x_stimulation.(fname_stat{typ}){comp,5} = mean(y);
                    stat_struct_ttest_anova.(ttest_fnames{rep}).choice_x_stimulation.(fname_stat{typ}){comp,6} = nanstd(y)/sqrt(size(y,2)); % SEM
                    stat_struct_ttest_anova.(ttest_fnames{rep}).choice_x_stimulation.(fname_stat{typ}){comp,7} = side4ttest{comparisons{comp}{2}};
                    stat_struct_ttest_anova.(ttest_fnames{rep}).choice_x_stimulation.(fname_stat{typ}){comp,8} = H; % is this comparison significantly different?
                    stat_struct_ttest_anova.(ttest_fnames{rep}).choice_x_stimulation.(fname_stat{typ}){comp,9} = P; % p-value of ttest
                    stat_struct_ttest_anova.(ttest_fnames{rep}).choice_x_stimulation.(fname_stat{typ}){comp,10} = stats.tstat; % t-value of ttest
                    stat_struct_ttest_anova.(ttest_fnames{rep}).choice_x_stimulation.(fname_stat{typ}){comp,11} = stats.df; % degrees of freedom
                end
                
                %%% FDR CORRECTION %%%
                if Settings.FDR
                    PValues = [stat_struct_ttest_anova.(ttest_fnames{rep}).choice_x_stimulation.(fname_stat{typ}){:,9}];
                    
                    [H_adj, ~, ~, P_adj]=fdr_bh(PValues,Settings.alpha,'pdep','no');
                    H_adj(isnan(H_adj)) = 0; % if x & y are vectors of zeros or if there are not enough data H will be NaN which cannot be converted into logicals in the plot_bargraph function, so this will be set to 0, i.e. not significant; P is still NaN
                    
                    stat_struct_ttest_anova.(ttest_fnames{rep}).choice_x_stimulation.(fname_stat{typ})(:,12) = num2cell(H_adj'); % corrected indicators of significance
                    stat_struct_ttest_anova.(ttest_fnames{rep}).choice_x_stimulation.(fname_stat{typ})(:,13) = num2cell(P_adj'); % corrected p-values of ttest
                    
                    stat_struct_ttest_anova.(ttest_fnames{rep}).choice_x_stimulation.(fname_stat{typ}) = ... % header if FDR
                    [{'target position','mean','SEM','stimulation','mean','SEM','compared to stimulation','significant','p','t','df','signif_FDR','p_FDR'};...
                    stat_struct_ttest_anova.(ttest_fnames{rep}).choice_x_stimulation.(fname_stat{typ})];
                    
                else
                    stat_struct_ttest_anova.(ttest_fnames{rep}).choice_x_stimulation.(fname_stat{typ}) = ... % header if no FDR
                    [{'target position','mean','SEM','stimulation','mean','SEM','compared to stimulation','significant','p','t','df'};...
                    stat_struct_ttest_anova.(ttest_fnames{rep}).choice_x_stimulation.(fname_stat{typ})];
                end
            end
        end
        
        %%% PERFORM ANOVA %%%
        if perform_ANOVA
            varnames = {'Stimulus' 'Choice' 'Stimulation'};
            if rep == 1 % latency
                levnames = {{'Single Targ' 'Difficult single distr' 'Easy single distr'}, choices_stat, {'None' num2str(-80) 'Go' num2str(+80)}};
            elseif rep == 2 % latency differences
                levnames = {{'Single Targ' 'Difficult single distr' 'Easy single distr'}, choices_stat, {num2str(-80) 'Go' num2str(+80)}};
            end
            % assignment = [1 1 1; 1 1 2; 1 1 3; 1 1 4; 1 2 1; 1 2 2; 1 2 3; 1 2 4;...
            %       2 1 1; 2 1 2; 2 1 3; 2 1 4; 2 2 1; 2 2 2; 2 2 3; 2 2 4];
            dep_var = fieldnames(stat_data);
            for dv = 1:numel(fieldnames(stat_data))
                if rep == 1 % latency
                    stat_results_anova.latency_2stimuli.(dep_var{dv}) = threeway_rmanova(stat_data.(dep_var{dv}), assignment, varnames, levnames, Settings.alpha);
                    stat_results_anova.latency_2stimuli.([dep_var{dv},'_stat_data']) = stat_data.(dep_var{dv});
                    stat_results_anova.latency_2stimuli.([dep_var{dv},'_stat_data_assignment']) = assignment;
                elseif rep == 2 % latency differences
                    stat_results_anova.latency_diff_2stimuli.(dep_var{dv}) = threeway_rmanova(stat_data.(dep_var{dv}), assignment, varnames, levnames, Settings.alpha);
                    stat_results_anova.latency_diff_2stimuli.([dep_var{dv},'_stat_data']) = stat_data.(dep_var{dv});
                    stat_results_anova.latency_diff_2stimuli.([dep_var{dv},'_stat_data_assignment']) = assignment;
                end
            end
        end
    end
end

%%
%-----------------------------------------------
        %%% SESSION BY SESSION LATENCY %%%
%-----------------------------------------------

%%%  2-WAY REAPEATED MEASURES ANOVA OF SINGLE DIFFICULT DISTRACTOR TRIALS WITH SIDE AND STIMULATION CONDITIONS AS FACTORS %%%
%%%  2-WAY REAPEATED MEASURES ANOVA OF SINGLE        TARGET        TRIALS WITH SIDE AND STIMULATION CONDITIONS AS FACTORS %%%
if Settings.latency_anova2way_SingleStimDisplay_Side_Stimulation_ses_by_ses
    for ses = 1:12
        stat_data.hits_latency = [];
        stat_data.error_latency = [];
        perform_ANOVA = 1; anova_con_targ = 0; anova_con_distr = 0;
        tr_con2stat = {'single_targ','single_distr'};
        num_tr_con2stat = numel(tr_con2stat);
        position = {'contra','ipsi'}; % refers to tr_con2stat
        difficulty = [1,1]; % first row in single_targ trials (no distr color), 1st row in single_distr (diff color), no easy distractor included because there are not enough trials in all conditions
        i = 2;
        for tr_con = 1:num_tr_con2stat
            for pos = 1:numel(position); % for all positions in this trial condition
                for stim = i:num_stim_con
                    for col = difficulty(tr_con) % for the respective color (rows in data_struct) per trial condition
                        tmp_data = [data_struct_per_session.(tr_con2stat{tr_con}).(position{pos}).(stim_con{stim}){col,:}]; % tmp_data with all sessions concatenated
                        tmp_data_control = [data_struct_per_session.(tr_con2stat{tr_con}).(position{pos}).stim_off{col,:}]; % tmp_data with all sessions concatenated
                        if isequal(tr_con2stat{tr_con},'single_targ') && isequal(unique([tmp_data.G_R_ratios]), 0) % if it's really single_targ trials without distractors
                            fname = 'latencies_hits_per_session';
                            fname_control = 'mean_latency_hits_per_session';
                            num_saccades = numel([tmp_data(ses).(fname)]);
                            for lat = 1: num_saccades
                                anova_con_targ = anova_con_targ + 1; % data for different conditions are concatenated vertically
                                
                                % calculate the difference of any single saccade latency in this (stimulation) ondition to the mean of control saccade latencies (stim_off) in this condition
                                
                                stat_data.hits_latency.session_by_session{1,ses}.values(1,anova_con_targ) = tmp_data(ses).(fname)(lat);% - tmp_data_control(ses).(fname_control); % latency of hits in stim condition minus latency of hits in control cond
                                
                                % stat_data_tmp.hits_latency{1,ses}.sessions(:,anova_con_targ) = [1:1:numel(tmp_data)]'; % subjects, i.e. sessions
                                stat_data.hits_latency.session_by_session{1,ses}.side(1,anova_con_targ) = pos*ones(numel(num_saccades),1); % factor 1, i.e. side
                                stat_data.hits_latency.session_by_session{1,ses}.stimulation(1,anova_con_targ) = stim*ones(numel(num_saccades),1); % factor 2, i.e. stimulation
                            end
                        elseif isequal(tr_con2stat{tr_con},'single_distr') && isequal(unique([tmp_data.G_R_ratios]), [0.1796875;1]) % if it's really single_distr trials with the correct distractor colors
                            fname = ['latencies_incorrect_',position{pos},'_choice_per_session'];
                            fname_control = ['mean_latency_incorrect_',position{pos},'_choice_per_session'];
                            num_saccades = numel([tmp_data(ses).(fname)]);
                            for lat = 1: num_saccades
                                anova_con_distr = anova_con_distr + 1; % data for different conditions are concatenated vertically
                                stat_data.error_latency.session_by_session{1,ses}.values(1,anova_con_distr) = tmp_data(ses).(fname)(lat);% - tmp_data_control(ses).(fname_control); % latency of errors in stim condition minus latency of errors in control cond
                                
                                % stat_data_tmp.error_latency.sessions(:,anova_con_distr) = [1:1:numel(tmp_data)]'; % subjects, i.e. sessions
                                stat_data.error_latency.session_by_session{1,ses}.side(1,anova_con_distr) = pos*ones(numel(num_saccades),1); % factor 1, i.e. side
                                stat_data.error_latency.session_by_session{1,ses}.stimulation(1,anova_con_distr) = stim*ones(numel(num_saccades),1); % factor 2, i.e. stimulation
                                
                            end
                        else
                            warning(sprintf('The distractor colors do not correspond to the easy and the difficult distractor. \nIn 3 way repeated measures ANOVA: single stimulus trials - side - stimulation \nANOVA is not performed!'));
                            perform_ANOVA = 0;
                        end
                    end
                end
                if Settings.latency_kruskal_wallis_Session_by_session
                    if isequal(tr_con2stat{tr_con},'single_targ') && isequal(unique([tmp_data.G_R_ratios]), 0) % if it's really single_targ trials without distractors
                        Y = stat_data.hits_latency.session_by_session{1,ses}.values;
                        group = stat_data.hits_latency.session_by_session{1,ses}.stimulation; % factor 1, i.e. trial condition
                        [p,tbl] = kruskalwallis(Y,group,'off');
                        stat_results_anova.kruskalwallis.(tr_con2stat{tr_con}).(position{pos}).tables{1,ses} = tbl;
                        stat_results_anova.kruskalwallis.(tr_con2stat{tr_con}).(position{pos}).p_values.(['Session_',num2str(ses)]) = p;
                    elseif isequal(tr_con2stat{tr_con},'single_distr') && isequal(unique([tmp_data.G_R_ratios]), [0.1796875;1]) % if it's really single_distr trials with the correct distractor colors
                        Y = stat_data.error_latency.session_by_session{1,ses}.values;
                        group = stat_data.error_latency.session_by_session{1,ses}.stimulation; % factor 1, i.e. trial condition
                        [p,tbl] = kruskalwallis(Y,group,'off');
                        stat_results_anova.kruskalwallis.(tr_con2stat{tr_con}).(position{pos}).tables{1,ses} = tbl;
                        stat_results_anova.kruskalwallis.(tr_con2stat{tr_con}).(position{pos}).p_values.(['Session_',num2str(ses)]) = p;
                    end
                   
                end
            end
        end
        if perform_ANOVA
            num_stim = 3; clear y;
            % Single target design - only hit latencies %
            y.hits_latency = stat_data.hits_latency.session_by_session{1,ses}.values;
            F1 = stat_data.hits_latency.session_by_session{1,ses}.side; % factor 1, i.e. side
            F2 = stat_data.hits_latency.session_by_session{1,ses}.stimulation; % factor 2, i.e. stimulation
            varnames = {'Side' 'Stimulation'};
            [p,tbl] = anovan(y.hits_latency,{F1,F2},'varnames',varnames,'model','full','display','off');
            stat_results_anova.lat_1targ_ses_by_ses.hits_latency_tables{1,ses} = tbl;
            stat_results_anova.lat_1targ_ses_by_ses.hits_latency_summary.(['Session_',num2str(ses)]) = p;
            
            % Single distractor design - only error latencies %
            y.error_latency = stat_data.error_latency.session_by_session{1,ses}.values;
            F1 = stat_data.error_latency.session_by_session{1,ses}.side; % factor 1, i.e. side
            F2 = stat_data.error_latency.session_by_session{1,ses}.stimulation; % factor 2, i.e. stimulation
            varnames = {'Side' 'Stimulation'};
            [p,tbl] = anovan(y.error_latency,{F1,F2},'varnames',varnames,'model','full','display','off');
            stat_results_anova.lat_1hard_distr_ses_by_ses.error_latency_tables{1,ses} = tbl;
            stat_results_anova.lat_1hard_distr_ses_by_ses.error_latency_summary.(['Session_',num2str(ses)]) = p;
        end
    end
end

%%%  3-WAY REAPEATED MEASURES ANOVA OF DOUBLE STIMULUS TRIALS WITH CHOICE AND STIMULATION CONDITIONS AS 2ND AND 3RD FACTORS %%%
if Settings.latency_anova3way_DoubleStimDisplay_Choice_Stimulation_ses_by_ses

    for ses = 1:12
        stat_data = [];
        perform_ANOVA = 1;
        tr_con2stat = {{'targ_distr','targ_targ'},{'targ_distr','targ_targ'};... % hits: difficult; easy distractor
            {'targ_distr','distr_distr'},{'targ_distr','distr_distr'}};      % errors: diff;easy
        position = {{{'contra','ipsi'},{'combined','combined'}},{{'contra','ipsi'},{'combined','combined'}};... % refers to tr_con2stat, i.e. contra target - ipsi difficult distractor, ipsi target - contra diff distr, ...
            {{'ipsi','contra'},{'combined','combined'}},{{'ipsi','contra'},{'combined','combined'}}}; % ... error saccades because position (of target) and choice are inverted, contra position to ipsi choice and vice versa
        difficulty = {{{1,1},{1,1}},{{2,2},{1,1}};...
            {{1,1},{1,1}},{{2,2},{2,2}}}; % difficult: first row in targ-distr trials, easy: 2nd row in targ-distr, 1st row in targ-targ and distr-distr trials (no distr color)
        distr_G_R_ratio = [0.1796875;1];
        choices_stat = {'contra','ipsi'};
        fname_stat = {'correctSaccades_difficultDistractor_session_by_session','correctSaccades_easyDistractor_session_by_session';...
            'errorSaccades_difficultDistractor_session_by_session','errorSaccades_easyDistractor_session_by_session'};
        
        i = 1; % start with control stimulation condition
        for h_e = 1:size(tr_con2stat,1) % targ_targ or distr_distr as 2nd trial condition
            for col = 1:size(difficulty,2) % difficulty - easy or difficult distractor color
                num_tr_con2stat = numel(tr_con2stat{h_e,col});
                anova_con = 0; assignment = [];
                for tr_con = 1:num_tr_con2stat
                    for pos = 1:numel(position{h_e,col}{tr_con}); % for the respective position in this trial condition
                        for stim = i:num_stim_con
                            ch = pos; % number choice options
                            
                            tmp_data = [data_struct_per_session.(tr_con2stat{h_e,col}{tr_con}).(position{h_e,col}{tr_con}{pos}).(stim_con{stim}){difficulty{h_e,col}{tr_con}{ch},:}]; % tmp_data with all sessions concatenated
                            tmp_data_control = [data_struct_per_session.(tr_con2stat{h_e,col}{tr_con}).(position{h_e,col}{tr_con}{pos}).stim_off{difficulty{h_e,col}{tr_con}{ch},:}]; % tmp_data with all sessions concatenated
                            if isequal(unique([tmp_data.G_R_ratios]),distr_G_R_ratio) || isequal(unique([tmp_data.G_R_ratios]),0)% if distr colors are the easy and difficult distr OR if no distr is present, i.e. targ-targ
                                fname = ['latencies_',choices_stat{ch},'_choice_per_session']; % fieldname for tmp_data structure
                                fname_control = ['mean_latency_',choices_stat{ch},'_choice_per_session'];
                                num_saccades = numel([tmp_data(ses).(fname)]);
                                if num_saccades == 0;
                                    stat_data.(fname_stat{h_e,col}){1,ses}.values = NaN; % latency in stim condition minus latency in control cond
                                    stat_data.(fname_stat{h_e,col}){1,ses}.tr_con = NaN; % factor 1, i.e. side
                                    stat_data.(fname_stat{h_e,col}){1,ses}.choice = NaN; % factor 1, i.e. side
                                    stat_data.(fname_stat{h_e,col}){1,ses}.stimulation = NaN; % factor 2, i.e. stimulation
                                else
                                    for lat = 1: num_saccades
                                        anova_con = anova_con + 1; % data for different conditions are concatenated vertically
                                        
                                        stat_data.(fname_stat{h_e,col}){1,ses}.values(1,anova_con) = tmp_data(ses).(fname)(lat) - tmp_data_control(ses).(fname_control); % latency in stim condition minus latency in control cond
                                        stat_data.(fname_stat{h_e,col}){1,ses}.tr_con(1,anova_con) = tr_con*ones(numel(num_saccades),1); % factor 1, i.e. side
                                        stat_data.(fname_stat{h_e,col}){1,ses}.choice(1,anova_con) = ch*ones(numel(num_saccades),1); % factor 1, i.e. side
                                        stat_data.(fname_stat{h_e,col}){1,ses}.stimulation(1,anova_con) = stim*ones(numel(num_saccades),1); % factor 2, i.e. stimulation
                                    end
                                end
                                
                            else
                                warning(sprintf('The distractor colors do not correspond to the easy and the difficult distractor. \nIn 3 way repeated measures ANOVA: single stimulus trials - side - stimulation \nANOVA is not performed!'));
                                perform_ANOVA = 0;
                            end
                        end
                    end
                end
            end
        end
        
        %%% POST HOC MANN-WHITNEY-U TEST %%%
        if perform_ANOVA && Settings.posthoc_ttest % NOT WORKING
            % Choice x Stimulation interaction -> take mean across stimulus               
            assignment_averaged = assignment(1:8,:);
            mean_idcs = {{1,9},{2,10},{3,11},{4,12},{5,13},{6,14},{7,15},{8,16}}; % which rows should be taken together and averaged across
            comparisons = {{1,2},{1,3},{1,4},{5,6},{5,7},{5,8}};
            stim_con4ttest = [stim_con;stim_con];
            side4ttest = {'contra','contra','contra','ipsi','ipsi','ipsi'};
            fname_stat = {'correctSaccades_difficultDistractor','correctSaccades_easyDistractor'}; % fieldname for stat_data structure
            for col = 1:numel(fname_stat)
                for i = 1:numel(mean_idcs)
                    for ses = 1:length(stat_data.(fname_stat{col}){1,mean_idcs{i}{1}})
                        ttest_data.(fname_stat{col}){i}(1,ses) = mean([stat_data.(fname_stat{col}){mean_idcs{i}{1}}(ses),stat_data.(fname_stat{col}){mean_idcs{i}{2}}(ses)]);
                        
                    end
                end
                for comp = 1:numel(comparisons) % which elements should be compared in a t test?
                    x = ttest_data.(fname_stat{col}){comparisons{comp}{1}};
                    y = ttest_data.(fname_stat{col}){comparisons{comp}{2}};
                    [H,P,~,stats] = ttest(x,y,'tail','both','alpha',Settings.alpha);
                    stat_struct_ttest_anova.latency_rmanova3way_DoubleStimDisplay_Choice_Stimulation.choice_x_stimulation.(fname_stat{col}){comp,1} = side4ttest{comp};
                    stat_struct_ttest_anova.latency_rmanova3way_DoubleStimDisplay_Choice_Stimulation.choice_x_stimulation.(fname_stat{col}){comp,2} = stim_con4ttest{comparisons{comp}{1}};
                    stat_struct_ttest_anova.latency_rmanova3way_DoubleStimDisplay_Choice_Stimulation.choice_x_stimulation.(fname_stat{col}){comp,3} = stim_con4ttest{comparisons{comp}{2}};
                    stat_struct_ttest_anova.latency_rmanova3way_DoubleStimDisplay_Choice_Stimulation.choice_x_stimulation.(fname_stat{col}){comp,4} = H; % is this comparison significantly different?
                    stat_struct_ttest_anova.latency_rmanova3way_DoubleStimDisplay_Choice_Stimulation.choice_x_stimulation.(fname_stat{col}){comp,5} = P; % p-value of ttest
                    stat_struct_ttest_anova.latency_rmanova3way_DoubleStimDisplay_Choice_Stimulation.choice_x_stimulation.(fname_stat{col}){comp,6} = stats.tstat; % t-value of ttest
                    stat_struct_ttest_anova.latency_rmanova3way_DoubleStimDisplay_Choice_Stimulation.choice_x_stimulation.(fname_stat{col}){comp,7} = stats.df; % degrees of freedom
                end
                stat_struct_ttest_anova.latency_rmanova3way_DoubleStimDisplay_Choice_Stimulation.choice_x_stimulation.(fname_stat{col}) = ...
                    [{'target position','stimulation','compared to stimulation','significant','p','t','df'};...
                    stat_struct_ttest_anova.latency_rmanova3way_DoubleStimDisplay_Choice_Stimulation.choice_x_stimulation.(fname_stat{col})];
                
            end
        end
        
        %%% PERFORM ANOVA %%%
        if perform_ANOVA
            dep_var = fieldnames(stat_data);
            for dv = 1:numel(fieldnames(stat_data))
                y = stat_data.(dep_var{dv}){1,ses}.values;
                F1 = stat_data.(dep_var{dv}){1,ses}.tr_con; % factor 1, i.e. trial condition
                F2 = stat_data.(dep_var{dv}){1,ses}.choice; % factor 2, i.e. choice
                F3 = stat_data.(dep_var{dv}){1,ses}.stimulation; % factor 3, i.e. stimulation
                varnames = {'Stimulus' 'Choice' 'Stimulation'};
                if ~isnan(y)
                    [p,tbl] = anovan(y,{F1,F2,F3},'varnames',varnames,'model','full','display','off');
                    stat_results_anova.lat_2stimuli_ses_by_ses.([dep_var{dv},'_tables']){1,ses} = tbl;
                    stat_results_anova.lat_2stimuli_ses_by_ses.([dep_var{dv},'_p_values']).(['Session_',num2str(ses)]) = p;
                end
                
                effects = {'Stimulus';'Choice';'Stimulation';'Stimulus*Choice';'Stimulus*Stimulation';'Choice*Stimulation';'Stimulus*Choice*Stimulation'};
            end
        end
    end
end

%%

% ------------------
%%% Kruskal-Wallis tests on absolute latencies with Posthoc Mann-Whitney-U tests %%%
% ------------------

%%%  3-WAY REAPEATED MEASURES ANOVA OF DOUBLE STIMULUS TRIALS WITH CHOICE AND STIMULATION CONDITIONS AS 2ND AND 3RD FACTORS %%%
if Settings.latency_kruskal_wallis_Session_by_session
    
    lat_ses_by_ses2plot= struct('minus_significant',struct('minus_80',0,'Go',0,'plus_80',0),...
        'plus_significant',struct('minus_80',0,'Go',0,'plus_80',0),...
        'minus',struct('minus_80',0,'Go',0,'plus_80',0),...
        'plus',struct('minus_80',0,'Go',0,'plus_80',0));
    min_num_trials_for_kruskal_wallis = 6;
    for ses = 1:12
        stat_data = [];
        perform_ANOVA = 1;
        tr_con2stat = {{'targ_distr','targ_targ'},{'targ_distr','targ_targ'};... % hits: difficult; easy distractor
            {'targ_distr','distr_distr'},{'targ_distr','distr_distr'}};      % errors: diff;easy
        position = {{{'contra','ipsi'},{'combined','combined'}},{{'contra','ipsi'},{'combined','combined'}};... % refers to tr_con2stat, i.e. contra target - ipsi difficult distractor, ipsi target - contra diff distr, ...
            {{'ipsi','contra'},{'combined','combined'}},{{'ipsi','contra'},{'combined','combined'}}}; % ... error saccades because position (of target) and choice are inverted, contra position to ipsi choice and vice versa
        difficulty = {{{1,1},{1,1}},{{2,2},{1,1}};...
            {{1,1},{1,1}},{{2,2},{2,2}}}; % difficult: first row in targ-distr trials, easy: 2nd row in targ-distr, 1st row in targ-targ and distr-distr trials (no distr color)
        distr_G_R_ratio = [0.1796875;1];
        choices_stat = {'contra','ipsi'};
        fname_stat = {'correctSaccades_difficultDistractor_session_by_session','correctSaccades_easyDistractor_session_by_session';...
            'errorSaccades_difficultDistractor_session_by_session','errorSaccades_easyDistractor_session_by_session'};
        
        i = 1; % start with control stimulation condition
        for h_e = 1:size(tr_con2stat,1) % targ_targ or distr_distr as 2nd trial condition
            for col = 1:size(difficulty,2) % difficulty - easy or difficult distractor color
                num_tr_con2stat = numel(tr_con2stat{h_e,col});
                anova_con = 0; assignment = [];
                for tr_con = 1:num_tr_con2stat
                    for pos = 1:numel(position{h_e,col}{tr_con}); % for the respective position in this trial condition
                        kruwal_con = 0;
                        lat_ses_by_ses2plot.(fname_stat{h_e,col}).(tr_con2stat{h_e,col}{tr_con}).(position{h_e,col}{tr_con}{pos}) = ...
                            struct('minus_significant',struct('minus_80',0,'Go',0,'plus_80',0),...
                            'plus_significant',struct('minus_80',0,'Go',0,'plus_80',0),...
                            'minus',struct('minus_80',0,'Go',0,'plus_80',0),...
                            'plus',struct('minus_80',0,'Go',0,'plus_80',0));
                        for stim = i:num_stim_con
                            ch = pos; % number choice options
                            
                            tmp_data = [data_struct_per_session.(tr_con2stat{h_e,col}{tr_con}).(position{h_e,col}{tr_con}{pos}).(stim_con{stim}){difficulty{h_e,col}{tr_con}{ch},:}]; % tmp_data with all sessions concatenated
                            tmp_data_control = [data_struct_per_session.(tr_con2stat{h_e,col}{tr_con}).(position{h_e,col}{tr_con}{pos}).stim_off{difficulty{h_e,col}{tr_con}{ch},:}]; % tmp_data with all sessions concatenated
                            if isequal(unique([tmp_data.G_R_ratios]),distr_G_R_ratio) || isequal(unique([tmp_data.G_R_ratios]),0)% if distr colors are the easy and difficult distr OR if no distr is present, i.e. targ-targ
                                fname = ['latencies_',choices_stat{ch},'_choice_per_session']; % fieldname for tmp_data structure
                                fname_control = ['mean_latency_',choices_stat{ch},'_choice_per_session'];
                                num_saccades = numel([tmp_data(ses).(fname)]);
                                if num_saccades == 0;
                                    stat_data.(fname_stat{h_e,col}){1,ses}.values = NaN; % latency in stim condition minus latency in control cond
                                    stat_data.(fname_stat{h_e,col}){1,ses}.tr_con = NaN; % factor 1, i.e. side
                                    stat_data.(fname_stat{h_e,col}){1,ses}.choice = NaN; % factor 1, i.e. side
                                    stat_data.(fname_stat{h_e,col}){1,ses}.stimulation = NaN; % factor 2, i.e. stimulation
                                    
                                    kruwal_con = kruwal_con + 1; % concatentation of data for kruskal-wallis test
                                    kruwal_data.(fname_stat{h_e,col}){1,ses}.values(1,kruwal_con) = NaN;% - tmp_data_control(ses).(fname_control); % latency in stim condition minus latency in control cond
                                    kruwal_data.(fname_stat{h_e,col}){1,ses}.tr_con(1,kruwal_con) = tr_con*ones(numel(num_saccades),1); % factor 1, i.e. side
                                    kruwal_data.(fname_stat{h_e,col}){1,ses}.choice(1,kruwal_con) = ch*ones(numel(num_saccades),1); % factor 1, i.e. side
                                    kruwal_data.(fname_stat{h_e,col}){1,ses}.stimulation(1,kruwal_con) = stim*ones(numel(num_saccades),1); % factor 2, i.e. stimulation
                                else
                                    for lat = 1: num_saccades
                                        kruwal_con = kruwal_con + 1; % concatentation of data for kruskal-wallis test
                                        anova_con = anova_con + 1; % data for different conditions are concatenated vertically
                                        
                                        stat_data.(fname_stat{h_e,col}){1,ses}.values(1,anova_con) = tmp_data(ses).(fname)(lat);% - tmp_data_control(ses).(fname_control); % latency in stim condition minus latency in control cond
                                        stat_data.(fname_stat{h_e,col}){1,ses}.tr_con(1,anova_con) = tr_con*ones(numel(num_saccades),1); % factor 1, i.e. side
                                        stat_data.(fname_stat{h_e,col}){1,ses}.choice(1,anova_con) = ch*ones(numel(num_saccades),1); % factor 1, i.e. side
                                        stat_data.(fname_stat{h_e,col}){1,ses}.stimulation(1,anova_con) = stim*ones(numel(num_saccades),1); % factor 2, i.e. stimulation
                                        % assignment(anova_con,:) = [tr_con,ch,stim];
                                        kruwal_data.(fname_stat{h_e,col}){1,ses}.values(1,kruwal_con) = tmp_data(ses).(fname)(lat);% - tmp_data_control(ses).(fname_control); % latency in stim condition minus latency in control cond
                                        kruwal_data.(fname_stat{h_e,col}){1,ses}.tr_con(1,kruwal_con) = tr_con*ones(numel(num_saccades),1); % factor 1, i.e. side
                                        kruwal_data.(fname_stat{h_e,col}){1,ses}.choice(1,kruwal_con) = ch*ones(numel(num_saccades),1); % factor 1, i.e. side
                                        kruwal_data.(fname_stat{h_e,col}){1,ses}.stimulation(1,kruwal_con) = stim*ones(numel(num_saccades),1); % factor 2, i.e. stimulation
                                    end
                                end
                                
                            else
                                warning(sprintf('The distractor colors do not correspond to the easy and the difficult distractor. \nIn 3 way repeated measures ANOVA: single stimulus trials - side - stimulation \nANOVA is not performed!'));
                                perform_ANOVA = 0;
                            end
                        end
                        if Settings.latency_kruskalwallis
                            Y = kruwal_data.(fname_stat{h_e,col}){1,ses}.values;
                            group = kruwal_data.(fname_stat{h_e,col}){1,ses}.stimulation; % factor 1, i.e. trial condition
                            [p,tbl] = kruskalwallis(Y,group,'off');
                            stat_results_anova.kruskalwallis.(tr_con2stat{h_e,col}{tr_con}).(position{h_e,col}{tr_con}{pos})...
                                .([choices_stat{ch},'_choice']).tables{difficulty{h_e,col}{tr_con}{ch},ses} = tbl;
                            stat_results_anova.kruskalwallis.(tr_con2stat{h_e,col}{tr_con}).(position{h_e,col}{tr_con}{pos})...
                                .([choices_stat{ch},'_choice']).p_values{difficulty{h_e,col}{tr_con}{ch},1}(ses,1) = ses; % .(['Session_',num2str(ses)])
                            stat_results_anova.kruskalwallis.(tr_con2stat{h_e,col}{tr_con}).(position{h_e,col}{tr_con}{pos})...
                                .([choices_stat{ch},'_choice']).p_values{difficulty{h_e,col}{tr_con}{ch},1}(ses,2) = p;
                            sessions = stat_results_anova.kruskalwallis.(tr_con2stat{h_e,col}{tr_con}).(position{h_e,col}{tr_con}{pos})...
                                .([choices_stat{ch},'_choice']).p_values{difficulty{h_e,col}{tr_con}{ch},1}(:,1);
                            p_values = stat_results_anova.kruskalwallis.(tr_con2stat{h_e,col}{tr_con}).(position{h_e,col}{tr_con}{pos})...
                                .([choices_stat{ch},'_choice']).p_values{difficulty{h_e,col}{tr_con}{ch},1}(:,2);
                            stat_results_anova.kruskalwallis.(tr_con2stat{h_e,col}{tr_con}).(position{h_e,col}{tr_con}{pos})...
                                .([choices_stat{ch},'_choice']).N_sign_sessions{difficulty{h_e,col}{tr_con}{ch},1} = sum(p_values < Settings.alpha);
                            stat_results_anova.kruskalwallis.(tr_con2stat{h_e,col}{tr_con}).(position{h_e,col}{tr_con}{pos})...
                                .([choices_stat{ch},'_choice']).Nonsign_sessions{difficulty{h_e,col}{tr_con}{ch},1} = sessions(p_values > Settings.alpha);
                            
                            %%% POST HOC MANN-WHITNEY-U TEST %%%
                            if ~isnan(p) && Settings.latency_posthoc_Mann_Whitney_U % p is NaN if kruskal wallis could not be calculated, e.g. if not all stim condition are present in this session
                                clear x; clear y; clear p_kw; perform_kw = 0;
                                stim_idx = unique(group);
                                if any(stim_idx == 1) % if there is a control condition
                                    if ~stim_idx(1) == 1 % control should always be the first condition
                                        error('Control condition is not in the first position');
                                    elseif sum(group==stim_idx(1)) >= min_num_trials_for_kruskal_wallis % both control conditions...
                                        for stim = 2:numel(stim_idx)
                                            if sum(group==stim_idx(stim)) >= min_num_trials_for_kruskal_wallis % ... and stimulation condition must have at least 6 trials per condition to be tested
                                                
                                                control = 1;
                                                x{stim-1} = Y(group==control); % control stimulation condition = stim_off
                                                y{stim-1} = Y(group==stim_idx(stim));
                                                if ~ isempty(x) || ~isempty(y)
                                                    [p_kw(stim-1),h(stim-1)] = ranksum(x{stim-1},y{stim-1}); % perform the MANN-WHITNEY-U Test
                                                end
                                                perform_kw = 1;
                                            end
                                        end
                                        %%% FDR CORRECTION %%%
                                        if perform_kw
                                            condition_subfields = {'stim_off','minus_80','Go','plus_80'};
                                            [H_adj, crit_p, adj_ci_cvrg, P_adj]=fdr_bh(p_kw,Settings.alpha,'pdep','no');
                                            H_adj(isnan(H_adj)) = 0; % if x & y are vectors of zeros or if there are not enough data H will be NaN which cannot be converted into logicals in the plot_bargraph function, so this will be set to 0, i.e. not significant; P is still NaN
                                            for stim = 1:numel(P_adj)
                                                if mean(x{stim}) > mean(y{stim})
                                                    direction = 'minus';
                                                    if H_adj(stim) % increase count of significant sessions
                                                        lat_ses_by_ses2plot.(tr_con2stat{h_e,col}{tr_con}).([direction,'_significant']).(condition_subfields{stim_idx(stim)})(stim) = ...
                                                            lat_ses_by_ses2plot.(tr_con2stat{h_e,col}{tr_con}).(direction).(condition_subfields{stim_idx(stim)})(stim) + 1; % count number of sessions
                                                    end
                                                elseif mean(x{stim}) < mean(y{stim})
                                                    direction = 'plus';
                                                    if H_adj(stim) % increase count of significant sessions
                                                        lat_ses_by_ses2plot.(tr_con2stat{h_e,col}{tr_con}).([direction,'_significant']).(condition_subfields{stim_idx(stim)})(stim) = ...
                                                            lat_ses_by_ses2plot.(tr_con2stat{h_e,col}{tr_con}).(direction).(condition_subfields{stim_idx(stim)})(stim) + 1; % count number of sessions
                                                    end
                                                end
                                                % increase count of total sessions
                                                lat_ses_by_ses2plot.(tr_con2stat{h_e,col}{tr_con}).(direction).(condition_subfields{stim_idx(stim)})(stim) = ...
                                                    lat_ses_by_ses2plot.(tr_con2stat{h_e,col}{tr_con}).(direction).(condition_subfields{stim_idx(stim)})(stim) + 1; % count number of sessions
                                                
                                            end
                                        end
                                    end
                                end
                                
                            end
                        end
                    end
                end
            end
        end
        
        %%% POST HOC MANN-WHITNEY-U TEST %%%
        if perform_ANOVA && Settings.posthoc_ttest % NOT WORKING
            % Choice x Stimulation interaction -> take mean across stimulus
            %[H_adj, crit_p, adj_ci_cvrg, P_adj]=fdr_bh(PValues,Settings.alpha,'pdep','no');
            %H_adj(isnan(H_adj)) = 0; % if x & y are vectors of zeros or if there are not enough data H will be NaN which cannot be converted into logicals in the plot_bargraph function, so this will be set to 0, i.e. not significant; P is still NaN
                            
            assignment_averaged = assignment(1:8,:);
            mean_idcs = {{1,9},{2,10},{3,11},{4,12},{5,13},{6,14},{7,15},{8,16}}; % which rows should be taken together and averaged across
            comparisons = {{1,2},{1,3},{1,4},{5,6},{5,7},{5,8}};
            stim_con4ttest = [stim_con;stim_con];
            side4ttest = {'contra','contra','contra','ipsi','ipsi','ipsi'};
            fname_stat = {'correctSaccades_difficultDistractor','correctSaccades_easyDistractor'}; % fieldname for stat_data structure
            for col = 1:numel(fname_stat)
                for i = 1:numel(mean_idcs)
                    for ses = 1:length(stat_data.(fname_stat{col}){1,mean_idcs{i}{1}})
                        ttest_data.(fname_stat{col}){i}(1,ses) = mean([stat_data.(fname_stat{col}){mean_idcs{i}{1}}(ses),stat_data.(fname_stat{col}){mean_idcs{i}{2}}(ses)]);
                        
                    end
                end
                for comp = 1:numel(comparisons) % which elements should be compared in a t test?
                    x = ttest_data.(fname_stat{col}){comparisons{comp}{1}};
                    y = ttest_data.(fname_stat{col}){comparisons{comp}{2}};
                    [H,P,~,stats] = ttest(x,y,'tail','both','alpha',Settings.alpha);
                    stat_struct_ttest_anova.latency_rmanova3way_DoubleStimDisplay_Choice_Stimulation.choice_x_stimulation.(fname_stat{col}){comp,1} = side4ttest{comp};
                    stat_struct_ttest_anova.latency_rmanova3way_DoubleStimDisplay_Choice_Stimulation.choice_x_stimulation.(fname_stat{col}){comp,2} = stim_con4ttest{comparisons{comp}{1}};
                    stat_struct_ttest_anova.latency_rmanova3way_DoubleStimDisplay_Choice_Stimulation.choice_x_stimulation.(fname_stat{col}){comp,3} = stim_con4ttest{comparisons{comp}{2}};
                    stat_struct_ttest_anova.latency_rmanova3way_DoubleStimDisplay_Choice_Stimulation.choice_x_stimulation.(fname_stat{col}){comp,4} = H; % is this comparison significantly different?
                    stat_struct_ttest_anova.latency_rmanova3way_DoubleStimDisplay_Choice_Stimulation.choice_x_stimulation.(fname_stat{col}){comp,5} = P; % p-value of ttest
                    stat_struct_ttest_anova.latency_rmanova3way_DoubleStimDisplay_Choice_Stimulation.choice_x_stimulation.(fname_stat{col}){comp,6} = stats.tstat; % t-value of ttest
                    stat_struct_ttest_anova.latency_rmanova3way_DoubleStimDisplay_Choice_Stimulation.choice_x_stimulation.(fname_stat{col}){comp,7} = stats.df; % degrees of freedom
                end
                stat_struct_ttest_anova.latency_rmanova3way_DoubleStimDisplay_Choice_Stimulation.choice_x_stimulation.(fname_stat{col}) = ...
                    [{'target position','stimulation','compared to stimulation','significant','p','t','df'};...
                    stat_struct_ttest_anova.latency_rmanova3way_DoubleStimDisplay_Choice_Stimulation.choice_x_stimulation.(fname_stat{col})];
                
            end
        end
        
        %%% PERFORM ANOVA %%%
        if perform_ANOVA
            dep_var = fieldnames(stat_data);
            for dv = 1:numel(fieldnames(stat_data))
                y = stat_data.(dep_var{dv}){1,ses}.values;
                F1 = stat_data.(dep_var{dv}){1,ses}.tr_con; % factor 1, i.e. trial condition
                F2 = stat_data.(dep_var{dv}){1,ses}.choice; % factor 2, i.e. choice
                F3 = stat_data.(dep_var{dv}){1,ses}.stimulation; % factor 3, i.e. stimulation
                varnames = {'Stimulus' 'Choice' 'Stimulation'};
                if ~isnan(y)
                    [p,tbl] = anovan(y,{F1,F2,F3},'varnames',varnames,'model','full','display','off');
                    stat_results_anova.lat_2stimuli_ses_by_ses.([dep_var{dv},'_tables']){1,ses} = tbl;
                    stat_results_anova.lat_2stimuli_ses_by_ses.([dep_var{dv},'_p_values']).(['Session_',num2str(ses)]) = p;
                end
                
                effects = {'Stimulus';'Choice';'Stimulation';'Stimulus*Choice';'Stimulus*Stimulation';'Choice*Stimulation';'Stimulus*Choice*Stimulation'};
            end
        end
    end
    
    %PLOT
    condition_subfields = {'minus_80','Go','plus_80'};
    condition_titles = {'-80','Go','+80'};
    indexes_to_plot = 1:3;
    all_subplots = [2,3];
    xlim = 0; ylim = 0; ytick = 0; text_option = 0; fontsizes = 0;
    for tr_con = 1:num_tr_con2stat
        tr_conditions = fieldnames(lat_ses_by_ses2plot);
        subplot_assignment = tr_con;
        subplot_number = significance_plots_audv_ls(lat_ses_by_ses2plot.(tr_conditions(tr_con)), condition_subfields, condition_titles, indexes_to_plot, all_subplots,subplot_assignment, xlim, ylim, ytick, text_option, fontsizes);
    end
end

%%

% ------------------
%%% LATENCY DIFFERENCES - COMPARING CONTRA AND IPSI CHOICES PER STIMULATION CONDITION - WILCOXON %%%
% ------------------

if Settings.latency_differences_wilcoxon
    conditions = fieldnames(data_struct_per_session);
    stim_con_no_ctrl = stim_con(2:4); 
    clear Latency_differences; clear stat_data
    for tr_con = 1:numel(conditions)
        positions = fieldnames(data_struct_per_session.(conditions{tr_con}));
        for pos = 1:numel(positions)
            if strcmp(conditions{tr_con},'targ_distr') && strcmp(positions{pos},'combined') || ...
                    strcmp(conditions{tr_con},'single_targ') && strcmp(positions{pos},'combined') || ...
                    strcmp(conditions{tr_con},'single_distr') && strcmp(positions{pos},'combined')
                continue
            else
                choices = {'contra','ipsi'};
                for ch = 1:numel(choices)
                    if strcmp(conditions{tr_con},'single_targ') && ~isequal(positions{pos},choices{ch}) || ... % calculate latencies only for choices towards the target/distr position if no alternative is there
                            strcmp(conditions{tr_con},'single_distr') && ~isequal(positions{pos},choices{ch})
                        continue
                    else
                        fname = ['mean_latency_',choices{ch},'_choice_per_session'];
                        for stim = 1:numel(stim_con_no_ctrl)
                            for col = 1:size(data_struct_per_session.(conditions{tr_con}).(positions{pos}).stim_off,1);
                                for ses = 1:size(data_struct_per_session.(conditions{tr_con}).(positions{pos}).stim_off,2);
                                    tmp_data = data_struct_per_session.(conditions{tr_con}).(positions{pos}).(stim_con_no_ctrl{stim}){col,ses}.(fname);
                                    tmp_data_control = data_struct_per_session.(conditions{tr_con}).(positions{pos}).stim_off{col,ses}.(fname);
                                    if strcmp(conditions{tr_con},'single_targ') || strcmp(conditions{tr_con},'single_distr')
                                        Latency_differences.(conditions{tr_con}).combined.(stim_con_no_ctrl{stim})(col).(choices{ch})(ses) = tmp_data - tmp_data_control;
                                    elseif strcmp(conditions{tr_con},'targ_distr')
                                        if strcmp(positions{pos},choices{ch}) % hits
                                            Latency_differences.(conditions{tr_con}).hit_saccades.(stim_con_no_ctrl{stim})(col).(choices{ch})(ses) = tmp_data - tmp_data_control;
                                        else % errors
                                            Latency_differences.(conditions{tr_con}).error_saccades.(stim_con_no_ctrl{stim})(col).(choices{ch})(ses) = tmp_data - tmp_data_control;
                                        end
                                    else
                                        Latency_differences.(conditions{tr_con}).(positions{pos}).(stim_con_no_ctrl{stim})(col).(choices{ch})(ses) = tmp_data - tmp_data_control;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        
        %%% PERFORM WILCOXON %%%
        positions = fieldnames(Latency_differences.(conditions{tr_con}));
        for pos = 1:numel(positions)
            if strcmp(conditions{tr_con},'single_targ') || strcmp(conditions{tr_con},'single_distr')
                position_fname = 'combined';
            else
                position_fname = positions{pos};
            end
            for stim = 1:numel(stim_con_no_ctrl)
                for col = 1:size(Latency_differences.(conditions{tr_con}).(positions{pos}).(stim_con_no_ctrl{stim}),2);
                    choices = fieldnames(Latency_differences.(conditions{tr_con}).(positions{pos}).(stim_con_no_ctrl{stim})(col));
                    for ch = 1:numel(choices);
                        stat_data(:,ch) = Latency_differences.(conditions{tr_con}).(positions{pos}).(stim_con_no_ctrl{stim})(col).(choices{ch});
                    end
                    if all(isnan(stat_data(:,1))) || all(isnan(stat_data(:,2)))
                        p = NaN;
                    else
                        [p] = signrank(stat_data(:,1),stat_data(:,2));
                    end
                    stat_results_anova.latency_differences.(conditions{tr_con}).(position_fname)(col).(stim_con_no_ctrl{stim}) = p;
                    
                end
            end
            %%% FDR CORRECTION %%%
            for col = 1:size(stat_results_anova.latency_differences.(conditions{tr_con}).(position_fname),1)
                for stim = 1:numel(stim_con_no_ctrl)
                    p2correct(stim) = stat_results_anova.latency_differences.(conditions{tr_con}).(position_fname)(col).(stim_con_no_ctrl{stim}); % get relevant p-values for correction
                end
                [H_adj, crit_p, adj_ci_cvrg, P_adj]=fdr_bh(p2correct,Settings.alpha,'pdep','no');
                for stim = 1:numel(stim_con_no_ctrl)
                    stat_results_anova.latency_differences.(conditions{tr_con}).(position_fname)(col).(stim_con_no_ctrl{stim}) = P_adj(stim); % change p value
                end
            end
        end
    end
end


end
%%

% ------------------
%%% SUBFUNCTIONS %%%
% ------------------

function stat_struct_ttest = get_stat_struct_ttest (mode, stat_struct_ttest, stim, tr_con2stat, position, fname_stat, col, tmp_data, stim_con)
switch mode 
    case 'all_comparisons'
        if stim == 1
            % stat_data_ttest and stat_result_ttest organization of 1st cell level: rows: distractor colors
            % stat_results_ttest organization of the 2nd cell level: row 1: contains vectors to compare (rowwise)
            % row 2: contains stimulation conditions in same organization (1st and 2nd column), if significantly different from eacht other (3rd column) and p-value (4th column)
            stat_struct_ttest.(tr_con2stat).(position).(fname_stat){col,1}{1,1} = tmp_data;
            stat_struct_ttest.(tr_con2stat).(position).(fname_stat){col,1}{1,1+2} = stim_con{stim};
            stat_struct_ttest.(tr_con2stat).(position).(fname_stat){col,1}{2,1} = tmp_data;
            stat_struct_ttest.(tr_con2stat).(position).(fname_stat){col,1}{2,1+2} = stim_con{stim};
            stat_struct_ttest.(tr_con2stat).(position).(fname_stat){col,1}{3,1} = tmp_data;
            stat_struct_ttest.(tr_con2stat).(position).(fname_stat){col,1}{3,1+2} = stim_con{stim};
        elseif stim == 2
            stat_struct_ttest.(tr_con2stat).(position).(fname_stat){col,1}{1,2} = tmp_data;
            stat_struct_ttest.(tr_con2stat).(position).(fname_stat){col,1}{1,2+2} = stim_con{stim};
            stat_struct_ttest.(tr_con2stat).(position).(fname_stat){col,1}{4,1} = tmp_data;
            stat_struct_ttest.(tr_con2stat).(position).(fname_stat){col,1}{4,1+2} = stim_con{stim};
            stat_struct_ttest.(tr_con2stat).(position).(fname_stat){col,1}{5,1} = tmp_data;
            stat_struct_ttest.(tr_con2stat).(position).(fname_stat){col,1}{5,1+2} = stim_con{stim};
        elseif stim == 3
            stat_struct_ttest.(tr_con2stat).(position).(fname_stat){col,1}{2,2} = tmp_data;
            stat_struct_ttest.(tr_con2stat).(position).(fname_stat){col,1}{2,2+2} = stim_con{stim};
            stat_struct_ttest.(tr_con2stat).(position).(fname_stat){col,1}{4,2} = tmp_data;
            stat_struct_ttest.(tr_con2stat).(position).(fname_stat){col,1}{4,2+2} = stim_con{stim};
            stat_struct_ttest.(tr_con2stat).(position).(fname_stat){col,1}{6,1} = tmp_data;
            stat_struct_ttest.(tr_con2stat).(position).(fname_stat){col,1}{6,1+2} = stim_con{stim};
        elseif stim == 4
            stat_struct_ttest.(tr_con2stat).(position).(fname_stat){col,1}{3,2} = tmp_data;
            stat_struct_ttest.(tr_con2stat).(position).(fname_stat){col,1}{3,2+2} = stim_con{stim};
            stat_struct_ttest.(tr_con2stat).(position).(fname_stat){col,1}{5,2} = tmp_data;
            stat_struct_ttest.(tr_con2stat).(position).(fname_stat){col,1}{5,2+2} = stim_con{stim};
            stat_struct_ttest.(tr_con2stat).(position).(fname_stat){col,1}{6,2} = tmp_data;
            stat_struct_ttest.(tr_con2stat).(position).(fname_stat){col,1}{6,2+2} = stim_con{stim};
        end
    case 'compare_to_control'
        if stim == 1
            % stat_data_ttest and stat_result_ttest organization of 1st cell level: rows: distractor colors
            % stat_results_ttest organization of the 2nd cell level: row 1: contains vectors to compare (rowwise)
            % row 2: contains stimulation conditions in same organization (1st and 2nd column), if significantly different from eacht other (3rd column) and p-value (4th column)
            stat_struct_ttest.(tr_con2stat).(position).(fname_stat){col,1}{1,1} = tmp_data;
            stat_struct_ttest.(tr_con2stat).(position).(fname_stat){col,1}{1,1+2} = stim_con{stim};
            %stat_struct_ttest.(tr_con2stat).(position).(fname_stat){col,1}{2,1} = tmp_data;
            %stat_struct_ttest.(tr_con2stat).(position).(fname_stat){col,1}{2,1+2} = stim_con{stim};
            %stat_struct_ttest.(tr_con2stat).(position).(fname_stat){col,1}{3,1} = tmp_data;
            %stat_struct_ttest.(tr_con2stat).(position).(fname_stat){col,1}{3,1+2} = stim_con{stim};
        elseif stim == 2
            stat_struct_ttest.(tr_con2stat).(position).(fname_stat){col,1}{1,2} = tmp_data;
            stat_struct_ttest.(tr_con2stat).(position).(fname_stat){col,1}{1,2+2} = stim_con{stim};
        elseif stim == 3
            stat_struct_ttest.(tr_con2stat).(position).(fname_stat){col,1}{2,2} = tmp_data;
            stat_struct_ttest.(tr_con2stat).(position).(fname_stat){col,1}{2,2+2} = stim_con{stim};
        elseif stim == 4
            stat_struct_ttest.(tr_con2stat).(position).(fname_stat){col,1}{3,2} = tmp_data;
            stat_struct_ttest.(tr_con2stat).(position).(fname_stat){col,1}{3,2+2} = stim_con{stim};
        end
        
end
end
