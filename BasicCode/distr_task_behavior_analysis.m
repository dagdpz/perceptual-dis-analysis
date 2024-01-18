function [data_struct_per_run, data_struct_per_session, data_struct_across_sessions, plot_struct_double, plot_struct_single, errors, stat_results_anova, stat_struct_ttest, stat_struct_ttest_anova, stat_struct_chi2] = ...
    distr_task_behavior_analysis(Settings)

% DISTR_TASK_BEHAVIOR_ANALYSIS This function analyzes several aspects of
%   saccadic performance in a target-distractor task with color defining
%   the target and the distractors:
%       - Hitrate
%       - Contra- vs ipsiversive choices
%       - Latency
% Input: Call this function via the script 'distr_task_behavior_callscript'
%
% General analysis: Data is always analyzed per session. Additionall y or
%   alternatively it can be analyzed across sessions
%   (Settings.analyze_across_sessions = 1).
%
% Hitrate analysis: The hitrate can be calculated from successful trials
%   (monkey chose and held the correct target and was rewarded) and from
%   trials where the correct target was selected but not necessarily hold
%   for the trial to be successfully ended.
%   The hitrate is calculated per run and / or per session.
% - When hitrate is plotted per session, it is calculated per run and the
%   mean hitrate across runs in one session is displayed. Error bars
%   indicate the SEM of hitrate across runs in one session.
% - When hitrate is plotted across sessions, it is calculated per session,
%   i.e. number of all hits of all runs of one session are divided by
%   number of all trials (where a target was selected). The mean hitrate
%   across sessions is then calculated from the single hitrates per
%   session. Error bars indicate the SEM of hitrate across sessions.
% - In any case, the fitted function is calculated directly from the number
%   of hits and the number of trials of the session/ of all sessions and
%   not on any mean values.
% - PSF: A psychometric function can be fitted to the hitrate. If this
%   setting is turned on, the hitrate plots will not show mean hitrates
%   across runs / sessions but one hitrate value calculated per session /
%   per all sessions (and distractor color). Thus, no SEM can be calculated
%   and no SEM will be shown in these plots. The PSF will only be plotted
%   if 5 distractor colors are to be plotted in a plot, i.e. are present in
%   a session or across sessions.
%
% Choices analysis: The proportion of choices is calculated equally as the
%   hitrate. Session wise: Proportions are calculated from all respective
%   trials run in that run, mean and SEM of proportion is calculated from
%   these proportion values per run. Across sessions: Proportions are
%   calculated from all respective trials in that session; Mean and SEM are
%   calculated from these values.
%
%
% Latency analysis:
% - 'Mean of Means' option
%   The average of the mean latencies per run in a plot per session (i.e.
%   mean across runs) and of mean latencies per session in a plot across
%   sessions (i.e. mean across sessions) is plotted. The SEM is the SEM
%   across all means per run (per session) in a plot per session (plot
%   across sessions).
% - 'True mean of concatenated trials' option
%   Latency values per respective trial are concatenated
%   and analyzed session wise or across sessions, i.e. latency values of all
%   respective trials per session or across all sessions are concatenated,
%   mean and SEM are calculated from these values.
%
% Plots
% Data per sessions and across sessions can be plotted. If data across
%   sessions is plotted Settings.analyze_across_sessions is automatically
%   turned on.
% - Line plots - plot choice proportions, hitrate and latency on a scale of
%   distractor color ratio (G/R ratio, ratio of green to red color channel
%   value; targ-targ and single-targ trials are 0)
% - Curve fitting - fits are only done if all 5 distractor colors are
%   present in a session (or in several sessions when plotting across
%   sessions)
% - Bar plots - plot choice proportions (= hitrate in single stimulus trial
%   conditions) or latencies together with stimulation condition. Separate
%   plots for each trial condition. In  single stimulus trials bars are
%   grouped by the position of the target / distractor in single-targ /
%   single-distr trials, respectively and the hitrate is plotted. In double
%   stimulus trials the x-axis shows the position of the chosen stimulus
%   (being either selected_StimulusRight to stimulation side, ipsiversive or the
%   central fixation spot). By selecting analyze_choice, -_hitrate or
%   -_latency you can selectively plot only choice, hitrate or latency.
%   If only latency shall be plotted, both analyze_choice and
%   analyze_hitrate need to be turned off.
% - Reaction Time Histograms - Are plotted in the same organization as the
%   bar plots, but separately for each choice side (contra vs. ipsi) and
%   stimulation condition.
%
% Saving plots
% - Plots are saved to a predefined directory as .png and/or .ai file.
%   A separate directory can be defined for each dataset. If a
%   session is included in one of the datasets the plots of this session
%   will be saved to the respective directory of this dataset.
%
% Statistics
% Statistics across sessions
% - 3 way repeated measures ANOVA with stimulus type (contra-targ
%   ipsi-distr, contra-distr ipsi-targ, contra-distr ipsi-distr),
%   stimulation (off, -80, Go +80) and difficulty (easy and
%   difficult distr color) as factors. Although it calculates the
%   statistics across sessions, the input for this test are the respective
%   means per session, so Settings.analyze_across_sessions does not need to
%   be on.
%
% Output
% - Data_struct_per_run -_per_session and -_across_sessions contains the
%   respective data which is used to plot the line plots. Structure
%   fieldnames are self-explaining.
% - The bar plot structures plot_struct_double and plot_struct_single
%   contain all relevant data for the bar plots and the Reaction Time
%   histograms. _double and _single refer to the presentation of two
%   stimuli (targ-targ, distr-distr or targ-distr trials) or one single
%   stimulus (single-targ, single-distr). Plots of double stimulus trials
%   show the position of the chosen stimulus on the x-axis while plots of
%   the single-stimulus trials show the position of the presented
%   targ/distr on the x-axis and the hitrate on the y-axis.
%   The structure organization is as follows:
%   1st level: Rows: Distractor color
%              Columns: Trial condition (targ-distr consists of 2 cond.)
%   2nd level: 1st row: Plots per session (if true in Settings)
%              Columns: Sessions
%              2nd row: Plots across all sessions (if true in Settings)
%   3rd level: see structure fieldname
%   4th level: Rows: Choice position / position of presented target/distr.
%                    as described in .groupnames
%              Columns: Stimulation condition as described in .legend
% - additional output

% add to code:
% - if 'combined' position is turned off, combined position should still be
%   plotted in conditions where contra and ipsi cannot be plotted (i.e.
%   targ-targ trials
% - maybe: if no condition for one plot is met it is still plotted: e.g.
%   20160720 single targ - single distr -> don't plot these
% - !!! latency of incorrect contra and ipsi choices don't seem to be
%   plotted properly; try 20160720 or across all sessions for the curve
%   fitting (especially colors of lines, two lines should not be the same,
%   legend)


%%% CHECK INPUT SETTINGS
% turn plot_across_sessions off if only one session is to be analyzed
num_sessions = numel(Settings.folders);
if num_sessions == 1
    Settings.plot_across_sessions = 0;
end
% turn on analyze_across_sessions if plot_across_sessions is turned on in Settings
if Settings.plot_across_sessions
    Settings.analyze_across_sessions = 1;
end
% turn latency calculation on if an ANOVA is wanted but was forgotten to be turned on in Settings
error_statistics = [Settings.error_rate_rmanova2way_DistrDistrTargTarg_Stimulation,Settings.error_rate_rmanova3way_TargDistr_Difficulty_Stimulation,Settings.error_rate_rmanova3way_SingleStimDisplay_Side_Stimulation];
if Settings.statistics_on && any(error_statistics)
    Settings.analyze_error_rates = 1;
end

latency_statistics = [Settings.latency_rmanova2way_SingleStimDisplay_Side_Stimulation,Settings.latency_rmanova3way_DoubleStimDisplay_Choice_Stimulation];
if Settings.statistics_on && any(latency_statistics)
    Settings.analyze_latency = 1;
end

tr_con2ana = [];
for con = 1:numel(Settings.analyze_trial_conditions)
    tr_con2ana = horzcat(tr_con2ana,Settings.analyze_trial_conditions{con});
end
tr_con2ana = unique(tr_con2ana);
tr_con2ana = sortifmember(tr_con2ana,{'targ_targ_1HF','distr_distr_1HF','targ_targ_2HF','distr_distr_2HF','targ_distr_1HF', 'targ_distr_2HF','single_targ','single_distr'});
num_tr_con2ana = numel(tr_con2ana);
if isequal(Settings.analyze_distr_colors,'any')
    num_distr_colors2ana = 5; %Color of distractors is hard coded
else
    num_distr_colors2ana = size(Settings.analyze_distr_colors,1);
end
if isempty(Settings.analyze_by_position)
    pos2ana = {'combined'};
else
    pos2ana = Settings.analyze_by_position;
end
num_pos2ana = numel(pos2ana);

data_struct_per_session = []; data_struct_across_sessions = [];
stat_struct_ttest_anova = []; stat_struct_chi2 = [];
stat_results_anova = []; stat_struct_ttest = [];
errors.by_position = []; errors.by_condition = [];

%%% LOAD DATA FROM DATASET %%%
if ~isempty(Settings.load_data)
    load(Settings.load_data)
else
    
    %%% CALCULATE DATA FROM RUN FILES %%%
    for fol = 1:num_sessions
        disp(['Analyze session: ',Settings.folders{fol}]);
        complete_dir = char(fullfile(Settings.directory,Settings.folders(fol)));
        files = what(complete_dir);
        
        Runs = [];
        for i_Run = 1: size(files.mat,1)
            Name =   strsplit(char(files.mat(i_Run)),'_');
            Name =  char(Name(2));
            Runs{i_Run} = Name(1:end-4);
        end
        num_runs = size(files.mat,1);

        if ~isequal(num_runs, numel(Runs))
            disp('Number of Runs are not consistent')
        end
        
        delete_special_runs = false;
        data_struct_per_run = [];
        if ismember(Settings.folders{fol},Settings.exclude_runs_of_sessions(:,1),'rows') % check if this session is a session where some runs should be excluded
            [~,del_idx] = ismember(Settings.folders{fol},Settings.exclude_runs_of_sessions(:,1));
            delete_special_runs = true;
        end
        for file = 1: num_runs %change to 1 %What happens if we miss a run??
            file_to_load = char(fullfile(complete_dir,files.mat(file)));
            load(file_to_load);
            
         
            %%
            num_trials = numel(trial);
            
            if delete_special_runs && ismember(num2str(file),Settings.exclude_runs_of_sessions(del_idx,2:end)) % if some runs should be deleted and if the current run is one of them
                continue;
            end
               %% change  trial(i).microstim = 1 if the Run is an PostInjection Run
            % the change is based on the given Variable: Settings.PreInjection_Runs
            if  strcmp(Settings.Experiment , 'Inactivation')
                Name =  strsplit(char(files.mat(file)),'_');
                Name =  char(Name(2));
                
                if ~any(ismember(Settings.PreInjection_Runs{fol},str2num(Name(1:end-4))) | isnan(Settings.PreInjection_Runs{fol})) % the loaded run belongs NOT to PreInjection run
                    disp(['changed to PostInjection  ', Name(1:end-4)])
                    for i = 1:length([trial.microstim])
                        trial(i).microstim = 1;
                    end
                else
                   disp(['PreInjection  ', Name(1:end-4)])
                    for i = 1:length([trial.microstim])
                        trial(i).microstim = 0;
                    end
                    
                end
            end
            if isequal([trial.type],Settings.task_type*ones(1,num_trials)) % check type
                if num_trials >= Settings.min_num_trials_per_run % check run size                    
                    % distractor colors per trial
                    all_colors_per_trial = arrayfun(@(x) vertcat(x.task.eye.tar.color_dim),trial,'uni',0);
                    distr_color_per_trial = cellfun(@(x) x(~ismember(x,Settings.fixation_color_dim,'rows') & ~ismember(x,Settings.target_color_dim,'rows'),:),all_colors_per_trial,'UniformOutput',0);
                    
                    red_targets_per_trial = cellfun(@(x) find(ismember(x,Settings.target_color_dim,'rows')),all_colors_per_trial,'UniformOutput',0);
                    distr_per_trial = cellfun(@(x) find(~ismember(x,Settings.fixation_color_dim,'rows') & ~ismember(x,Settings.target_color_dim,'rows')),all_colors_per_trial,'UniformOutput',0);
                    %which target was the fixation target?
                    fix_per_trial = cellfun(@(x) find(ismember(x,Settings.fixation_color_dim,'rows')),all_colors_per_trial,'UniformOutput',0);
                    red_targets_per_trial_empty = cellfun(@(x) isempty(x),red_targets_per_trial,'UniformOutput',1);
                    
                    Idx_two_red_targets_per_trial                           = cellfun(@(x) isequal(size(x,1),2),red_targets_per_trial,'UniformOutput',1);
                    Idx_two_distr_per_trial                                 = cellfun(@(x) isequal(size(x,1),2),distr_per_trial,'UniformOutput',1);
                    correct_target_per_trial                                = red_targets_per_trial;
                    correct_target_per_trial(red_targets_per_trial_empty)   = fix_per_trial(red_targets_per_trial_empty);
                    
                    red_targ_or_distr_per_trial = red_targets_per_trial;
                    red_targ_or_distr_per_trial(red_targets_per_trial_empty) = distr_per_trial(red_targets_per_trial_empty);
                    non_fix_stimuli_per_trial = red_targ_or_distr_per_trial; % also includes targ_targ and distr_distr trials
                    red_targ_or_distr_per_trial(Idx_two_red_targets_per_trial) = {[]};
                    red_targ_or_distr_per_trial(Idx_two_distr_per_trial) = {[]}; % contains indices of the target (if present and only if it is a single target or a target-distractor trial)
                   
                    %% DISTRACTOR COLORS
                    % in each trial. If no target is present it contains the index of the distractor (only if it is a single
                    % distractor trial).
                    distr_colors_all = [NaN,NaN,NaN];
                    idx = 1;
                    distr_ID = 0;
                    distractors = nan(num_distr_colors2ana,num_trials);
                    for trial_no = 1:num_trials
                        curr_distr_colors = distr_color_per_trial{trial_no};
                        if size(distr_color_per_trial{trial_no},1) > 1
                            if ~isequal(curr_distr_colors(1,:),curr_distr_colors(2,:));
                                warning(char(strcat(sprintf('There were two distractors with different colors in trial %d in',trial_no),{' '},file_to_load,{'. The one with index 1 was used for analysis!'})));
                            end
                            curr_distr_colors = curr_distr_colors(1,:);
                        end
                        % This will overwrite logicals in "distractors" with the
                        % last distractor color per trial in "distr_color_per_trial"
                        
                        if ~ismember(curr_distr_colors,distr_colors_all,'rows') % "distr_colors_all" is the sequence of appearance of the distractor colors in one run
                            distr_colors_all(idx,:) = curr_distr_colors;
                            idx = idx + 1;
                            distr_ID = distr_ID + 1;
                            distractors(distr_ID,trial_no) = trial_no; % "distractors": rows represent distractor colors in order of appearance as in "distr_colors_all"
                        else
                            [~,distr_ID_member] = ismember(curr_distr_colors,distr_colors_all,'rows');
                            distractors(distr_ID_member,trial_no) = trial_no;
                        end
                    end;
                    %                 distractors = distractors(~isnan(distractors));
                    
                    distr_colors2ana = distr_colors_all;
                    num_distr_colors2ana = size(distr_colors2ana,1);
                    
                    % calculate G/R ratio of RGB values
                    G_R_ratio = nan(num_distr_colors2ana,1);
                    for col = 1:num_distr_colors2ana
                        G_R_ratio(col,1) = distr_colors2ana(col,2) / distr_colors2ana(col,1);
                    end
                    
                    % sort data according to G/R ratio: from yellow (small ratio) to red (high ratio)
                    [G_R_ratio_sorted,ind] = sort(G_R_ratio);
                    distr_colors2ana_sorted = distr_colors2ana(ind,:);
                    distractor_color = distractors(ind,:); % indices for each distractor color(rows) in each trial (columns)
                    
                    % Get logical indices of trial condition for each trial


                    num_targets_per_trial = cellfun(@(x) sum(ismember(x,Settings.target_color_dim,'rows')),all_colors_per_trial,'UniformOutput',1);
                    num_distr_per_trial = cellfun(@(x) sum(~ismember(x,Settings.fixation_color_dim,'rows') & ~ismember(x,Settings.target_color_dim,'rows')),all_colors_per_trial,'UniformOutput',1);                   
                    
                    x_pos1HF_stimuli = NaN(1,num_trials) ; 
                    for trial_no = 1:num_trials
                        if trial(trial_no).task.eye.tar(1).x  == trial(trial_no).task.eye.tar(2).x
                            x_pos1HF_stimuli(trial_no) = 1;
                        elseif trial(trial_no).task.eye.tar(1).x  == 0 || trial(trial_no).task.eye.tar(2).x  == 0
                            x_pos1HF_stimuli(trial_no) = 2 ;
                        elseif trial(trial_no).task.eye.tar(1).x  ~= trial(trial_no).task.eye.tar(2).x
                            x_pos1HF_stimuli(trial_no) = 3 ;
                        else
                            x_pos1HF_stimuli(trial_no) = 0 ;
                        end
                    end
                    trial_cond.targ_targ_1HF    = find(num_targets_per_trial == 2 & x_pos1HF_stimuli == 1 ); % target-target trials
                    trial_cond.targ_targ_2HF    = find(num_targets_per_trial == 2 & x_pos1HF_stimuli == 3 ); % target-target trials
                   
                    trial_cond.distr_distr_1HF  = find(num_distr_per_trial == 2 & x_pos1HF_stimuli == 1 ); % distr-distr trials
                    trial_cond.distr_distr_2HF  = find(num_distr_per_trial == 2 & x_pos1HF_stimuli == 3 ); % distr-distr trials
                   
                    trial_cond.single_targ  = find(num_targets_per_trial == 1 & num_distr_per_trial == 0); % target only trials
                    trial_cond.single_distr = find(num_targets_per_trial == 0 & num_distr_per_trial == 1); % distr only trials

                    %%How to differentiate where the stimuli are
                    trial_cond.targ_distr_1HF   = find(num_targets_per_trial == 1 & num_distr_per_trial == 1 & x_pos1HF_stimuli == 1 ); % target-distr trials
                    trial_cond.targ_distr_2HF   = find(num_targets_per_trial == 1 & num_distr_per_trial == 1 & x_pos1HF_stimuli == 3 ); % target-distr trials

                    if ismember(trial_cond.targ_distr_1HF,trial_cond.targ_distr_2HF), disp('Problem in separating stimulus conditions'),end
                    if ismember(trial_cond.targ_distr_1HF,trial_cond.targ_targ_1HF), disp('Problem in separating stimulus conditions'),end



                    % Count trials per condition to check if all necessary trial conditions are present in this run, the actual CHECK is below within the BIG LOOP
                    num_trials_per_cond = zeros(size(Settings.necessary_trial_conditions,2),1);
                    for nec_con = 1:numel(Settings.necessary_trial_conditions)
                        num_trials_per_cond(nec_con,1) = numel(trial_cond.(Settings.necessary_trial_conditions{nec_con})); % number of trials per necessary trial condition as defined in 'Settings.necessary_trial_conditions'
                    end
                    
                    % Get indices of positions of target or distr for each trial
                    x_pos_correct       = nan(1,num_trials); % will not contain distr-distr and targ-targ trials since the position is irrelevant for these
                    y_pos_correct       = nan(1,num_trials);
                    comb_pos_correct    = nan(1,num_trials);
                    contra_pos_correct  = nan(1,num_trials);
                    ipsi_pos_correct    = nan(1,num_trials);
                    x_pos_selected      = nan(1,num_trials);
                    contra_pos_selected = nan(1,num_trials);
                    ipsi_pos_selected   = nan(1,num_trials);
                    fix_pos_selected    = nan(1,num_trials);
                    y_pos_selected      = nan(1,num_trials);
                    contra_up_pos_selected      = nan(1,num_trials);
                    contra_down_pos_selected    = nan(1,num_trials);
                    ipsi_up_pos_selected        = nan(1,num_trials);
                    ipsi_down_pos_selected      = nan(1,num_trials);
                    contra_up_pos_correct       = nan(1,num_trials);
                    contra_down_pos_correct     = nan(1,num_trials);
                    ipsi_up_pos_correct         = nan(1,num_trials);
                    ipsi_down_pos_correct       = nan(1,num_trials);
                    
                    
                    for trial_no = 1:num_trials
                        % sort trials by position of correct target
                        if ismember(trial_no,trial_cond.targ_distr_1HF) || ...  if this trial is a targ-distr, single targ or single distr trial
                                ismember(trial_no,trial_cond.targ_distr_2HF) || ... 
                                ismember(trial_no,trial_cond.single_targ) || ...
                                ismember(trial_no,trial_cond.targ_targ_1HF) || ... 
                                ismember(trial_no,trial_cond.distr_distr_1HF) || ... 
                                ismember(trial_no,trial_cond.single_distr)
                            
                            %get the position of the stimuli - !But 
                           if  ismember(trial_no,trial_cond.targ_targ_1HF)
                            x_pos_correct(trial_no) = trial(trial_no).task.eye.tar(Idx_two_red_targets_per_trial(trial_no)).x; % x position values of targ or distr per trial
                            y_pos_correct(trial_no) = trial(trial_no).task.eye.tar(Idx_two_red_targets_per_trial(trial_no)).y; % y position values of targ or distr per trial 
                           
                           elseif ismember(trial_no,trial_cond.distr_distr_1HF)
                            x_pos_correct(trial_no) = trial(trial_no).task.eye.tar(Idx_two_distr_per_trial(trial_no)).x; % x position values of targ or distr per trial
                            y_pos_correct(trial_no) = trial(trial_no).task.eye.tar(Idx_two_distr_per_trial(trial_no)).y; % y position values of targ or distr per trial 
                           else %single simuli & distr-targ stimuli: get the position (ipsi vs contra) of the target 
                           x_pos_correct(trial_no) = trial(trial_no).task.eye.tar(red_targ_or_distr_per_trial{trial_no}).x; % x position values of targ or distr per trial
                           y_pos_correct(trial_no) = trial(trial_no).task.eye.tar(red_targ_or_distr_per_trial{trial_no}).y; % y position values of targ or distr per trial
                           end
                           %%%!!! KK:the offset_x is here important
                            if strcmp('right',Settings.hemisphere_of_stimulation) % if stimulation is in the right hemisphere
                                if x_pos_correct(trial_no) <  trial(trial_no).task.eye.fix.x %target is one the left side of the screen
                                    contra_pos_correct(trial_no) = trial_no;
                                    if ismember(trial_no,trial_cond.targ_distr_1HF)
%                                         if  y_pos_correct(trial_no) < trial(trial_no).task.eye.fix.y
%                                             contra_down_pos_correct(trial_no) = trial_no;
%                                         elseif y_pos_correct(trial_no) > trial(trial_no).task.eye.fix.y
%                                             contra_up_pos_correct(trial_no) = trial_no;
%                                         end
                                    end
                              
                                elseif x_pos_correct(trial_no) > trial(trial_no).task.eye.fix.x
                                    ipsi_pos_correct(trial_no) = trial_no;
                                     if ismember(trial_no,trial_cond.targ_distr_1HF)
%                                         if  y_pos_correct(trial_no) < trial(trial_no).task.eye.fix.y
%                                             ipsi_down_pos_correct(trial_no) = trial_no;
%                                         elseif y_pos_correct(trial_no) > trial(trial_no).task.eye.fix.y
%                                             ipsi_up_pos_correct(trial_no) = trial_no;
%                                         end
                                    end
                                end
                            elseif strcmp('left',Settings.hemisphere_of_stimulation) % if stimulation is in the left hemisphere
                                if x_pos_correct(trial_no) <  trial(trial_no).task.eye.fix.x %target is one the left side of the screen
                                    ipsi_pos_correct(trial_no) = trial_no;
                                    
                                elseif x_pos_correct(trial_no) > trial(trial_no).task.eye.fix.x
                                    contra_pos_correct(trial_no) = trial_no;
                                end
                                
                            end
                            
                            comb_pos_correct(trial_no) = trial_no; %
                            %adding all trials in here?
                            
                        elseif  ismember(trial_no,trial_cond.targ_targ_2HF) || ...  if this trial is targ-targ or a distr-distr trial for both HF
                                ismember(trial_no,trial_cond.distr_distr_2HF)
                            comb_pos_correct(trial_no) = trial_no; 
                        end
                        
                        % sort trials by position of selected target for choice analysis
                        if ~isnan(trial(trial_no).target_selected(1)) % if no target was selected there is NaN
                            x_pos_selected(trial_no) = trial(trial_no).task.eye.tar(trial(trial_no).target_selected(1)).x; % x position values of selected target per trial
                            y_pos_selected(trial_no) = trial(trial_no).task.eye.tar(trial(trial_no).target_selected(1)).y; % x position values of selected target per trial
                        end
                        
                        
             %Here is defined which trials count for ipsi, contra or combined           
                        if strcmp('right',Settings.hemisphere_of_stimulation) % if stimulation is in the right hemisphere
                            %%!!! KK:the offset_x is here important
                            if x_pos_selected(trial_no) < trial(trial_no).task.eye.fix.x
                                contra_pos_selected(trial_no) = trial_no;
                              if ismember(trial_no,trial_cond.targ_distr_1HF) 
                                if  y_pos_selected(trial_no) < trial(trial_no).task.eye.fix.y
                                 contra_down_pos_selected(trial_no) = trial_no;
                                elseif y_pos_selected(trial_no) > trial(trial_no).task.eye.fix.y
                                 contra_up_pos_selected(trial_no) = trial_no;
                                end
                              end
                            elseif x_pos_selected(trial_no) > trial(trial_no).task.eye.fix.x
                                ipsi_pos_selected(trial_no) = trial_no;
                                           if ismember(trial_no,trial_cond.targ_distr_1HF) 
                                 if  y_pos_selected(trial_no) < trial(trial_no).task.eye.fix.y
                                 ipsi_down_pos_selected(trial_no) = trial_no;
                                elseif y_pos_selected(trial_no) > trial(trial_no).task.eye.fix.y
                                 ipsi_up_pos_selected(trial_no) = trial_no;
                                 end
                                           end
                            elseif x_pos_selected(trial_no) == trial(trial_no).task.eye.fix.x
                                fix_pos_selected(trial_no) = trial_no;
                            end
                            
                            
                        elseif strcmp('left',Settings.hemisphere_of_stimulation) % if stimulation is in the left hemisphere
                            if x_pos_selected(trial_no) < trial(trial_no).task.eye.fix.x %left side of the screen
                                ipsi_pos_selected(trial_no) = trial_no;
                                if  y_pos_selected(trial_no) < trial(trial_no).task.eye.fix.y
                                 ipsi_down_pos_selected(trial_no) = trial_no;
                                elseif y_pos_selected(trial_no) > trial(trial_no).task.eye.fix.y
                                 ipsi_up_pos_selected(trial_no) = trial_no;
                                end                                  
                            elseif x_pos_selected(trial_no) > trial(trial_no).task.eye.fix.x
                                contra_pos_selected(trial_no) = trial_no; %contra_pos_selected
                                if  y_pos_selected(trial_no) < trial(trial_no).task.eye.fix.y
                                 contra_down_pos_selected(trial_no) = trial_no;
                                elseif y_pos_selected(trial_no) > trial(trial_no).task.eye.fix.y
                                 contra_up_pos_selected(trial_no) = trial_no;
                                end
                            elseif x_pos_selected(trial_no) == trial(trial_no).task.eye.fix.x
                                fix_pos_selected(trial_no) = trial_no;
                            end
                        end
                    end
                    
                    %% correct trials KK changed contra ipsi combined
                    pos_correct_per_trial.contra    = contra_pos_correct(~isnan(contra_pos_correct)); % indices of correct trials with contra targ or distr
                    pos_correct_per_trial.ipsi      = ipsi_pos_correct(~isnan(ipsi_pos_correct)); % indices of trials with ipsi targ or distr

                    pos_correct_per_trial.contra_up     = contra_up_pos_correct(~isnan(contra_up_pos_correct)); % indices of correct trials with contra targ or distr
                    pos_correct_per_trial.ipsi_up       = ipsi_up_pos_correct(~isnan(ipsi_up_pos_correct)); % indices of trials with ipsi targ or distr
                    pos_correct_per_trial.contra_down   = contra_down_pos_correct(~isnan(contra_down_pos_correct)); % indices of correct trials with contra targ or distr
                    pos_correct_per_trial.ipsi_down     = ipsi_down_pos_correct(~isnan(ipsi_down_pos_correct)); % indices of trials with ipsi targ or distr
                    pos_correct_per_trial.combined      = comb_pos_correct(~isnan(comb_pos_correct)); % indices of trials with either contra or ipsi targ or distr (-> all trials)
                    % monkey selected
                    pos_selected_per_trial.contra       = contra_pos_selected(~isnan(contra_pos_selected));
                    pos_selected_per_trial.contra_up    = contra_up_pos_selected(~isnan(contra_up_pos_selected));
                    pos_selected_per_trial.contra_down  = contra_down_pos_selected(~isnan(contra_down_pos_selected));

                    pos_selected_per_trial.ipsi         = ipsi_pos_selected(~isnan(ipsi_pos_selected)); %ipsi_pos_selected
                    pos_selected_per_trial.ipsi_up      = ipsi_up_pos_selected(~isnan(ipsi_up_pos_selected)); %ipsi_pos_selected
                    pos_selected_per_trial.ipsi_down    = ipsi_down_pos_selected(~isnan(ipsi_down_pos_selected)); %ipsi_pos_selected
                    
                    pos_selected_per_trial.fixation = fix_pos_selected(~isnan(fix_pos_selected));
                    
                    % Call monkeypsych_analyze to calculate latencies
                    if Settings.analyze_latency || Settings.analyze_error_rates
                        [out_comp,variable_to_test,counter]= monkeypsych_analyze_working(...
                            {file_to_load},...% change 3 to file
                            {'summary',0,'display',0,'keep_raw_data',1,'correct_offset',0}); %, 'success',1,'n_targets', 2,
                        if Settings.analyze_latency
                            latency = [out_comp{1,1}.saccades.lat];
                        end
                        if Settings.analyze_error_rates
                            eye_traces = [out_comp{1, 1}.raw];
                        end
                    end
                    
                    %%% GET STIMULATION CONDITIONS %%%
                    if strcmp(Settings.Experiment , 'Microstimulation')
                        idx_all_trials          = [trial.n]; % numerical indeces of all trials in this run
                        trial_stim_timing       = arrayfun(@(x) horzcat(x.task.microstim.start{1}),trial,'uni',1);
                        trial_stim_on           = arrayfun(@(x) horzcat(x.task.microstim.stim_on),trial,'uni',1);
                        stim_condition.stim_off = idx_all_trials(trial_stim_on==0); % stim is off
                        if strcmp(Settings.Monkey, 'Cornelius')
                            disp('Stimulation timing- Cornelius')
                        stim_condition.minus_250ms   = idx_all_trials(trial_stim_on == 1 & trial_stim_timing == -0.25); % stim needs to be on and timing correct
                        stim_condition.go_signal    = idx_all_trials(trial_stim_on == 1 & trial_stim_timing == -0.1);
                        
                       
                        
                        if  ismember(str2num(Settings.folders{fol}),[20180123,20180122, 20180121, 20180119,20180114, 20180113,20180111,20180109])
                        unique(trial_stim_timing)
                        stim_condition.minus_50ms    = idx_all_trials(trial_stim_on == 1 & trial_stim_timing == -0.05);
                         else
                        stim_condition.minus_50ms    = idx_all_trials(trial_stim_on == 1 & trial_stim_timing == +0.05);
                        end
                       % end
                        elseif strcmp(Settings.Monkey, 'Curius')
                           disp('Stimulation timing- Curius')
                        stim_condition.minus_80ms = idx_all_trials(trial_stim_on == 1 & trial_stim_timing == -0.08); % stim needs to be on and timing correct
                        stim_condition.go_signal = idx_all_trials(trial_stim_on == 1 & trial_stim_timing == 0);
                        stim_condition.plus_80ms = idx_all_trials(trial_stim_on == 1 & trial_stim_timing == +0.08);
                        end
                        num_stim_con = numel(fieldnames(stim_condition));
                        stim_con = fieldnames(stim_condition);
                    elseif strcmp(Settings.Experiment , 'Inactivation')
                        idx_all_trials              = [trial.n]; % numerical indeces of all trials in this run
                        trial_stim_on               = arrayfun(@(x) horzcat(x.microstim),trial,'uni',1);
                        stim_condition.stim_off     = idx_all_trials(trial_stim_on==0); %  preinjection
                        stim_condition.go_signal    = idx_all_trials(trial_stim_on == 1 ); % stim needs to be on and timing correct
                        num_stim_con                = numel(fieldnames(stim_condition));
                        stim_con                    = fieldnames(stim_condition);
                        
                     elseif strcmp(Settings.Experiment , 'Electrophysiology')
                        idx_all_trials              = [trial.n]; % numerical indeces of all trials in this run
                        trial_stim_on               = arrayfun(@(x) horzcat(x.microstim),trial,'uni',1);
                        stim_condition.stim_off     = idx_all_trials(trial_stim_on==0); %  preinjection
                        stim_condition.go_signal    = idx_all_trials(trial_stim_on == 1 ); % stim needs to be on and timing correct
                        num_stim_con                = numel(fieldnames(stim_condition));
                        stim_con                    = fieldnames(stim_condition);
                        
                    end
                    %%% GET NUMBER OF ALL INITIATED TRIALS %%%
                    % all trials with abort_state = 2, i.e. abort eye fix acq
                    initiated_trials = idx_all_trials([trial.aborted_state]~=2);
                    
                    
                    %%% CALCULATE ERROR RATES %%%
                    if Settings.analyze_error_rates
                        %                         cond = fieldnames(trial_cond);
                        %                         session = Settings.folders{fol};
                        %                         run = [session,'_',num2str(file)];
                        %                         for con = 1:numel(cond)
                        %                             tr_condition = cond{con};
                        %                             [errors.by_position] = count_errors_by_position(Settings,errors.by_position,trial,non_fix_stimuli_per_trial,trial_cond.(cond{con}),tr_condition,session,run);
                        %
                        %                         end
                        fix_break_error = arrayfun(@(x) horzcat(strcmp(x.abort_code,'ABORT_EYE_FIX_HOLD_STATE')),trial,'uni',1);
                        trial_fix_break_error = idx_all_trials(fix_break_error);
                        tar_acq_error = arrayfun(@(x) horzcat(strcmp(x.abort_code,'ABORT_EYE_TAR_ACQ_STATE')),trial,'uni',1);
                        trial_tar_acq_error = idx_all_trials(tar_acq_error);
                        tar_hold_error = arrayfun(@(x) horzcat(strcmp(x.abort_code,'ABORT_EYE_TAR_HOLD_STATE')),trial,'uni',1);
                        trial_tar_hold_error = idx_all_trials(tar_hold_error);
                        if ~exist('error_data2plot', 'var') % created empty variable
                            for i = 1:numel(stim_con)
                                error_data2plot.(stim_con{i}) = struct('x',NaN,'y',NaN);
                            end
                        end
                        [errors.by_condition,error_data2plot] = count_errors_by_condition(Settings,errors.by_condition,trial,pos_correct_per_trial,trial_cond,stim_condition,...
                            trial_fix_break_error,trial_tar_hold_error,trial_tar_acq_error,distractor_color,eye_traces, error_data2plot,fol,file);
                        
                    end
                    
                    %%%% START BIG LOOP & GET DATA FOR EVERY COMBINATION OF CONDITIONS %%%%
                    warning_count = 0;
                    empty_conditions = [];
                    for tr_con = 1:num_tr_con2ana
                        for pos = 1: num_pos2ana
                            for stim = 1:num_stim_con
                                
                                for col = 1:num_distr_colors2ana
                                    if (strcmp(tr_con2ana{tr_con},'targ_targ_2HF') || strcmp(tr_con2ana{tr_con},'distr_distr_2HF') )&& (strcmp(pos2ana{pos},'contra') || strcmp(pos2ana{pos},'ipsi')|| strcmp(pos2ana{pos},'ipsi_up')|| strcmp(pos2ana{pos},'ipsi_down') || strcmp(pos2ana{pos},'contra_up')|| strcmp(pos2ana{pos},'contra_down'))
                                        continue % don't do calculations for trial conditions with targ/distr positions contra or ipsi
                                        
                                    else % no contra or ipsi in targ_targ or distr_distr
                                        if ~all(num_trials_per_cond) % CHECK if all values are not nonzero, i.e. if there is at least one necessary condition with zero trials
                                            data_struct_per_run.(tr_con2ana{tr_con}).(pos2ana{pos}).(stim_con{stim}){col,file} = [];
                                        else % if all values are nonzero, i.e. if there is at least one trial per necessary trial conditions
                                            trials2ana_tmp1.(tr_con2ana{tr_con}).(pos2ana{pos}).(stim_con{stim}){col,file} = ... % intermediate steps of trials2ana
                                                intersect(trial_cond.(tr_con2ana{tr_con}), ...
                                                pos_correct_per_trial.(pos2ana{pos}));
                                            
                                            trials2ana_tmp2.(tr_con2ana{tr_con}).(pos2ana{pos}).(stim_con{stim}){col,file} = ... % intermediate steps of trials2ana
                                                intersect(trials2ana_tmp1.(tr_con2ana{tr_con}).(pos2ana{pos}).(stim_con{stim}){col,file}, ...
                                                stim_condition.(stim_con{stim}));
                                            if strcmp(tr_con2ana{tr_con},'targ_targ_2HF') || strcmp(tr_con2ana{tr_con},'targ_targ_1HF') || strcmp(tr_con2ana{tr_con},'single_targ') % don't get intersect with any distractor color trials because there are no distractors in these trials
                                                col = 1;%only one colume because no distractor
                                                data_struct_per_run.(tr_con2ana{tr_con}).(pos2ana{pos}).(stim_con{stim}){col,file}.trials2ana = ... % put trials2ana into data_struct as if there was only 1 col (there are not distr colors here)
                                                    trials2ana_tmp2.(tr_con2ana{tr_con}).(pos2ana{pos}).(stim_con{stim}){col,file};
                                                data_struct_per_run.(tr_con2ana{tr_con}).(pos2ana{pos}).(stim_con{stim}){col,file}.G_R_ratios = 0; % set x_values for line graph plotting in the data_struct to zero
                                            else
                                                data_struct_per_run.(tr_con2ana{tr_con}).(pos2ana{pos}).(stim_con{stim}){col,file}.trials2ana = ... % put trials2ana into data_struct dependent on difficulty level
                                                    intersect(trials2ana_tmp2.(tr_con2ana{tr_con}).(pos2ana{pos}).(stim_con{stim}){col,file}, ...
                                                    distractor_color(col,:));
                                                data_struct_per_run.(tr_con2ana{tr_con}).(pos2ana{pos}).(stim_con{stim}){col,file}.G_R_ratios = G_R_ratio_sorted; % put x_values for line graph plotting into data_struct
                                            end
                                            data_struct_per_run.(tr_con2ana{tr_con}).(pos2ana{pos}).(stim_con{stim}){col,file}.num_trials2ana = ... % put num_trials2ana into data_struct
                                                size(data_struct_per_run.(tr_con2ana{tr_con}).(pos2ana{pos}).(stim_con{stim}){col,file}.trials2ana,2);
                                            trials2ana = data_struct_per_run.(tr_con2ana{tr_con}).(pos2ana{pos}).(stim_con{stim}){col,file}.trials2ana; % current trials2ana for these specific conditions
                                            num_trials2ana = data_struct_per_run.(tr_con2ana{tr_con}).(pos2ana{pos}).(stim_con{stim}){col,file}.num_trials2ana; % current num_trials2ana
                                            
                                            trial_struct2ana = trial(trials2ana);
                                            correct_target_per_trial2ana = correct_target_per_trial(trials2ana);
                                            initiated_trials2ana = intersect(initiated_trials,trials2ana); % get only the initiated trials for this condition
                                            num_initiated_trials2ana = numel(initiated_trials2ana); % relevant for error rate calculation
                                            
                                            %%% HITS AND TOTAL TRIALS  %%%
                                            % if strcmp(stim_con{stim},'stim_off');
                                            hits = nan(1,num_trials2ana);
                                            total_trials = nan(1,num_trials2ana);
                                            targ_selected = nan(1,num_trials2ana);
                                            switch Settings.hitrate_mode
                                                case 'success'
                                                    for trial_no = 1:num_trials2ana
                                                        if ~isnan(trial_struct2ana(trial_no).target_selected(1))
                                                            targ_selected(trial_no) = trial_struct2ana(trial_no).target_selected(1);
                                                            total_trials(trial_no) = trials2ana(trial_no);
                                                            if trial_struct2ana(trial_no).success
                                                                hits(trial_no) = trials2ana(trial_no);
                                                            end
                                                        end
                                                    end
                                                case 'correct target selected'
                                                    for trial_no = 1:num_trials2ana
                                                        if ~isnan(trial_struct2ana(trial_no).target_selected(1))
                                                            targ_selected(trial_no) = trial_struct2ana(trial_no).target_selected(1);
                                                            total_trials(trial_no) = trials2ana(trial_no);
                                                            if ismember(trial_struct2ana(trial_no).target_selected(1),correct_target_per_trial2ana{trial_no})
                                                                hits(trial_no) = trials2ana(trial_no);
                                                            end
                                                        end
                                                    end
                                            end
                                            % else
                                            %     hits = [];
                                            %     total_trials = [];
                                            %     targ_selected = [];
                                            % end
                                            hits = hits(~isnan(hits));
                                            total_trials = total_trials(~isnan(total_trials));
                                            num_hits = numel(hits);
                                            num_total_trials = numel(total_trials);
                                            data_struct_per_run.(tr_con2ana{tr_con}).(pos2ana{pos}).(stim_con{stim}){col,file}.run              = file;
                                            data_struct_per_run.(tr_con2ana{tr_con}).(pos2ana{pos}).(stim_con{stim}){col,file}.hits             = hits; % put hits into data_struct
                                            data_struct_per_run.(tr_con2ana{tr_con}).(pos2ana{pos}).(stim_con{stim}){col,file}.num_hits         = num_hits; % put num_hits into data_struct
                                            data_struct_per_run.(tr_con2ana{tr_con}).(pos2ana{pos}).(stim_con{stim}){col,file}.total_trials     = total_trials; % put total_trials into data_struct
                                            data_struct_per_run.(tr_con2ana{tr_con}).(pos2ana{pos}).(stim_con{stim}){col,file}.num_total_trials     = num_total_trials; % put num_total_trials into data_struct
                                            data_struct_per_run.(tr_con2ana{tr_con}).(pos2ana{pos}).(stim_con{stim}){col,file}.num_initiated_trials = num_initiated_trials2ana;
                                            %%% HITRATE %%%
                                            
                                            if Settings.analyze_hitrate
                                                if num_hits ~= 0 && num_total_trials ~= 0
                                                    data_struct_per_run.(tr_con2ana{tr_con}).(pos2ana{pos}).(stim_con{stim}){col,file}.hitrate = num_hits/num_total_trials; % put hitrate into data_struct
                                                else
                                                    data_struct_per_run.(tr_con2ana{tr_con}).(pos2ana{pos}).(stim_con{stim}){col,file}.hitrate = NaN; % no hits and/or no total_trials for these conditions
                                                    warning_count = warning_count + 1;
                                                    empty_conditions{warning_count,1} = ...
                                                        ['Trial condition: ',tr_con2ana{tr_con},' | Position: ',pos2ana{pos},' | Stim condition: ',stim_con{stim},' | color: ',num2str(col),' | run: ',num2str(file)];
                                                end
                                            end
                                            
                                            %%% CHOICES %%%
                                            
                                            if Settings.analyze_choice || Settings.analyze_latency % pos_selected_per_trials2ana is needed for latencies to contra or ipsi stimuli
                                                
                                                correct_pos = fieldnames(pos_correct_per_trial);
                                                % for targ_targ & dist_distr the selection of right & left position are correct
                                                for copos = 1:numel(correct_pos)
                                                    pos_correct_per_trials2ana.(correct_pos{copos}) = intersect(pos_correct_per_trial.(correct_pos{copos}),trials2ana); % trials2ana sorted by position of correct target/distr
                                                end
                                                
                                                choice_pos = {'contra', 'ipsi', 'fixation'} ; %fieldnames(pos_selected_per_trial);
                                                for chpos = 1:numel(choice_pos)

                                                    % which trials are in trials2ana and belong to category targ_targ & single_targ
                                                    if strcmp(tr_con2ana{tr_con},'targ_targ_2HF') || strcmp(tr_con2ana{tr_con},'single_targ') % if targ-targ or distr-distr trial
                                                        pos_selected_per_trials2ana.(choice_pos{chpos}) = intersect(pos_selected_per_trial.(choice_pos{chpos}),trials2ana); % get choice for all trials2ana
                                                    elseif  strcmp(tr_con2ana{tr_con},'targ_distr_1HF')  % if targ-targ or distr-distr trial
                                                       % pos_selected_per_trials2ana.(choice_pos{chpos}) = intersect(pos_selected_per_trial.(choice_pos{chpos}),pos_correct_per_trials2ana.(pos2ana{pos})); % get choice for all trials where the correct position is contra / ipsi / center
                                                        if strcmp (choice_pos{chpos}, 'fixation')
                                                          pos_selected_per_trials2ana.(choice_pos{chpos}) = intersect(pos_selected_per_trial.(choice_pos{chpos}),pos_correct_per_trials2ana.(pos2ana{pos})); % get choice for all trials where the correct position is contra / ipsi / center  
                                                        elseif ~strcmp (pos2ana{pos}, 'combined')  %
                                                        pos_selected_per_trials2ana.([choice_pos{chpos},'_up']) = intersect(pos_selected_per_trial.([choice_pos{chpos},'_up']),pos_correct_per_trials2ana.([pos2ana{pos},'_up'])); % get choice for all trials where the correct position is contra / ipsi / center
                                                        pos_selected_per_trials2ana.([choice_pos{chpos},'_down']) = intersect(pos_selected_per_trial.([choice_pos{chpos},'_down']),pos_correct_per_trials2ana.([pos2ana{pos},'_down'])); % get choice for all trials where the correct position is contra / ipsi / center
                                                        pos_selected_per_trials2ana.(choice_pos{chpos}) =   [pos_selected_per_trials2ana.([choice_pos{chpos},'_up']) , pos_selected_per_trials2ana.([choice_pos{chpos},'_down'])]; 
                                                        end
                                                    else % if targ-distr, single-targ or single-distr
                                                        pos_selected_per_trials2ana.(choice_pos{chpos}) = intersect(pos_selected_per_trial.(choice_pos{chpos}),pos_correct_per_trials2ana.(pos2ana{pos})); % get choice for all trials where the correct position is contra / ipsi / center
                                                    end

                                                    if strcmp(tr_con2ana{tr_con},'targ_distr_1HF') 
                                                                                 % KK: todo: important for varibles up or down!!!
                                                    fname_num = ['num_choice_',choice_pos{chpos}];
                                                    fname_prop = ['prop_choice_',choice_pos{chpos}];
                                                    num_choice = numel(pos_selected_per_trials2ana.(choice_pos{chpos})); %Number of trials
                                                    prop_choice= numel(pos_selected_per_trials2ana.(choice_pos{chpos}))/num_total_trials;
                                                    data_struct_per_run.(tr_con2ana{tr_con}).(pos2ana{pos}).(stim_con{stim}){col,file}.(fname_num) = num_choice;
                                                    data_struct_per_run.(tr_con2ana{tr_con}).(pos2ana{pos}).(stim_con{stim}){col,file}.(fname_prop) = prop_choice;
                                                    else
                                                    fname_num = ['num_choice_',choice_pos{chpos}];
                                                    fname_prop = ['prop_choice_',choice_pos{chpos}];
                                                    num_choice = numel(pos_selected_per_trials2ana.(choice_pos{chpos})); %Number of trials
                                                    prop_choice= numel(pos_selected_per_trials2ana.(choice_pos{chpos}))/num_total_trials;
                                                    data_struct_per_run.(tr_con2ana{tr_con}).(pos2ana{pos}).(stim_con{stim}){col,file}.(fname_num) = num_choice;
                                                    data_struct_per_run.(tr_con2ana{tr_con}).(pos2ana{pos}).(stim_con{stim}){col,file}.(fname_prop) = prop_choice;
                                                    end
                                                end
                                            end
                                            
                                            %%% LATENCY %%%
                                            
                                            if Settings.analyze_latency
                                                
                                                % latency if any peripheral target is selected
                                                peripheral_targ_sel = trials2ana(~isnan(targ_selected) & targ_selected~=3); % calculate mean latency of trials2ana with a target selected that was not the central fix spot
                                                
                                                % latency of hits
                                                if strcmp(tr_con2ana{tr_con},'distr_distr_2HF') || strcmp(tr_con2ana{tr_con},'distr_distr_1HF') || strcmp(tr_con2ana{tr_con},'single_distr')
                                                    data_struct_per_run.(tr_con2ana{tr_con}).(pos2ana{pos}).(stim_con{stim}){col,file}.latencies_hits_per_run = NaN; % no latency for maintaining fixation can be measured -> set it to NaN
                                                else
                                                    lat_hits = latency(hits);
                                                    data_struct_per_run.(tr_con2ana{tr_con}).(pos2ana{pos}).(stim_con{stim}){col,file}.latencies_hits_per_run = lat_hits; % all trials in this run concatenated
                                                end
                                                mean_latency_hits_per_run           = nanmean(lat_hits);
                                                SEM_latency_hits_per_run            = nanstd(lat_hits) / sqrt(size(lat_hits,2));
                                                data_struct_per_run.(tr_con2ana{tr_con}).(pos2ana{pos}).(stim_con{stim}){col,file}.mean_latency_hits_per_run = mean_latency_hits_per_run; % mean per run
                                                data_struct_per_run.(tr_con2ana{tr_con}).(pos2ana{pos}).(stim_con{stim}){col,file}.SEM_latency_hits_per_run = SEM_latency_hits_per_run;
                                                
                                                % latency if an incorrect target is selected
                                                incor_periph_targ_sel       = setdiff(peripheral_targ_sel,hits); % peripheral target selected but no hit, i.e. wrong peripheral target selected
                                                incor_contra_choice         = intersect(pos_selected_per_trials2ana.contra,incor_periph_targ_sel); % trials with incorrect left choice
                                                incor_ipsi_choice           = intersect(pos_selected_per_trials2ana.ipsi,incor_periph_targ_sel); % trials with incorrect ipsi choice
                                                lat_incor_contra_choice_per_run                     = single(latency(incor_contra_choice));
                                                lat_incor_ipsi_choice_per_run                       = single(latency(incor_ipsi_choice));
                                                data_struct_per_run.(tr_con2ana{tr_con}).(pos2ana{pos}).(stim_con{stim}){col,file}.latencies_incorrect_contra_choice_per_run = lat_incor_contra_choice_per_run;
                                                data_struct_per_run.(tr_con2ana{tr_con}).(pos2ana{pos}).(stim_con{stim}){col,file}.latencies_incorrect_ipsi_choice_per_run = lat_incor_ipsi_choice_per_run;
                                                
                                                mean_latency_incorrect_contra_choice_per_run        = single(nanmean(lat_incor_contra_choice_per_run));
                                                mean_latency_incorrect_ipsi_choice_per_run          = single(nanmean(lat_incor_ipsi_choice_per_run));
                                                data_struct_per_run.(tr_con2ana{tr_con}).(pos2ana{pos}).(stim_con{stim}){col,file}.mean_latency_incorrect_contra_choice_per_run = mean_latency_incorrect_contra_choice_per_run;
                                                data_struct_per_run.(tr_con2ana{tr_con}).(pos2ana{pos}).(stim_con{stim}){col,file}.mean_latency_incorrect_ipsi_choice_per_run = mean_latency_incorrect_ipsi_choice_per_run;
                                                
                                                % latency of contra and ipsi stimuli independent whether response is correct or incorrect
                                                lat_contra_choice_per_run                           =  single(latency(pos_selected_per_trials2ana.contra));
                                                lat_ipsi_choice_per_run                             =  single(latency(pos_selected_per_trials2ana.ipsi));
                                                data_struct_per_run.(tr_con2ana{tr_con}).(pos2ana{pos}).(stim_con{stim}){col,file}.latencies_contra_choice_per_run = lat_contra_choice_per_run; % for the bar plots (targ-targ especially)
                                                data_struct_per_run.(tr_con2ana{tr_con}).(pos2ana{pos}).(stim_con{stim}){col,file}.latencies_ipsi_choice_per_run = lat_ipsi_choice_per_run; % for the bar plots (targ-targ especially)
                                                
                                                mean_latency_contra_choice_per_run                  = single(nanmean(lat_contra_choice_per_run));
                                                mean_latency_ipsi_choice_per_run                    = single(nanmean(lat_ipsi_choice_per_run));
                                                data_struct_per_run.(tr_con2ana{tr_con}).(pos2ana{pos}).(stim_con{stim}){col,file}.mean_latency_contra_choice_per_run = mean_latency_contra_choice_per_run;
                                                data_struct_per_run.(tr_con2ana{tr_con}).(pos2ana{pos}).(stim_con{stim}){col,file}.mean_latency_ipsi_choice_per_run = mean_latency_ipsi_choice_per_run;
                                                
                                            end
                                            %%% CORRECT TRIAL EYE TRACES %%%
                                            
                                            if Settings.analyze_eye_traces
                                                state = 4;
                                                traces2plot.(stim_con{stim}) = struct('x',NaN,'y',NaN);
                                                traces2plot = get_traces(Settings,traces2plot,eye_traces,trial,trials2ana, state, stim_con{stim});
                                                data_struct_per_run.(tr_con2ana{tr_con}).(pos2ana{pos}).(stim_con{stim}){col,file}.(['traces_',pos2ana{pos},'_choice']) = traces2plot;
                                                
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                    clear trials2ana_tmp1
                    clear trials2ana_tmp2
                    if warning_count ~= 0
                        warning(['For ',num2str(warning_count),' condition combinations the specified conditions the number of hits and/or number of trials with a selected target is zero!'])
                    end
                end
            end % end if -condition to check task-type
        end
        if isempty(data_struct_per_run) % the variable data_struct_per_run is only not empty if all criteria are met in at least one run
            warning(['In session ',Settings.folders{fol},' no trials met the specified criteria.'])
            continue;
        end
        
        
        
%         %% KK - change the number of runs related to inactivation - Pre and Postinjection
%         if strcmp(Settings.Experiment , 'Inactivation')
%             for tr_con = 1:num_tr_con2ana
%                 position = fieldnames(data_struct_per_run.(tr_con2ana{tr_con}));
%                 for pos = 1: numel(position);
%                     for stim = 1:num_stim_con
%                         %stim = 1 -> stim_off                     
%                        if stim == 1    %stim = 1 -> stim_off
%                             for run_PreInjection = 1:  numel(data_struct_per_run.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}))
%                                 if ~ismember(run_PreInjection, Settings.Datasets.PreInjection_Runs{fol})
%                                     data_struct_per_run.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){1,run_PreInjection} = [];
%                                 end
%                             end
%                         elseif  stim == 2
%                             for run_PreInjection = 1:  numel(Settings.Datasets.PreInjection_Runs{fol})
%                                 data_struct_per_run.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){1,run_PreInjection} = [];
%                             end
%                         end
%                         
%                     end
%                 end
%             end            
%         end
        %%
        %%% PER SESSION %%%
        
        for tr_con = 1:num_tr_con2ana % different task conditions
            position = fieldnames(data_struct_per_run.(tr_con2ana{tr_con}));
            for pos = 1: numel(position); % for all positions in this trial condition
                for stim = 1:num_stim_con
                    for col = 1:size(data_struct_per_run.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}),1); % for all colors (rows in data_struct) in this trial condition
                        tmp_data = [data_struct_per_run.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,:}]; % get data out of the structure and make it accessible
                        if isempty(tmp_data) % if tmp_data is empty, i.e. in this session and respective conditions no data was analyzed (all runs were discarded by Setting criteria)
                            data_struct_per_session.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,fol} = [];
                        else
                            %KK
                            if strcmp(Settings.Experiment , 'Inactivation')
                                indEmpty = NaN;  ind = 1; 

                                if strcmp(stim_con{stim}, 'stim_off')%PreInjection
                                    %How many runs of task?
                                    if isnan(Settings.PreInjection_Runs{fol}) % All Runs will be PreInjection because there it is training without clear devision
                                        num_runs =  size(tmp_data,2);
                                    else
                                        total_numRuns_preInjection = numel(Settings.PreInjection_Runs{fol});
                                        for Runs = 1: total_numRuns_preInjection
                                            indEmpty(Runs) = ~isempty(data_struct_per_run.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){1,Runs});
                                        end
                                        num_runs =sum(indEmpty);
                                    end
                                elseif strcmp(stim_con{stim}, 'go_signal')%PostInjection
                                    total_numRuns_preInjection = numel(Settings.PreInjection_Runs{fol});
                                    for Runs = (total_numRuns_preInjection +1): size(data_struct_per_run.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}),2)
                                        indEmpty(ind) = ~isempty(data_struct_per_run.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){1,Runs});
                                        ind = ind +1;
                                    end
                                    num_runs =sum(indEmpty);
                                else
                                    num_runs = size(tmp_data,2);
                                end
                            elseif strcmp(Settings.Experiment , 'Microstimulation')
                                num_runs = size(tmp_data,2);
                            end
                            
                            data_struct_per_session.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,fol}.data_per_run = tmp_data; % rows are runs
                            data_struct_per_session.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,fol}.session      = Settings.folders(fol); % which session is this
                            data_struct_per_session.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,fol}.G_R_ratios   = tmp_data.G_R_ratios; % rows are runs
                            data_struct_per_session.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,fol}.num_runs     = num_runs;
                            num_hits = nansum([tmp_data.num_hits]);
                            num_total_trials = nansum([tmp_data.num_total_trials]);
                            num_initiated_trials = nansum([tmp_data.num_initiated_trials]);
                            data_struct_per_session.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,fol}.num_hits             = num_hits;
                            data_struct_per_session.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,fol}.num_total_trials     = num_total_trials;
                            data_struct_per_session.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,fol}.num_initiated_trials = num_initiated_trials;
                            %%% HITRATE %%%
                            if Settings.analyze_hitrate
                                
                                mean_hitrate_across_runs = nanmean([tmp_data.hitrate]);
                                data_struct_per_session.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,fol}.mean_hitrate_across_runs = mean_hitrate_across_runs; % mean of hitrates per run
                                data_struct_per_session.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,fol}.hitrate_per_session = num_hits / num_total_trials; % this is used to calculate mean hitrate across sessions
                                data_struct_per_session.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,fol}.SEM_hitrate = ...
                                    nanstd([tmp_data.hitrate]) / sqrt(num_runs); % std of hitrate / sqrt of number of runs per session
                            end
                            
                            %%% CHOICE %%%
                            if Settings.analyze_choice
                                for chpos = 1:numel(choice_pos)
                                    
                                    fname_num = ['num_choice_',choice_pos{chpos}];
                                    fname_prop = ['prop_choice_',choice_pos{chpos}];
                                    fname_mean_prop_across_runs = ['mean_prop_choice_across_runs_',choice_pos{chpos}];
                                    fname_SEM_choice = ['SEM_choice_',choice_pos{chpos}];
                                    fname_prop_per_session = ['true_prop_choice_per_session_',choice_pos{chpos}];
                                    
                                    num_choice = nansum([tmp_data.(fname_num)]);
                                    mean_prop_choice_across_runs = nanmean([tmp_data.(fname_prop)]);
                                    prop_choice_per_session = num_choice / num_total_trials; % this is used to calculate mean choices across sessions and for the plot of true proportions per session
                                    
                                    data_struct_per_session.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,fol}.(fname_num) = num_choice;
                                    data_struct_per_session.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,fol}.(fname_mean_prop_across_runs) = mean_prop_choice_across_runs;
                                    data_struct_per_session.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,fol}.(fname_SEM_choice)= ...
                                        nanstd([tmp_data.(fname_prop)]) / sqrt(num_runs); % std of choice proportion / sqrt of number of runs per session
                                    data_struct_per_session.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,fol}.(fname_prop_per_session) = prop_choice_per_session;
                                end
                            end
                            
                            %%% ERROR RATE %%%
                            if Settings.analyze_error_rates
                                % errors are not calculated for positions 'left' and 'right' in distr-distr and targ-targ trials & not for 'combined' in the other trial conditions
                                if isfield(errors.by_condition.per_run.(tr_con2ana{tr_con}),position{pos})
                                    tmp_errors = errors.by_condition.per_run.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,fol};
                                    % take sum of errors in all runs per session and divide by total number of initiated trials
                                    data_struct_per_session.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,fol}.error_rate_fixation_break_per_session = ... % Fixation break rate per session
                                        asin(sqrt(sum(tmp_errors.fbre) / ...
                                        data_struct_per_session.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,fol}.num_initiated_trials));
                                    
                                    data_struct_per_session.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,fol}.error_rate_target_acquis_per_session = ... % Target hold rate per session
                                        asin(sqrt(sum(tmp_errors.tacq) / ...
                                        data_struct_per_session.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,fol}.num_initiated_trials));
                                    
                                    data_struct_per_session.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,fol}.error_rate_target_hold_per_session = ... % Target hold rate per session
                                        asin(sqrt(sum(tmp_errors.thol) / ...
                                        data_struct_per_session.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,fol}.num_initiated_trials));
                                end
                            end
                            
                            %%% LATENCY %%%
                            if Settings.analyze_latency
                                
                                %%% latency of correct saccades (no fixation trials) %%%
                                latencies_hits_per_session = [tmp_data.latencies_hits_per_run];
                                mean_latency_hits_per_session = nanmean(latencies_hits_per_session); % all trials of this session concatenated, i.e. mean per session
                                SEM_latency_hits_per_session = nanstd(latencies_hits_per_session) / sqrt(size(latencies_hits_per_session,2)); % all trials of this session concatenated, i.e. SEM per session
                                data_struct_per_session.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,fol}.latencies_hits_per_session       = latencies_hits_per_session; % all trials of this session conatenated
                                data_struct_per_session.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,fol}.mean_latency_hits_per_session    = mean_latency_hits_per_session;
                                data_struct_per_session.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,fol}.SEM_latency_hits_per_session     = SEM_latency_hits_per_session;
                                
                                mean_latency_hits_across_runs = nanmean([tmp_data.mean_latency_hits_per_run]); % mean of mean latencies per run, i.e. mean across runs
                                SEM_latency_hits_across_runs = nanstd([tmp_data.mean_latency_hits_per_run]) / sqrt(size([tmp_data.mean_latency_hits_per_run],2));
                                data_struct_per_session.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,fol}.mean_latency_hits_across_runs    = mean_latency_hits_across_runs;
                                data_struct_per_session.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,fol}.SEM_latency_hits_across_runs     = SEM_latency_hits_across_runs;
                                
                                %%% latency of incorrect saccades to left or right (no fixation trials) %%%
                                num_incorrect_contra_choices_per_session = size([tmp_data.latencies_incorrect_contra_choice_per_run],2);
                                num_incorrect_ipsi_choices_per_session = size([tmp_data.latencies_incorrect_ipsi_choice_per_run],2);
                                latencies_incor_contra_per_session = [tmp_data.latencies_incorrect_contra_choice_per_run];
                                latencies_incor_ipsi_per_session = [tmp_data.latencies_incorrect_ipsi_choice_per_run];
                                mean_latency_incor_contra_per_session = nanmean(latencies_incor_contra_per_session);
                                mean_latency_incor_ipsi_per_session = nanmean(latencies_incor_ipsi_per_session);
                                SEM_latency_incor_contra_per_session = nanstd(latencies_incor_contra_per_session) / sqrt(size(latencies_incor_contra_per_session,2));
                                SEM_latency_incor_ipsi_per_session = nanstd(latencies_incor_ipsi_per_session) / sqrt(size(latencies_incor_ipsi_per_session,2));
                                
                                % all trials of this session concatenated
                                data_struct_per_session.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,fol}.num_incorrect_contra_choices_per_session = num_incorrect_contra_choices_per_session;
                                data_struct_per_session.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,fol}.latencies_incorrect_contra_choice_per_session = latencies_incor_contra_per_session;
                                data_struct_per_session.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,fol}.mean_latency_incorrect_contra_choice_per_session = mean_latency_incor_contra_per_session;
                                data_struct_per_session.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,fol}.SEM_latency_incorrect_contra_choice_per_session = SEM_latency_incor_contra_per_session;
                                
                                % mean of means per run
                                mean_latency_incor_contra_across_runs = nanmean([tmp_data.mean_latency_incorrect_contra_choice_per_run]);
                                SEM_latency_incor_contra_across_runs = nanstd([tmp_data.mean_latency_incorrect_contra_choice_per_run]) / sqrt(size([tmp_data.mean_latency_incorrect_contra_choice_per_run],2));
                                data_struct_per_session.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,fol}.mean_latency_incorrect_contra_choice_across_runs = mean_latency_incor_contra_across_runs;
                                data_struct_per_session.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,fol}.SEM_latency_incorrect_contra_choice_across_runs = SEM_latency_incor_contra_across_runs;
                                
                                % all trials of this session concatenated
                                data_struct_per_session.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,fol}.num_incorrect_ipsi_choices_per_session = num_incorrect_ipsi_choices_per_session;
                                data_struct_per_session.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,fol}.latencies_incorrect_ipsi_choice_per_session = latencies_incor_ipsi_per_session;
                                data_struct_per_session.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,fol}.mean_latency_incorrect_ipsi_choice_per_session = mean_latency_incor_ipsi_per_session;
                                data_struct_per_session.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,fol}.SEM_latency_incorrect_ipsi_choice_per_session = SEM_latency_incor_ipsi_per_session;
                                
                                % mean of means per run
                                mean_latency_incor_ipsi_across_runs = nanmean([tmp_data.mean_latency_incorrect_ipsi_choice_per_run]);
                                SEM_latency_incor_ipsi_across_runs = nanstd([tmp_data.mean_latency_incorrect_ipsi_choice_per_run]) / sqrt(size([tmp_data.mean_latency_incorrect_ipsi_choice_per_run],2));
                                data_struct_per_session.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,fol}.mean_latency_incorrect_ipsi_choice_across_runs = mean_latency_incor_ipsi_across_runs;
                                data_struct_per_session.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,fol}.SEM_latency_incorrect_ipsi_choice_across_runs = SEM_latency_incor_ipsi_across_runs;
                                
                                %%% latency of all saccades to contra or ipsi (no fixation trials) for bar plots (targ-targ especially) %%%
                                num_contra_choices_per_session = size([tmp_data.latencies_contra_choice_per_run],2);
                                num_ipsi_choices_per_session = size([tmp_data.latencies_ipsi_choice_per_run],2);
                                latencies_contra_per_session = [tmp_data.latencies_contra_choice_per_run];
                                latencies_ipsi_per_session = [tmp_data.latencies_ipsi_choice_per_run];
                                mean_latency_contra_per_session = nanmean(latencies_contra_per_session);
                                mean_latency_ipsi_per_session = nanmean(latencies_ipsi_per_session);
                                SEM_latency_contra_per_session = nanstd(latencies_contra_per_session) / sqrt(size(latencies_contra_per_session,2));
                                SEM_latency_ipsi_per_session = nanstd(latencies_ipsi_per_session) / sqrt(size(latencies_ipsi_per_session,2));
                                
                                % all trials of this session concatenated
                                data_struct_per_session.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,fol}.num_contra_choices_per_session = num_contra_choices_per_session;
                                data_struct_per_session.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,fol}.latencies_contra_choice_per_session = latencies_contra_per_session;
                                data_struct_per_session.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,fol}.mean_latency_contra_choice_per_session = mean_latency_contra_per_session;
                                data_struct_per_session.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,fol}.SEM_latency_contra_choice_per_session = SEM_latency_contra_per_session;
                                
                                % mean of means per run
                                mean_latency_contra_across_runs = nanmean([tmp_data.mean_latency_contra_choice_per_run]);
                                SEM_latency_contra_across_runs = nanstd([tmp_data.mean_latency_contra_choice_per_run]) / sqrt(size([tmp_data.mean_latency_contra_choice_per_run],2));
                                data_struct_per_session.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,fol}.mean_latency_contra_choice_across_runs = mean_latency_contra_across_runs;
                                data_struct_per_session.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,fol}.SEM_latency_contra_choice_across_runs = SEM_latency_contra_across_runs;
                                
                                % all trials of this session concatenated
                                data_struct_per_session.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,fol}.num_ipsi_choices_per_session = num_ipsi_choices_per_session;
                                data_struct_per_session.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,fol}.latencies_ipsi_choice_per_session = latencies_ipsi_per_session;
                                data_struct_per_session.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,fol}.mean_latency_ipsi_choice_per_session = mean_latency_ipsi_per_session;
                                data_struct_per_session.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,fol}.SEM_latency_ipsi_choice_per_session = SEM_latency_ipsi_per_session;
                                
                                % mean of means per run
                                mean_latency_ipsi_across_runs = nanmean([tmp_data.mean_latency_ipsi_choice_per_run]);
                                SEM_latency_ipsi_across_runs = nanstd([tmp_data.mean_latency_ipsi_choice_per_run]) / sqrt(size([tmp_data.mean_latency_ipsi_choice_per_run],2));
                                data_struct_per_session.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,fol}.mean_latency_ipsi_choice_across_runs = mean_latency_ipsi_across_runs;
                                data_struct_per_session.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,fol}.SEM_latency_ipsi_choice_across_runs = SEM_latency_ipsi_across_runs;
                                
                            end
                            
                            %%% CORRECT TRIAL EYE TRACES %%%
                            
                            if Settings.analyze_eye_traces
                                state = 4;
                                traces2plot.(stim_con{stim}) = struct('x',NaN,'y',NaN);
                                traces2plot = get_traces(Settings,traces2plot,eye_traces,trial,trials2ana, state, stim_con{stim});
                                data_struct_per_session.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,fol}.(['traces_',pos2ana{pos},'_choice']) = [];
                                data_struct_per_run.(tr_con2ana{tr_con}).(pos2ana{pos}).(stim_con{stim}){col,file}.(['traces_',pos2ana{pos},'_choice']) = traces2plot;
                                
                            end
                        end
                    end
                end
            end
        end
        
    end
    
    %%% ACROSS SESSIONS %%%
    
    if Settings.analyze_across_sessions
        for tr_con = 1:num_tr_con2ana
            position = fieldnames(data_struct_per_session.(tr_con2ana{tr_con}));
            for pos = 1: numel(position); % for all positions in this trial condition
                for stim = 1:num_stim_con
                    sessions = {};
                    for col = 1:size(data_struct_per_session.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}),1); % for all colors (rows in data_struct) in this trial condition
                        tmp_data = [data_struct_per_session.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,:}]; % get data out of the structure and make it accessible
                        if isempty(tmp_data) % if tmp_data is empty, i.e. in this session and respective conditions no data was analyzed (all runs were discarded by Setting criteria)
                            data_struct_per_session.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,1} = [];
                        else
                            num_sessions = size(tmp_data,2);
                            sessions{end+1} = [tmp_data.session];
                            data_struct_across_sessions.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,1}.data_per_session   = tmp_data; % rows are runs
                            data_struct_across_sessions.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,1}.session            = sessions{1,1}; % which sessions went into this data
                            data_struct_across_sessions.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,1}.G_R_ratios         = tmp_data.G_R_ratios; % rows are runs
                            data_struct_across_sessions.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,1}.num_sessions       = num_sessions;
                            num_hits                                                    = nansum([tmp_data.num_hits]);
                            num_total_trials                                            = nansum([tmp_data.num_total_trials]);
                            data_struct_across_sessions.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,1}.num_hits           = num_hits;
                            data_struct_across_sessions.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,1}.num_total_trials   = num_total_trials;
                            
                            %%% HITRATE %%%
                            if Settings.analyze_hitrate
                                mean_hitrate_across_sessions = nanmean([tmp_data.hitrate_per_session]);
                                data_struct_across_sessions.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,1}.mean_hitrate_across_sessions = mean_hitrate_across_sessions; % mean of hitrates per run
                                data_struct_across_sessions.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,1}.hitrate_per_all_sessions = num_hits / num_total_trials; % mean of hitrates per run
                                data_struct_across_sessions.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,1}.SEM_hitrate = ...
                                    nanstd([tmp_data.hitrate_per_session]) / sqrt(num_sessions); % std of hitrate / sqrt of number of runs per session
                            end
                            
                            %%% CHOICE %%%
                            if Settings.analyze_choice
                                for chpos = 1:numel(choice_pos)
                                    
                                    fname_num                       = ['num_choice_',choice_pos{chpos}];
                                    fname_prop_per_session          = ['true_prop_choice_per_session_',choice_pos{chpos}];
                                    fname_prop_across_sessions      = ['true_prop_choice_across_sessions_',choice_pos{chpos}];
                                    fname_mean_prop_across_sessions = ['mean_prop_choice_across_sessions_',choice_pos{chpos}];
                                    fname_SEM_choice                = ['SEM_choice_',choice_pos{chpos}];
                                    
                                    num_choice                          = nansum([tmp_data.(fname_num)]);% 20 6 30
                                    
                                    
                                    
                                    mean_prop_choice_across_sessions    = nanmean([tmp_data.(fname_prop_per_session)]); %0.5000    0.4000    0.5000
                                    prop_choice_across_sessions         = num_choice / num_total_trials; % this is used to calculate mean choices across sessions and for the plot of true proportions per session
                                    
                                    data_struct_across_sessions.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,1}.(fname_num)                        = num_choice;
                                    data_struct_across_sessions.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,1}.(fname_mean_prop_across_sessions)  = mean_prop_choice_across_sessions;
                                    data_struct_across_sessions.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,1}.(fname_SEM_choice)= ...
                                        nanstd([tmp_data.(fname_prop_per_session)]) / sqrt(num_sessions); % std of choice proportion / sqrt of number of runs per session
                                    % std will be zero if there is only one session in this condition
                                    data_struct_across_sessions.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,1}.(fname_prop_across_sessions) = prop_choice_across_sessions;
                                end
                            end
                            
                            %%% LATENCY %%%
                            if Settings.analyze_latency
                                
                                %%% latency of correct saccades (no fixation trials) %%%
                                % all trials of all sessions concatenated
                                latencies_hits_all_sessions = [tmp_data.latencies_hits_per_session];
                                mean_latency_hits_all_sessions = nanmean(latencies_hits_all_sessions);
                                SEM_latency_hits_all_sessions = nanstd(latencies_hits_all_sessions) / sqrt(size(latencies_hits_all_sessions,2));
                                data_struct_across_sessions.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,1}.latencies_hits_all_sessions = latencies_hits_all_sessions; % all latencies of all sessions of this condition concatenated
                                data_struct_across_sessions.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,1}.mean_latency_hits_all_sessions = mean_latency_hits_all_sessions; % all latencies of all sessions of this condition concatenated and one mean created across all sessions
                                data_struct_across_sessions.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,1}.SEM_latency_hits_all_sessions = SEM_latency_hits_all_sessions;
                                
                                % mean of means per session, i.e. mean across session
                                mean_latency_hits_across_sessions = nanmean([tmp_data.mean_latency_hits_per_session]); % mean of mean latencies per session, i.e. mean across sessions
                                SEM_latency_hits_across_sessions = nanstd([tmp_data.mean_latency_hits_per_session]) / sqrt(size([tmp_data.mean_latency_hits_per_session],2));
                                data_struct_across_sessions.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,1}.mean_latency_hits_across_sessions = mean_latency_hits_across_sessions;
                                data_struct_across_sessions.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,1}.SEM_latency_hits_across_sessions = SEM_latency_hits_across_sessions;
                                
                                %%% latency of incorrect saccades to contra or ipsi (no fixation trials) %%%
                                num_incorrect_contra_choices_all_sessions = size([tmp_data.latencies_incorrect_contra_choice_per_session],2);
                                num_incorrect_ipsi_choices_all_sessions = size([tmp_data.latencies_incorrect_ipsi_choice_per_session],2);
                                latencies_incor_contra_all_sessions = [tmp_data.latencies_incorrect_contra_choice_per_session];
                                latencies_incor_ipsi_all_sessions = [tmp_data.latencies_incorrect_ipsi_choice_per_session];
                                mean_latency_incor_contra_all_sessions = nanmean(latencies_incor_contra_all_sessions);
                                mean_latency_incor_ipsi_all_sessions = nanmean(latencies_incor_ipsi_all_sessions);
                                SEM_latency_incor_contra_all_sessions = nanstd(latencies_incor_contra_all_sessions) / sqrt(size(latencies_incor_contra_all_sessions,2));
                                SEM_latency_incor_ipsi_all_sessions = nanstd(latencies_incor_ipsi_all_sessions) / sqrt(size(latencies_incor_ipsi_all_sessions,2));
                                
                                % all trials of all sessions concatenated
                                data_struct_across_sessions.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,1}.num_incorrect_contra_choices_all_sessions = num_incorrect_contra_choices_all_sessions;
                                data_struct_across_sessions.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,1}.latencies_incorrect_contra_choice_all_sessions = latencies_incor_contra_all_sessions;
                                data_struct_across_sessions.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,1}.mean_latency_incorrect_contra_choice_all_sessions = mean_latency_incor_contra_all_sessions;
                                data_struct_across_sessions.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,1}.SEM_latency_incorrect_contra_choice_all_sessions = SEM_latency_incor_contra_all_sessions;
                                
                                % mean of means per session, i.e. mean across session
                                mean_latency_incor_contra_across_sessions = nanmean([tmp_data.mean_latency_incorrect_contra_choice_per_session]);
                                SEM_latency_incor_contra_across_sessions = nanstd([tmp_data.mean_latency_incorrect_contra_choice_per_session]) / sqrt(size([tmp_data.mean_latency_incorrect_contra_choice_per_session],2));
                                data_struct_across_sessions.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,1}.mean_latency_incorrect_contra_choice_across_sessions = mean_latency_incor_contra_across_sessions;
                                data_struct_across_sessions.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,1}.SEM_latency_incorrect_contra_choice_across_sessions = SEM_latency_incor_contra_across_sessions;
                                
                                % all trials of all sessions concatenated
                                data_struct_across_sessions.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,1}.num_incorrect_ipsi_choices_all_sessions = num_incorrect_ipsi_choices_all_sessions;
                                data_struct_across_sessions.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,1}.latencies_incorrect_ipsi_choice_all_sessions = latencies_incor_ipsi_all_sessions;
                                data_struct_across_sessions.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,1}.mean_latency_incorrect_ipsi_choice_all_sessions = mean_latency_incor_ipsi_all_sessions;
                                data_struct_across_sessions.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,1}.SEM_latency_incorrect_ipsi_choice_all_sessions = SEM_latency_incor_ipsi_all_sessions;
                                
                                % mean of means per session, i.e. mean across session
                                mean_latency_incor_ipsi_across_sessions = nanmean([tmp_data.mean_latency_incorrect_ipsi_choice_per_session]); % mean of mean latencies per session, i.e. mean across sessions
                                SEM_latency_incor_ipsi_across_sessions = nanstd([tmp_data.mean_latency_incorrect_ipsi_choice_per_session]) / sqrt(size([tmp_data.mean_latency_incorrect_ipsi_choice_per_session],2));
                                data_struct_across_sessions.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,1}.mean_latency_incorrect_ipsi_choice_across_sessions = mean_latency_incor_ipsi_across_sessions;
                                data_struct_across_sessions.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,1}.SEM_latency_incorrect_ipsi_choice_across_sessions = SEM_latency_incor_ipsi_across_sessions;
                                
                                % latency of all saccades to contra or ipsi (no fixation trials) for bar plots (targ-targ especially)
                                num_contra_choices_all_sessions = size([tmp_data.latencies_contra_choice_per_session],2);
                                num_ipsi_choices_all_sessions = size([tmp_data.latencies_ipsi_choice_per_session],2);
                                latencies_contra_all_sessions = [tmp_data.latencies_contra_choice_per_session];
                                latencies_ipsi_all_sessions = [tmp_data.latencies_ipsi_choice_per_session];
                                mean_latency_contra_all_sessions = nanmean(latencies_contra_all_sessions);
                                mean_latency_ipsi_all_sessions = nanmean(latencies_ipsi_all_sessions);
                                SEM_latency_contra_all_sessions = nanstd(latencies_contra_all_sessions) / sqrt(size(latencies_contra_all_sessions,2));
                                SEM_latency_ipsi_all_sessions = nanstd(latencies_ipsi_all_sessions) / sqrt(size(latencies_ipsi_all_sessions,2));
                                
                                % all trials of all sessions concatenated
                                data_struct_across_sessions.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,1}.num_contra_choices_all_sessions = num_contra_choices_all_sessions;
                                data_struct_across_sessions.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,1}.latencies_contra_choice_all_sessions = latencies_contra_all_sessions;
                                data_struct_across_sessions.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,1}.mean_latency_contra_choice_all_sessions = mean_latency_contra_all_sessions;
                                data_struct_across_sessions.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,1}.SEM_latency_contra_choice_all_sessions = SEM_latency_contra_all_sessions;
                                
                                % mean of means per session, i.e. mean across session
                                mean_latency_contra_across_sessions = nanmean([tmp_data.mean_latency_contra_choice_per_session]);
                                SEM_latency_contra_across_sessions = nanstd([tmp_data.mean_latency_contra_choice_per_session]) / sqrt(size([tmp_data.mean_latency_contra_choice_per_session],2));
                                data_struct_across_sessions.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,1}.mean_latency_contra_choice_across_sessions = mean_latency_contra_across_sessions;
                                data_struct_across_sessions.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,1}.SEM_latency_contra_choice_across_sessions = SEM_latency_contra_across_sessions;
                                
                                % all trials of all sessions concatenated
                                data_struct_across_sessions.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,1}.num_ipsi_choices_all_sessions = num_ipsi_choices_all_sessions;
                                data_struct_across_sessions.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,1}.latencies_ipsi_choice_all_sessions = latencies_ipsi_all_sessions;
                                data_struct_across_sessions.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,1}.mean_latency_ipsi_choice_all_sessions = mean_latency_ipsi_all_sessions;
                                data_struct_across_sessions.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,1}.SEM_latency_ipsi_choice_all_sessions = SEM_latency_ipsi_all_sessions;
                                
                                % mean of means per session, i.e. mean across session
                                mean_latency_ipsi_across_sessions = nanmean([tmp_data.mean_latency_ipsi_choice_per_session]);
                                SEM_latency_ipsi_across_sessions = nanstd([tmp_data.mean_latency_ipsi_choice_per_session]) / sqrt(size([tmp_data.mean_latency_ipsi_choice_per_session],2));
                                data_struct_across_sessions.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,1}.mean_latency_ipsi_choice_across_sessions = mean_latency_ipsi_across_sessions;
                                data_struct_across_sessions.(tr_con2ana{tr_con}).(position{pos}).(stim_con{stim}){col,1}.SEM_latency_ipsi_choice_across_sessions = SEM_latency_ipsi_across_sessions;
                            end
                            
                        end
                    end
                end
            end
        end
    end
end

%%% ERROR RATE CALCULATION ACROSS SESSIONS %%%
if Settings.analyze_error_rates && ~isempty(errors.by_condition) % if errors were calculated (which is not the case if a dataset without errors is loaded)
    tr_cond = fieldnames(errors.by_condition.per_run);
    for tr_con = 1:numel(tr_cond)
        position = fieldnames(errors.by_condition.per_run.(tr_cond{tr_con}));
        for pos = 1:numel(position)
            for st_con = 1:numel(stim_con)
                colors = size(errors.by_condition.per_run.(tr_cond{tr_con}).(position{pos}).(stim_con{st_con}),1);
                for col = 1:colors
                    if strcmp(tr_cond{tr_con},'targ_targ_2HF') || strcmp(tr_cond{tr_con},'single_targ')
                        col = 1; % there is no distractor color so write in first row
                    end
                    
                    %                     tmp_data = errors.by_condition.per_run.(tr_cond{tr_con}).(position{pos}).(stim_con{st_con}){col,fol};
                    %
                    %                     data_struct_per_session.(tr_cond{tr_con}).(position{pos}).(stim_con{st_con}){col,fol}.error_rate_fixation_break_per_session = ... % Fixation break rate per session
                    %                         sum(tmp_data.fbre) / ...
                    %                         data_struct_per_session.(tr_cond{tr_con}).(position{pos}).(stim_con{st_con}){col,fol}.num_total_trials;
                    %
                    %                     data_struct_per_session.(tr_cond{tr_con}).(position{pos}).(stim_con{st_con}){col,fol}.error_rate_target_acquis_per_session = ... % Target hold rate per session
                    %                         sum(tmp_data.tacq) / ...
                    %                         data_struct_per_session.(tr_cond{tr_con}).(position{pos}).(stim_con{st_con}){col,fol}.num_total_trials;
                    %
                    %                     data_struct_per_session.(tr_cond{tr_con}).(position{pos}).(stim_con{st_con}){col,fol}.error_rate_target_hold_per_session = ... % Target hold rate per session
                    %                         sum(tmp_data.thol) / ...
                    %                         data_struct_per_session.(tr_cond{tr_con}).(position{pos}).(stim_con{st_con}){col,fol}.num_total_trials;
                    
                    
                    
                    %                     data_struct_per_session.(tr_cond{tr_con}).(position{pos}).(stim_con{st_con}){col,fol}.mean_error_rate_fixation_break_across_sessions = ... % Fixation break mean error rate across sessions
                    %                         mean(data_struct_per_session.(tr_cond{tr_con}).(position{pos}).(stim_con{st_con}){col,fol}.error_rate_fixation_break);
                    %                     data_struct_per_session.(tr_cond{tr_con}).(position{pos}).(stim_con{st_con}){col,fol}.mean_error_rate_target_hold_across_sessions = ... % Target hold mean error rate across sessions
                    %                         mean(data_struct_per_session.(tr_cond{tr_con}).(position{pos}).(stim_con{st_con}){col,fol}.error_rate_target_hold);
                    %
                    %                     data_struct_per_session.(tr_cond{tr_con}).(position{pos}).(stim_con{st_con}){col,fol}.SEM_error_rate_target_hold = [];...
                    %                         %%%% SEM needs all error_rate_fixation_breaks of
                    %                         %%%% all sessions
                    %                 nanstd(latencies_hits_all_sessions) / sqrt(size(latencies_hits_all_sessions,2));
                    %
                    %                     errors.by_condition.across_sessions.(tr_cond{tr_con}).(position{pos}).(stim_con{st_con}){col,fol}.fbre = sum(sum(tmp_data.fbre));
                    %                     errors.by_condition.across_sessions.(tr_cond{tr_con}).(position{pos}).(stim_con{st_con}){col,fol}.thol = sum(sum(tmp_data.thol));
                    
                    %%% Call the plot function here or a statistics function
                end
            end
        end
    end
end

%%% SAVE DATA STRUCTURES TO DATASET %%%

if ~isempty(Settings.save_dataset) % save data structure to a dataset that can be loaded without analyzing all the run files again
    save_dir = get_save_dir(Settings.folders, Settings);
    if ~exist(save_dir, 'dir')
        mkdir(save_dir); % create directory if it does not exist yet
    end
    t = datetime('now','TimeZone','local','Format','yyyyMMdd-HHmm');
    t = datestr(t,'yyyymmdd-HHMM'); % convert to string with specific output format
    file_name = fullfile(save_dir,[Settings.save_dataset,'_',t]);
    save(file_name,'data_struct_per_run','data_struct_per_session','data_struct_across_sessions','stim_con') % save as .mat file
    disp(['Dataset saved: ',file_name]);
end
%%%%%%%%%%%%%%%%%%%%%
%%% STATISTICS %%%
%%%%%%%%%%%%%%%%%%
if Settings.statistics_on
    
    [stat_results_anova,stat_struct_ttest, stat_struct_ttest_anova] = distr_task_behavior_statistics(data_struct_per_session, Settings, stim_con);
    
    % chi-square tests for separate sessions %LG
    if ~Settings.analyze_across_sessions
        stat_struct_chi2 = distr_task_behavior_chisquare(data_struct_per_session, stim_con, Settings);
    else stat_struct_chi2 = [];
    end
    
    %%% SAVE STATISTICS OUTPUT %%%
    if Settings.save_stat_output
        save_dir = get_save_dir(Settings.folders, Settings);
        if ~exist(char(save_dir), 'dir')
            mkdir(save_dir); % create directory if it does not exist yet
        end
        %         t = datetime('now','TimeZone','local','Format','yyyyMMdd-HHmm');
        %         t = datestr(t,'yyyymmdd-HHMM'); % convert to string with specific output format
        t = date;
        if ~isempty(stat_struct_ttest)
            file_name = char(fullfile(save_dir,['stat_struct_ttest_',t]));
            save(file_name,'stat_struct_ttest') % save as .mat file
            disp(['T-test results saved: ',file_name]);
            fname1 = fieldnames(stat_struct_ttest);
            for ano = 1:numel(fname1) % anova
                if ~strcmp(fname1{ano},'legend') % exclude the 'legend' fieldname which does not contain anova results
                    fname2 = fieldnames(stat_struct_ttest.(fname1{ano}));
                    for int = 1:numel(fname2) % interactions
                        fname3 = fieldnames(stat_struct_ttest.(fname1{ano}).(fname2{int}));
                        for dv = 1:numel(fname3)
                            excel_file_name = [file_name,'.xlsx'];
                            file = stat_struct_ttest.(fname1{ano}).(fname2{int}).(fname3{dv});
                            file = [cell(2,size(file,2));file]; % add two empty line to 'file'
                            file{1,1} = fname2{int}; % writes tested interaction in 1st cell in 1st row
                            file{2,1} = fname3{dv}; % writes dependent variable into 1st cell in 2nd row
                            sheet = fname1{ano};
                            % save_in_excel(excel_file_name,sheet,file); % save as excel file
                        end
                    end
                end
            end
        end
        if ~isempty(stat_results_anova)
            file_name = char(fullfile(save_dir,['stat_results_anova_',t]));
            save(file_name,'stat_results_anova') % save as .mat file
            disp(['ANOVA results saved: ',file_name]);
            fname1 = fieldnames(stat_results_anova);
            %             for ano = 1:numel(fname1)
            %                 fname2 = fieldnames(stat_results_anova.(fname1{ano}));
            %                 for dv = 1:numel(fname2)
            %                     excel_file_name = [file_name,'.xlsx'];
            %                     file = stat_results_anova.(fname1{ano}).(fname2{dv});
            %                     file = df_to_txt(file);
            %                     file = [cell(2,size(file,2));file]; % add two empty line to 'file'
            %                     file{1,1} = ''''; % creates an empty line in excel that is not empty to matlab
            %                     file{2,1} = fname2{dv}; % writes dependent variable into 1st cell in 2nd row
            %                     sheet = fname1{ano};
            %                     save_in_excel(excel_file_name,sheet,file); % save as excel file
            %                     disp(['ANOVA results saved: ',excel_file_name]);
            %                 end
            %             end
            
        end
    end
end
%%%%%%%%%%%%%%%%%%
%%% PLOTTING %%%
%%%%%%%%%%%%%%%
set(0,'defaulttextinterpreter','none') % sets interpreter of all texts, legends and titles to none
% Change default axes fonts.
set(0,'DefaultAxesFontName', 'Arial')
set(0,'DefaultAxesFontSize', 12)

% Change default text fonts.
set(0,'DefaultTextFontname', 'Arial')
set(0,'DefaultTextFontSize', 12)

plot_per_or_across = [Settings.plot_per_session,Settings.plot_across_sessions];
plot_hitrate_choice_or_latency = [Settings.analyze_hitrate,...
    Settings.analyze_choice,Settings.analyze_choice,Settings.analyze_choice,...
    Settings.analyze_latency,Settings.analyze_latency,Settings.analyze_latency];
ydata_fname = {'mean_hitrate_across_runs','mean_hitrate_across_sessions';...
    'mean_prop_choice_across_runs_contra','mean_prop_choice_across_sessions_contra';...
    'mean_prop_choice_across_runs_ipsi','mean_prop_choice_across_sessions_ipsi';...
    'mean_prop_choice_across_runs_fixation','mean_prop_choice_across_sessions_fixation';...
    'mean_latency_hits','mean_latency_hits';...
    'mean_latency_incorrect_contra_choice','mean_latency_incorrect_contra_choice';...
    'mean_latency_incorrect_ipsi_choice','mean_latency_incorrect_ipsi_choice'};
SEM_fname = {'SEM_hitrate';...
    'SEM_choice_contra';...
    'SEM_choice_ipsi';...
    'SEM_choice_fixation';...
    'SEM_latency_hits';...
    'SEM_latency_incorrect_contra_choice';...
    'SEM_latency_incorrect_ipsi_choice'};
ydata_fname_PSF = {'hitrate_per_session','hitrate_per_all_sessions';...
    '','';...
    '','';...
    '','';...
    '','';...
    '','';...
    '',''};
fit_fname = {'num_hits';...
    '';...
    '';...
    '';...
    '';...
    '';...
    ''};
ylabl = {'Hitrate';...
    'contra choices';...
    'ipsi choices';...
    'Fixation choices';...
    'Latency of hits';...
    'Latency of incorrect contra choices';...
    'Latency of incorrect ipsi choices'};
ntrials = {'num_total_trials';...
    'num_total_trials';...
    'num_total_trials';...
    'num_total_trials';...
    'num_hits';...
    'num_incorrect_contra_choices';...
    'num_incorrect_ipsi_choices'};
n = {'num_runs','num_sessions'};

%%%%%%%%%%%%%%%%%%%%%%%%
%%% LINE PLOT %%%%%%%%%
if Settings.line_plot
    for p_a = 1:numel(plot_per_or_across)
        if plot_per_or_across(p_a) % if the current plot_per_session or plot_across_sessions is true
            for h_c_l = 1:numel(plot_hitrate_choice_or_latency)
                if plot_hitrate_choice_or_latency(h_c_l)
                    for fol = 1:numel(Settings.folders)
                        switch p_a
                            case 1 % if plot per session
                                data_struct = data_struct_per_session;
                                session_s = Settings.folders(fol);
                            case 2 % if plot across sessions
                                if fol > 1 % in data_struct_across_sessions there is only one column -> fol needs to be 1
                                    continue;
                                end
                                data_struct = data_struct_across_sessions;
                                session_s = Settings.folders;
                        end
                        plot_linegraph(Settings,data_struct,ydata_fname{h_c_l,p_a},SEM_fname{h_c_l,1},ydata_fname_PSF{h_c_l,p_a},fit_fname{h_c_l,1},ylabl{h_c_l,1},ntrials{h_c_l,1},n{1,p_a},fol,session_s)
                    end
                end
            end
        end
    end
end

%%%%%%%%%%%%%%%
%%% BAR PLOT %%% - the loop which generates all hitrate / RT  plots for all
%%% conditions
if Settings.bar_plot || Settings.bar_plot_II || Settings.hist_RTs % plot_struct_double and -_single is needed for the histogram of RTs
    plot_struct_double = []; plot_struct_single = [];
    plot_struct_double_II = []; plot_struct_single_II = [];
    for p_a = 1:numel(plot_per_or_across)
        if plot_per_or_across(p_a) % if the current plot_per_session or plot_across_sessions is true
            switch p_a
                case 1 % if plot per session
                    data_struct = data_struct_per_session;
                case 2 % if plot across sessions
                    data_struct = data_struct_across_sessions;
            end
            if Settings.bar_plot || Settings.hist_RTs
                [plot_struct_double, plot_struct_single] = getbarplotstruct (Settings,plot_struct_double,plot_struct_single,data_struct,p_a,stim_con); % get plot_struct
            end
            if Settings.bar_plot_II
                [plot_struct_double_II, plot_struct_single_II] = getbarplotstruct_II (Settings,plot_struct_double_II,plot_struct_single_II,data_struct,p_a,stim_con); % get plot_struct
            end
        end
    end
    if Settings.bar_plot
        plot_bargraph(Settings,plot_struct_double,stat_struct_ttest) % plot double stimuli trials
        plot_bargraph(Settings,plot_struct_single,stat_struct_ttest) % plot single stimulus trials
    end
    if Settings.bar_plot_II
        plot_bargraph_II(Settings,plot_struct_double_II,stat_struct_ttest) % plot double stimuli trials
        plot_bargraph_II(Settings,plot_struct_single_II,stat_struct_ttest) % plot single stimulus trials
    end
else
    plot_struct_double = []; plot_struct_single = [];
end

%%% HISTOGRAM OF REACTION TIMES %%%
if Settings.analyze_latency && Settings.hist_RTs
    hist_RTs(Settings,plot_struct_double)
    hist_RTs(Settings,plot_struct_single)
end

%%% PLOT TACQ ERROR TRACES OF IPSI_TARG - CONTRA EASY DISTRACTOR TRIALS %%%
if Settings.analyze_error_rates && Settings.plot_error_traces
    stim_con = fieldnames(error_data2plot);
    h_tar_pos_up = figure;
    h_tar_pos_zero = figure;
    h_tar_pos_down = figure;
    fig_handles = [h_tar_pos_up;h_tar_pos_zero;h_tar_pos_down];
    fig_handle_fnames = {'h_tar_pos_up','h_tar_pos_zero','h_tar_pos_down'};
    stim_colors = {[0.7 0.7 0.7];[0 1 1];[0 0.8 0];[0.8 0.5 0]};
    radius = 5;
    target_color = 'r'; % Settings.target_color_dim./255;
    distr_color = [0.9 0.7 0]; % [60,60,0]./255;
    for stim = 1:numel(stim_con)
        ntrials = size(error_data2plot.(stim_con{stim}),2);
        for tr = 2:ntrials % first column (first 'trial') is NaN
            if ~isequal(fieldnames(error_data2plot.(stim_con{stim})(tr)),{'x';'y'}) % other fieldnames have been added if there were error trials in this condition
                hold on;
                x = error_data2plot.(stim_con{stim})(tr).x;
                y = error_data2plot.(stim_con{stim})(tr).y;
                tar_pos_x = error_data2plot.(stim_con{stim})(tr).tar_pos_x; tar_pos_y = error_data2plot.(stim_con{stim})(tr).tar_pos_y;
                distr_pos_x = error_data2plot.(stim_con{stim})(tr).distr_pos_x; distr_pos_y = error_data2plot.(stim_con{stim})(tr).distr_pos_y;
                if tar_pos_y == 0 % ipsi target at
                    plt = 1;
                elseif tar_pos_y < 0 % ipsi target below horizontal midline
                    plt = 2;
                elseif tar_pos_y > 0 % right target above horizontal midline
                    plt = 3;
                end
                figure(fig_handles(plt))
                sh.(fig_handle_fnames{plt}) = subplot(1,2,1);
                circle(tar_pos_x,tar_pos_y,radius,'Color',[237, 28, 36]./255,'LineWidth',2); % plot target
                circle(distr_pos_x,distr_pos_y,radius,'Color',[255 212 38]./255,'LineWidth',2); % plot distractor
                ha.(stim_con{stim}) = plot(x,y,'Color',stim_colors{stim});
                subplot (2,2,2) % plot x position
                plot(x,'Color',stim_colors{stim})
                hold on;
                subplot (2,2,4) % plot x position
                plot(y,'Color',stim_colors{stim})
                hold on;
            end
        end
    end
    for plt = 1:3
        figure(fig_handles(plt))
        % FIGURE SIZE %
        fig_p = fig_handles(plt).Position;
        fig_p(3) = fig_p(3)*1.5;
        fig_handles(plt).Position = fig_p;
        
        subplot(1,2,1)
        fix_color = [0.7 0.7 0.7]; % Settings.fixation_color_dim./255;
        circle(0,0,radius,'Color',fix_color,'LineWidth',2); % plot distractor
        xlim([-25,25]);
        ylim([-25,25]);
        set(gcf,'color','w'); % sets background color to white
        title('Tacq errors - Ipsi_targ-contra_easy_distr');
        % LEGEND %
        ha_fnames = fieldnames(ha);
        num_ha_fnames = numel(ha_fnames);
        pl_ha = zeros(num_ha_fnames,1);
        [ha_fnames_new, idx] = sortifmember(ha_fnames,stim_con);
        for fn = 1:num_ha_fnames
            pl_ha(fn) = ha.(ha_fnames_new{fn}); % get all necessary handles for the legend
        end
        lh = legend(pl_ha,stim_con(idx),'Location','southeast');
        set(lh,'interpreter','none') % turns off latex interpreter for title, set to 'tex' to turn it back on
        % SUBPLOT SIZE %
        p = get(sh.(fig_handle_fnames{plt}),'position');
        p(3) = p(4)*0.5; p(1) = 0.1;
        set(sh.(fig_handle_fnames{plt}),'position',p);
        hold off;
        xlabel('Horizontal eye position');
        ylabel('Vertical eye position');
        % plot([1,5],[1,5],'Color',[0.7 0.7 0.7],'LineWidth',2)
        
        subplot(2,2,2)
        xlim([0,550]);
        ylim([-25,25]);
        xlabel('Time [ms]');
        ylabel('Horizontal eye position');
        
        subplot(2,2,4)
        xlim([0,550]);
        ylim([-25,25]);
        xlabel('Time [ms]');
        ylabel('Vertical eye position');
    end
    for plt = 1:3
        figure(fig_handles(plt))
        %%% SAVE PLOT %%%
        if Settings.save_figures
            file_name_part = {'Target_up_right','Target_mid_right','Target_down_right'};
            file_name = ['Errors - Tacq - Ipsi-targ - contra-easy-distr',file_name_part{plt}];
            save_dir = get_save_dir(Settings.folders, Settings);
            save_figure_as(fig_handles(plt), Settings, file_name, {'.ai','.png'}, save_dir) % save
            if Settings.close_after_saving
                close(gcf)
            end
        end
    end
end

disp('Done :)')
end


%% Subfunctions

function [out_sorted,idx] = sortifmember(in2sort,sortby)
sort_idx = 1;
for elem = 1:numel(sortby)
    if ismember(sortby{elem},in2sort)
        out_sorted{sort_idx} = sortby{elem};
        idx(sort_idx) = elem;
        sort_idx = sort_idx + 1;
    end
end
if sort_idx == 1 % if no sorting has been done
    out_sorted = in2sort;
    idx = [];
    warning('No sorting could be done. No element of the variable to sort were found in the variable to sort by.')
end
end

function [errors] = count_errors_by_position(Settings,errors,trial,relevant_stimulus,trial_idx,curr_tr_cond,curr_session,curr_run)
for trial_no = 1:numel(trial_idx)
    txt = [];
    for tar = 1:numel(relevant_stimulus{trial_idx(trial_no)})
        xy_pos = ['x','y'];
        for xy = 1:numel(xy_pos)
            posi = trial(trial_idx(trial_no)).task.eye.tar(relevant_stimulus{trial_idx(trial_no)}(tar)).(xy_pos(xy));
            if isempty(txt)
                txt = [txt,num2str(posi)];
            else
                txt = [txt,';',num2str(posi)];
            end
            all_positions{trial_no,1} = txt; % get x and y values of all relevant stimuli in every trial
        end
    end
end
positions = unique(all_positions);
abort_codes = {'STIMULUS POSITIONS','HITS','ABORT_EYE_FIX_ACQ_STATE','ABORT_EYE_FIX_HOLD_STATE','ABORT_EYE_TAR_ACQ_STATE','ABORT_EYE_TAR_HOLD_STATE','ABORT_WRONG_TARGET_SELECTED'};
errors_mtx= num2cell(zeros(numel(positions),numel(abort_codes)-1));
errors_mtx = [positions,errors_mtx]; % put position values in the first column
errors_mtx = [abort_codes;errors_mtx]; % put abort_codes as header in the first row
Fn_curr_ses = ['Ses_',curr_session];
Fn_curr_run = ['Run_',curr_run];
% per run
for trial_no = 1:numel(trial_idx)
    [~,pos_idx] = ismember(all_positions{trial_no,1},positions); % get index of where positions of relevant stimuli are in the current trial
    pos_idx = pos_idx + 1; % header is first row, so 1 needs to be added to the row index
    if trial(trial_idx(trial_no)).success
        errors_mtx{pos_idx,2} = errors_mtx{pos_idx,2} +1;
    else
        [~,abort_idx] = ismember(trial(trial_idx(trial_no)).abort_code,abort_codes);
        errors_mtx{pos_idx,abort_idx} = errors_mtx{pos_idx,abort_idx} +1;
    end
end
% per session
if isempty(errors)
    errors.(curr_tr_cond).(Fn_curr_ses) = errors_mtx; % data per session
else
    tr_cond = fieldnames(errors);
    [~,con_idx] = ismember(curr_tr_cond,tr_cond);
    if con_idx ~= 0 % if the current trial condition is already in the errors structure
        session = fieldnames(errors.(tr_cond{con_idx}));
        [~,ses_idx] = ismember(Fn_curr_ses,session);
        if ses_idx ~= 0 % if the current session is already in the errors structure
            %%%%%%% complete conditions
            for pos = 2: numel(errors.(curr_tr_cond).(Fn_curr_ses)(:,1)) % for every column that contains values of errors - position of stimuli
                for abo = 2:numel(errors.(curr_tr_cond).(Fn_curr_ses)(1,:)) % for every row that contains values of errors - abort-code of trial
                    errors.(curr_tr_cond).(Fn_curr_ses){pos,abo} = errors.(curr_tr_cond).(Fn_curr_ses){pos,abo} + errors_mtx{pos,abo}; % add errors of this run to errors of previous run of the same session
                end
            end
        else
            errors.(curr_tr_cond).(Fn_curr_ses) = errors_mtx;
        end
    else
        errors.(curr_tr_cond).(Fn_curr_ses) = errors_mtx;
    end
    
end

errors.(curr_tr_cond).(Fn_curr_run) = errors_mtx; % data per run; put abort_codes as header in the first row

for pos = 1:size(positions,1)
    
end
end

function [errors,error_data2plot] = count_errors_by_condition(Settings,errors,trial,pos_correct_per_trial,trial_cond,stim_cond,trial_fix_break_error,trial_tar_hold_error,trial_tar_acq_error,distractor_color,eye_traces,error_data2plot,fol,file)
st_cond = fieldnames(stim_cond);
tr_cond = fieldnames(trial_cond);
for tr_con = 1:numel(tr_cond)
    if strcmp(tr_cond{tr_con},'targ_distr_1HF') || strcmp(tr_cond{tr_con},'targ_distr_2HF') || strcmp(tr_cond{tr_con},'single_targ') || strcmp(tr_cond{tr_con},'single_distr')
        position = {'contra','ipsi'};
    elseif strcmp(tr_cond{tr_con},'targ_targ_1HF') || strcmp(tr_cond{tr_con},'targ_targ_2HF') || strcmp(tr_cond{tr_con},'distr_distr_1HF') || strcmp(tr_cond{tr_con},'distr_distr_2HF')
        position = {'combined'};
    end
    for pos = 1:numel(position)
        for st_con = 1:numel(st_cond)
            for col = 1:size(distractor_color,1)
                trials = intersect(trial_cond.(tr_cond{tr_con}),stim_cond.(st_cond{st_con})); % get trials of one trial and stimulation condition
                trials = intersect(trials,pos_correct_per_trial.(position{pos})); % get trials separately for target / distr position, respectively
                if strcmp(tr_cond{tr_con},'distr_distr_2HF') ||strcmp(tr_cond{tr_con},'distr_distr_1HF') || strcmp(tr_cond{tr_con},'single_distr') || strcmp(tr_cond{tr_con},'targ_distr') % if this is a trial with a distractor
                    trials = intersect(trials,distractor_color(col,:)); % get only trials with the respective distractor color
                end
                errors.per_run.(tr_cond{tr_con}).(position{pos}).(st_cond{st_con}){col,fol}.fbre(file) = numel(intersect(trials,trial_fix_break_error));
                tacq_errors = intersect(trials,trial_tar_acq_error);
                errors.per_run.(tr_cond{tr_con}).(position{pos}).(st_cond{st_con}){col,fol}.tacq(file) = numel(tacq_errors);
                errors.per_run.(tr_cond{tr_con}).(position{pos}).(st_cond{st_con}){col,fol}.thol(file) = numel(intersect(trials,trial_tar_hold_error));
                
                %%% error traces %%%
                if strcmp(tr_cond{tr_con},'targ_distr') && strcmp(position{pos},'ipsi') && col==2 % because only in this condition there is a significant effect of stimulation on tacq error rate
                    
                    error_data2plot = get_traces(Settings,error_data2plot,eye_traces,trial,tacq_errors, 4, st_cond{st_con});
                    
                    % this will do the same thing as the get_traces
                    % function, but manually; WORKS!!
                    %                     if ~isempty(tacq_errors)
                    %                         for tr = 1:numel(tacq_errors)
                    %                             count = size(error_data2plot.(st_cond{st_con}),2) + 1;
                    %                             error_data2plot.(st_cond{st_con})(1,count).x = error_traces(tacq_errors(tr)).x_eye(error_traces(tacq_errors(tr)).states==4); % x-values during the target acquisition state
                    %                             error_data2plot.(st_cond{st_con})(1,count).y = error_traces(tacq_errors(tr)).y_eye(error_traces(tacq_errors(tr)).states==4); % y-values during the target acquisition state
                    %                             clear red_target;
                    %                             for i = 1:size(trial(tacq_errors(tr)).task.eye.tar,2)
                    %                                 if isequal(trial(tacq_errors(tr)).task.eye.tar(i).color_dim,Settings.target_color_dim) % if target
                    %                                     red_target.x = trial(tacq_errors(tr)).task.eye.tar(i).x;
                    %                                     red_target.y = trial(tacq_errors(tr)).task.eye.tar(i).y;
                    %                                 end
                    %                                 if ~isequal(trial(tacq_errors(tr)).task.eye.tar(i).color_dim,Settings.target_color_dim) && ~isequal(trial(tacq_errors(tr)).task.eye.tar(i).color_dim,Settings.fixation_color_dim) % not target nor fixation spot
                    %                                     distr.x = trial(tacq_errors(tr)).task.eye.tar(i).x;
                    %                                     distr.y = trial(tacq_errors(tr)).task.eye.tar(i).y;
                    %                                 end
                    %                             end
                    %                             error_data2plot.(st_cond{st_con})(1,count).tar_pos_x = red_target.x;
                    %                             error_data2plot.(st_cond{st_con})(1,count).tar_pos_y = red_target.y;
                    %                             error_data2plot.(st_cond{st_con})(1,count).distr_pos_x = distr.x;
                    %                             error_data2plot.(st_cond{st_con})(1,count).distr_pos_y = distr.y;
                    %                             disp(['count:',num2str(count)])
                    %                         end
                    %                     end
                end
            end
        end
        
    end
end

end

function traces2plot = get_traces(Settings,traces2plot,eye_traces,trial,trials, state, st_con)
if ~isempty(trials)
    for tr = 1:numel(trials)
        count = size(traces2plot.(st_con),2) + 1;
        traces2plot.(st_con)(1,count).x = eye_traces(trials(tr)).x_eye(eye_traces(trials(tr)).states==state); % x-values during the target acquisition state
        traces2plot.(st_con)(1,count).y = eye_traces(trials(tr)).y_eye(eye_traces(trials(tr)).states==state); % y-values during the target acquisition state
        clear red_target; clear distr;
        for i = 1:size(trial(trials(tr)).task.eye.tar,2) % find the position of the red target ...
            if isequal(trial(trials(tr)).task.eye.tar(i).color_dim,Settings.target_color_dim) % if target
                red_target.x = trial(trials(tr)).task.eye.tar(i).x;
                red_target.y = trial(trials(tr)).task.eye.tar(i).y;
            end
            if ~isequal(trial(trials(tr)).task.eye.tar(i).color_dim,Settings.target_color_dim) && ~isequal(trial(trials(tr)).task.eye.tar(i).color_dim,Settings.fixation_color_dim) % not target nor fixation spot
                distr.x = trial(trials(tr)).task.eye.tar(i).x; % ... and the distractor
                distr.y = trial(trials(tr)).task.eye.tar(i).y;
            end
        end
        if ~exist('red_target','var') % if no red target exists in this trial, say position is NaN
            red_target.x = NaN;
            red_target.y = NaN;
        end
        if ~exist('distr','var') % if no distractor exists in this trial, say position is NaN
            distr.x = NaN;
            distr.y = NaN;
        end
        traces2plot.(st_con)(1,count).tar_pos_x = red_target.x;
        traces2plot.(st_con)(1,count).tar_pos_y = red_target.y;
        traces2plot.(st_con)(1,count).distr_pos_x = distr.x;
        traces2plot.(st_con)(1,count).distr_pos_y = distr.y;
        disp(['count:',num2str(count)])
    end
end
end

function plot_linegraph(Settings,data_struct,ydata_fname,SEM_fname,ydata_fname_PSF,fit_fname,ylabl,ntrials,n,fol,session_s)
stim = {'stim_off'};
num_stim_con = numel(stim);
plot_colors = {'ob','og','or','ok'};
% for fol = 1:numel(session_s)
for tr_con = 1:numel(Settings.analyze_trial_conditions) % for every condition to analyze as in Settings.analyze_trial_conditions
    fig = figure; hold on; clear ha
    keep_cur_fig = false;
    con_in_this_plot = {}; sessions_in_this_plot = {}; xlimits = [];
    for con = 1:numel(Settings.analyze_trial_conditions{tr_con}) % second level condition; plot targ-targ or single-targ trials first and then the trials that contained distractors
        position = fieldnames(data_struct.(Settings.analyze_trial_conditions{tr_con}{con}));
        num_pos = numel(position);
        for pos = 1:num_pos; % for all positions in this trial condition ot be analyzed
            switch position{pos}
                case 'contra'
                    data_col = plot_colors{1};
                case 'ipsi'
                    data_col = plot_colors{3};
                case'combined'
                    data_col = plot_colors{4};
            end
            for stim_con = 1:num_stim_con
                num_col = size([data_struct.(Settings.analyze_trial_conditions{tr_con}{con}).(position{pos}).(stim{stim_con}){:,fol}],2);
                x_data = nan(num_col,1); y_data = nan(num_col,1); SEM = nan(num_col,1); txt = {}; txt_y_pos = nan(num_col,1); NumPos = nan(1,num_col); OutOfNum = nan(1,num_col);
                for col = 1:num_col % for all colors (rows in data_struct) in this trial condition
                    
                    x_data(col) = data_struct.((Settings.analyze_trial_conditions{tr_con}{con})).(position{pos}).(stim{stim_con}){col,fol}.G_R_ratios(col);
                    
                    if Settings.fit_psychometric_curve && ~isempty(fit_fname) % Plot PSF for hitrate only: fit_fname is empty for choice and latency plots
                        % this will plot the same data as is fitted
                        % SEM(col) = NaN;
                        % txt_y_pos(col) = y_data(col) * 0.96; % y position to place text in plot
                        % y_data(col) = data_struct.((Settings.analyze_trial_conditions{tr_con}{con})).(position{pos}).(stim{stim_con}){col,fol}.(ydata_fname_PSF);
                        
                        % this will plot means of hitrate_per_session, i.e. mean across sessions,
                        SEM(col) = data_struct.((Settings.analyze_trial_conditions{tr_con}{con})).(position{pos}).(stim{stim_con}){col,fol}.(SEM_fname);
                        NumPos(col) = data_struct.((Settings.analyze_trial_conditions{tr_con}{con})).(position{pos}).(stim{stim_con}){col,fol}.(fit_fname);
                        OutOfNum(col) = data_struct.((Settings.analyze_trial_conditions{tr_con}{con})).(position{pos}).(stim{stim_con}){col,fol}.num_total_trials;
                        y_data(col) = data_struct.((Settings.analyze_trial_conditions{tr_con}{con})).(position{pos}).(stim{stim_con}){col,fol}.(ydata_fname);
                        txt_y_pos(col) = (y_data(col) - SEM(col)) * 0.96; % y position to place text in plot % same as in 'else' right below
                    else % don't plot PSF -> mean values per run/session and SEMs are plotted
                        y_data(col) = data_struct.((Settings.analyze_trial_conditions{tr_con}{con})).(position{pos}).(stim{stim_con}){col,fol}.(ydata_fname);
                        SEM(col) = data_struct.((Settings.analyze_trial_conditions{tr_con}{con})).(position{pos}).(stim{stim_con}){col,fol}.(SEM_fname);
                        txt_y_pos(col) = (y_data(col) - SEM(col)) * 0.96; % y position to place text in plot
                    end
                    if size(y_data,1) ~= 0 || all(isnan(y_data)) % if there is data to plot then add these sessions / this session to the list of sessions for title and filename
                        sessions_in_this_plot = [sessions_in_this_plot;data_struct.((Settings.analyze_trial_conditions{tr_con}{con})).(position{pos}).(stim{stim_con}){col,fol}.session];
                    end
                    
                    txt_tmp = sprintf('n=%g \nntrials=%g',...
                        data_struct.((Settings.analyze_trial_conditions{tr_con}{con})).(position{pos}).(stim{stim_con}){col,fol}.(n),...
                        data_struct.((Settings.analyze_trial_conditions{tr_con}{con})).(position{pos}).(stim{stim_con}){col,fol}.(ntrials));
                    txt = [txt;txt_tmp];
                    
                end
                if numel(y_data) == 0 || all(isnan(y_data)) % if there is no data to plot or if all data is NaN
                    continue
                else
                    %%% PLOT THE DATA %%%
                    keep_cur_fig = true;
                    ha.(position{pos}) = errorbar(x_data,y_data,SEM,data_col);
                    vers = version; % get current matlab version
                    if strcmp(vers(end-6:end-1),'R2014b')
                        % MATLAB 2014b to 2016a
                        % Create errorbar
                        X = x_data;
                        Y = y_data;
                        E = SEM;
                        % Width of the top and bottom lines of errorbar
                        xlength = 0.02;
                        % Make horizontal lines with 'line'
                        for k = 1:length(X)
                            x = [X(k) - xlength, X(k) + xlength];
                            y_h = [Y(k) + E(k), Y(k) + E(k)];
                            line(x, y_h,'Color',data_col(end),'LineWidth',0.5);
                            y_b = [Y(k) - E(k), Y(k) - E(k)];
                            line(x, y_b,'Color',data_col(end),'LineWidth',0.5);
                        end
                        
                    elseif strcmp(vers(end-6:end-1),'R2012a')
                        % up to MATLAB 2014a
                        hb = get(ha.(position{pos}),'children');
                        Xdata = get(hb(2),'Xdata');
                        temp = 4:3:length(Xdata);
                        temp(3:3:end) = [];
                        % xcontra and xipsi contain the indices of the contra and ipsi
                        %  endpoints of the horizontal lines
                        xcontra = temp; xipsi = temp+1;
                        % Increase line length by 0.2 units
                        Xdata(xcontra) = Xdata(xcontra) - .02;
                        Xdata(xipsi) = Xdata(xipsi) + .02;
                        set(hb(2),'Xdata',Xdata)
                        
                    end
                    % text(x_data,double(txt_y_pos),txt,'Color',data_col(length(data_col):end),'HorizontalAlignment','center','VerticalAlignment','top') % plaace text below data point
                    % text(x_data + 0.03,y_data,txt,'Color',data_col(3:end),'HorizontalAlignment','contra','VerticalAlignment','middle') % place text right of data point
                    
                    %%% CURVE FITTING %%%
                    if Settings.fit_psychometric_curve && ...
                            ~isequal(x_data,0) && ... % if this is not a target-target or a single-target condition which should not be fitted
                            any(OutOfNum) && ... % if there are any  trials for this combination of conditions
                            size(OutOfNum,2) > 2 && ... % Curius had only if there are 5 colors for this condition
                            ~isempty(fit_fname) % fit_fname is empty for latency plots
                        StimLevels = x_data';
                        
                        [pDev, p_converged] = plot_psychometric_curve(StimLevels, NumPos, OutOfNum, session_s, data_col, Settings);
                        disp(Settings.analyze_trial_conditions{tr_con}{con})
                        
                        %%% INFORMATION FOR TITLE AND FILENAME %%%
                        titl_pPSF = ' (PSF)';
                        file_pPSF = '_PSF';
                    else
                        titl_pPSF = '';
                        file_pPSF = '';
                        pDev = NaN;
                        p_converged = NaN;
                    end
                    
                    %%% INFORMATION FOR TITLE AND FILENAME %%%
                    con_in_this_plot{end+1} = Settings.analyze_trial_conditions{tr_con}{con}; % relevant conditions
                    xlimits = [xlimits;[min(x_data),max(x_data)]]; % saves the x values of plot data for the xlim of plots
                end
            end
        end
    end
    if ~keep_cur_fig % if there was no data to plot, otherwise the 'else' statements would create a new figure
        close(gcf) % close the current figure
    else
        %%% SET X AND Y LIM %%%
        % Change default axes fonts.
        set(0,'DefaultAxesFontName', 'Arial')
        set(0,'DefaultAxesFontSize', 12)
        
        % Change default text fonts.
        set(0,'DefaultTextFontname', 'Arial')
        set(0,'DefaultTextFontSize', 12)
        
        xlabel('Distractor G/R ratio','FontSize',12)
        ylabel(ylabl,'FontSize',12)
        xlim(round([min(xlimits(:,1)) - 0.1, max(xlimits(:,2)) + 0.1],1)) % finds min and max of plotted x values, substracts / adds 0.1 and rounds to first decimal
        if strcmp(ylabl,'Hitrate')
            ylim([0 1]) % KK changed
        elseif strcmp(ylabl(1:7),'Latency')
            ylim([0.15 0.3])
        elseif strfind(ylabl,'choice') && isempty(strfind(ylabl,'Latency')) % if 'choice' is in the ylabl but 'Latency' is not
            ylim([0 1])
        end
        for i = 1:numel(x_data)
            x_tick(i,1) = str2double(sprintf('%.2f', x_data(i)));
        end
        set(gca,'XTick',x_tick); % [0:0.2:1.0]
        set(gca,'YTick',[0:0.1:1.0],'FontSize',12);
        set(gcf,'color','w'); % sets background color to white
        
        %%% CREATE PLOT TITLE, FILENAME PARTS AND LEGEND %%%
        con_in_this_plot = sortifmember(unique(con_in_this_plot),Settings.analyze_trial_conditions{tr_con});
        sessions_in_this_plot = unique(sessions_in_this_plot);
        sessions_in_this_plot = strjoin(sessions_in_this_plot',' - ');
        title_pcon = strjoin(con_in_this_plot,' - ');
        title_pses = get_title_session_part(sessions_in_this_plot, Settings);
        % titl_pses = strjoin(sessions_in_this_plot,' - ');
        % titl_pses = [' - session: ',titl_pses];
        title_pdat = [ylabl,titl_pPSF];
        
        if ~isnan(pDev) && ~isnan(p_converged) % if goodness of fit has been successfully calculated
            title_pGOF = sprintf('\n Fit: pDev = %.3g; sim. exp. was fit successfully in %.1f%% simulations',pDev,p_converged);
            file_pGOF = '_GOF';
        else
            title_pGOF = '';
            file_pGOF = '';
        end
        
        titl = [title_pdat,title_pses,' - ',title_pcon,title_pGOF];
        th = title(titl);
        set(th,'interpreter','none') % turns off latex interpreter for title, set to 'tex' to turn it back on
        
        ha_fnames = fieldnames(ha);
        num_ha_fnames = numel(ha_fnames);
        pl_ha = zeros(num_ha_fnames,1);
        legend_order = {'contra' 'ipsi' 'combined'};
        [ha_fnames_new, idx] = sortifmember(ha_fnames,legend_order);
        for fn = 1:num_ha_fnames
            pl_ha(fn) = ha.(ha_fnames_new{fn}); % get all necessary handles for the legend
        end
        if strcmp(Settings.hemisphere_of_stimulation,'left')
            legend_entries = {'ipsi' 'contra' 'all'}; % contra = ipsi, ipsi = contra
        elseif strcmp(Settings.hemisphere_of_stimulation,'right')
            legend_entries = {'contra' 'ipsi' 'all'}; % contra = contra, ipsi = ipsi
        end
        ha_legend = legend(pl_ha,legend_entries(idx),'Location','southeast');
        set(ha_legend,'FontSize',12);
        hold off
        
        %%% SAVE PLOT %%%
        if Settings.save_figures
            % Filename
            if strcmp(Settings.hitrate_mode,'correct target selected')
                file_phit = '_Corr-targ-sel_';
            elseif strcmp(Settings.hitrate_mode,'success')
                file_phit = '_success_';
            end
            file_pses = sessions_in_this_plot; % strjoin(sessions_in_this_plot,'-');
            file_pcon = strjoin(con_in_this_plot,'-');
            file_name = [ylabl,'_',file_pses,file_pPSF,file_pGOF,file_phit,file_pcon];
            % Save
            save_dir = get_save_dir(Settings.folders{fol}, Settings);
            save_figure_as(fig, Settings, file_name, {'.ai','.png'}, save_dir)
            if Settings.close_after_saving
                close(gcf)
            end
        end
    end
end
% end
end

function [pDev, p_converged] = plot_psychometric_curve(StimLevels, NumPos, OutOfNum, session_s, data_col, Settings)

disp('Plotting psychometric function:');
searchGrid.alpha = -1:.01:1; % threshold -1:.01:1;
searchGrid.beta = 10.^[-1:.1:2]; % slope 10.^[-1:.1:2]
searchGrid.gamma = [0:.01:0.5]; % guess-rate [0.33:.01:0.5]
searchGrid.lambda = [0:0.01:0.3]; % lapse-rate [0:0.01:0.3]
lapseLimits = [0 0.3];
guessLimits = [0.33 0.5];
paramsFree = [1 1 1 1]; % logical index for which parameters are free or fixed
PF = @PAL_CumulativeNormal;
options = PAL_minimize('options');
% plot(StimLevels,NumPos,'o');

fit_col = ['--',data_col(end)];
[paramsValues, LL, exitflag, output] = PAL_PFML_Fit(StimLevels, ...
    NumPos, OutOfNum, searchGrid, paramsFree, PF,...
    'lapseLimits',lapseLimits,'guessLimits',guessLimits, ...
    'searchOptions',options); % ,'gammaEQlambda', 1

if exitflag
    fit_x_values = [0:0.01:1];
    fit_y_values = PAL_CumulativeNormal(paramsValues, fit_x_values);
    h = plot(fit_x_values,fit_y_values,fit_col);
    % set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % this turns off the display of this plot (the fitted curve) in the legend
end

if Settings.calc_goodness_of_fit
    % Goodness of Fit
    B = 1000;
    [Dev, pDev, DevSim, converged] = PAL_PFML_GoodnessOfFit(StimLevels, NumPos, OutOfNum, paramsValues, paramsFree, B, PF);%, ...
    %                 'searchGrid', searchGrid, 'lapseLimits',lapseLimits);
    p_converged = sum(converged)/B*100;
    
    %             descr = {'Goodness of Fit:';
    %                 'Proportion of the Deviance values from simulations that were';
    %                 sprintf('greater than Deviance value of data = %.3g',pDev);
    %                 };
    %             text(0.3,0.2,descr)
    %             ylim = get(gca,'YLim');
    %             xlim = get(gca,'XLim');
    %             text(xlim(2),ylim(2)*1/3,'descr', ...
    %                 'VerticalAlignment','bottom', ...
    %                 'HorizontalAlignment','ipsi');
else
    pDev = NaN;
    p_converged = NaN;
end
end

function [plot_struct_double, plot_struct_single] = getbarplotstruct(Settings,plot_struct_double,plot_struct_single,data_struct,p_a,stim_con)
if  strcmp(Settings.Experiment , 'Inactivation') 
    legend_entries = {'Pre','Post'};
elseif  strcmp(Settings.Experiment , 'Microstimulation') 
legend_entries = {'control','-250 ms','-100 ms','+50 ms'};
elseif strcmp(Settings.Experiment , 'Electrophysiology') 
    legend_entries = {''};
end

num_stim_con = numel(stim_con);
single_or_double_stimulus = {'double','single'};
for s_d = 1:numel(single_or_double_stimulus)
    switch single_or_double_stimulus{s_d} % is this going to be a plot of trials with two stimuli or a single stimulus presented
        case 'double' % double stimulus conditions
            hitrate_choice_or_latency = {'choice','latency'};
            plot_hitrate_choice_or_latency = [Settings.analyze_choice,Settings.analyze_latency];
            

            
            if strcmp(Settings.hemisphere_of_stimulation,'left')
             plot_groups = {'targ_targ','distr_distr','ipsi_targ-contra_distr','ipsi_distr-contra_targ'};
            elseif strcmp(Settings.hemisphere_of_stimulation,'right')
            plot_groups = {'targ_targ','distr_distr','contra_targ-ipsi_distr','contra_distr-ipsi_targ'};
                
            end
            
            tr_conditions = {'targ_targ','distr_distr','targ_distr','targ_distr'};
            positions = {{'combined'},{'combined'},{'contra'},{'ipsi'}};
            
        case 'single' % single stimulus conditions
            hitrate_choice_or_latency = {'hitrate','latency'};
            plot_hitrate_choice_or_latency = [Settings.analyze_hitrate,Settings.analyze_latency];
            plot_groups = {'single_targ','single_distr'};
            tr_conditions = {'single_targ','single_distr'};
            
            if strcmp(Settings.hemisphere_of_stimulation,'left') %KK092019
            positions = {{'ipsi','contra'},{'ipsi','contra'}};
            elseif strcmp(Settings.hemisphere_of_stimulation,'right')
            positions = {{'contra','ipsi'},{'contra','ipsi'}};
                
            end
    end
    for h_c_l = 1:numel(plot_hitrate_choice_or_latency)
        if plot_hitrate_choice_or_latency(h_c_l)
            for tr_pos_con = 1:numel(plot_groups) % for each separate plots as defined in plot_groups
                switch hitrate_choice_or_latency{h_c_l}
                    case 'choice'
                        choices = {'contra' 'fixation' 'ipsi'};
                        if strcmp(Settings.bar_plot_choice_mean_or_true_proportion,'mean_and_SEM') % plot the mean and SEM across runs (plot per session) or across sessions (plot across sessions)
                            ydata_fname = {'mean_prop_choice_across_runs_contra','mean_prop_choice_across_sessions_contra';...
                                'mean_prop_choice_across_runs_fixation','mean_prop_choice_across_sessions_fixation';...
                                'mean_prop_choice_across_runs_ipsi','mean_prop_choice_across_sessions_ipsi'};
                        elseif strcmp(Settings.bar_plot_choice_mean_or_true_proportion,'true_proportion') % plot the true proportions per session or across sessions
                            ydata_fname = {'true_prop_choice_per_session_contra','true_prop_choice_across_sessions_contra';...
                                'true_prop_choice_per_session_fixation','true_prop_choice_across_sessions_fixation';...
                                'true_prop_choice_per_session_ipsi','true_prop_choice_across_sessions_ipsi'};
                        end
                        SEM_fname = {'SEM_choice_contra';...
                            'SEM_choice_fixation';...
                            'SEM_choice_ipsi'};
                        ntrials = {'num_choice_contra';...
                            'num_choice_fixation';...
                            'num_choice_ipsi'};
                        n = {'num_runs','num_sessions'};
                    case 'hitrate'
                        choices = Settings.analyze_by_position(1:2);
                        if strcmp(Settings.bar_plot_choice_mean_or_true_proportion,'mean_and_SEM') % plot the mean and SEM
                            ydata_fname = {'mean_hitrate_across_runs','mean_hitrate_across_sessions'};
                        elseif strcmp(Settings.bar_plot_choice_mean_or_true_proportion,'true_proportion') % plot the true proportions per session or across sessions
                            ydata_fname = {'hitrate_per_session','hitrate_per_all_sessions'};
                        end
                        SEM_fname = {'SEM_hitrate'};
                        ntrials = {'num_hits'};
                        n = {'num_runs','num_sessions'};
                    case 'latency'
                        switch single_or_double_stimulus{s_d}
                            case 'double' % double stimulus conditions
                                if strcmp(Settings.bar_plot_latency_avg_mean_or_true_mean,'mean_of_means') % plot the average of mean latency per run (session)
                                    ydata_fname = {'mean_latency_contra_choice_across_runs','mean_latency_contra_choice_across_sessions';...
                                        'mean_latency_ipsi_choice_across_runs','mean_latency_ipsi_choice_across_sessions'};
                                    SEM_fname = {'SEM_latency_contra_choice_across_runs','SEM_latency_contra_choice_across_sessions';...
                                        'SEM_latency_ipsi_choice_across_runs','SEM_latency_ipsi_choice_across_sessions'};
                                elseif strcmp(Settings.bar_plot_latency_avg_mean_or_true_mean,'true_mean_of_concatenated_trials') % plot the mean of all concatenated trials of a session (of all sessions)
                                    ydata_fname = {'mean_latency_contra_choice_per_session','mean_latency_contra_choice_all_sessions';...
                                        'mean_latency_ipsi_choice_per_session','mean_latency_ipsi_choice_all_sessions'};
                                    SEM_fname = {'SEM_latency_contra_choice_per_session','SEM_latency_contra_choice_all_sessions';...
                                        'SEM_latency_ipsi_choice_per_session','SEM_latency_ipsi_choice_all_sessions'};
                                end
                                choices = Settings.analyze_by_position(1:2);
                                ntrials = {'num_contra_choices_per_session','num_contra_choices_all_sessions';...
                                    'num_ipsi_choices_per_session','num_ipsi_choices_all_sessions'};
                                latencies = {'latencies_contra_choice_per_session','latencies_contra_choice_all_sessions';...
                                    'latencies_ipsi_choice_per_session','latencies_ipsi_choice_all_sessions'};
                            case 'single' % single stimulus conditions
                                choices = Settings.analyze_by_position(1:2);
                                switch tr_conditions{tr_pos_con} %
                                    case  'single_targ'
                                        if strcmp(Settings.bar_plot_latency_avg_mean_or_true_mean,'mean_of_means') % plot the average of mean latency per run (session)
                                            ydata_fname = {'mean_latency_hits_across_runs','mean_latency_hits_across_sessions'};
                                            SEM_fname = {'SEM_latency_hits_across_runs','SEM_latency_hits_across_sessions'};
                                        elseif strcmp(Settings.bar_plot_latency_avg_mean_or_true_mean,'true_mean_of_concatenated_trials') % plot the mean of all concatenated trials of a session (of all sessions)
                                            ydata_fname = {'mean_latency_hits_per_session','mean_latency_hits_all_sessions'};
                                            SEM_fname = {'SEM_latency_hits_per_session','SEM_latency_hits_all_sessions'};
                                        end
                                        ntrials = {'num_hits','num_hits'};
                                        latencies = {'latencies_hits_per_session','latencies_hits_all_sessions'};
                                    case  'single_distr' % ch = num of rows = 1; pos defines choice position here
                                        if strcmp(Settings.bar_plot_latency_avg_mean_or_true_mean,'mean_of_means') % plot the average of mean latency per run (session)
                                            ydata_fname = {'mean_latency_incorrect_contra_choice_across_runs','mean_latency_incorrect_contra_choice_across_sessions';...
                                                'mean_latency_incorrect_ipsi_choice_across_runs','mean_latency_incorrect_ipsi_choice_across_sessions'};
                                            SEM_fname = {'SEM_latency_incorrect_contra_choice_across_runs','SEM_latency_incorrect_contra_choice_across_sessions';...
                                                'SEM_latency_incorrect_ipsi_choice_across_runs','SEM_latency_incorrect_ipsi_choice_across_sessions'};
                                        elseif strcmp(Settings.bar_plot_latency_avg_mean_or_true_mean,'true_mean_of_concatenated_trials') % plot the mean of all concatenated trials of a session (of all sessions)
                                            ydata_fname = {'mean_latency_incorrect_contra_choice_per_session','mean_latency_incorrect_contra_choice_all_sessions';...
                                                'mean_latency_incorrect_ipsi_choice_per_session','mean_latency_incorrect_ipsi_choice_all_sessions'};
                                            SEM_fname = {'SEM_latency_incorrect_contra_choice_per_session','SEM_latency_incorrect_contra_choice_all_sessions';...
                                                'SEM_latency_incorrect_ipsi_choice_per_session','SEM_latency_incorrect_ipsi_choice_all_sessions'};
                                        end
                                        ntrials = {'num_incorrect_contra_choices_per_session','num_incorrect_contra_choices_all_sessions';...
                                            'num_incorrect_ipsi_choices_per_session','num_incorrect_ipsi_choices_all_sessions'};
                                        latencies = {'latencies_incorrect_contra_choice_per_session','latencies_incorrect_contra_choice_all_sessions';...
                                            'latencies_incorrect_ipsi_choice_per_session','latencies_incorrect_ipsi_choice_all_sessions'};
                                end
                        end
                        if strcmp(Settings.bar_plot_latency_avg_mean_or_true_mean,'mean_of_means') % plot the average of mean latency per run (session)
                            txt_nruns_ses = {'n=nruns=%g \nntrials=%g','n=nses=%g \nntrials=%g'}; % text to eacht datapoint in the plot
                        elseif strcmp(Settings.bar_plot_latency_avg_mean_or_true_mean,'true_mean_of_concatenated_trials') % plot the mean of all concatenated trials of a session (of all sessions)
                            txt_nruns_ses = {'nruns=%g \nn=ntrials=%g','nses=%g \nn=ntrials=%g'};
                        end
                        nruns_ses = {'num_runs','num_sessions'};
                        
                end
                % find the sessions to plot
                if isempty(Settings.load_data)
                    sessions_idx = 1:numel(Settings.folders); % if no data is loaded all sessions in the settings are plotted
                else
                    tr_condition = fieldnames(data_struct);
                    side = fieldnames(data_struct.(tr_condition{1}));
                    stimulation = fieldnames(data_struct.(tr_condition{1}).(side{1}));
                    array = [data_struct.(tr_condition{1}).(side{1}).(stimulation{1}){1,:}]; %
                    all_loaded_sessions = [array.session];
                    [~,sessions_idx] = ismember(Settings.folders,all_loaded_sessions); % get indices of loaded sessions that are specified in Settings
                end
                for fol = 1:numel(sessions_idx)
                    fol_idx = sessions_idx(fol);
                    if p_a == 2 && fol > 1 % in data_struct_across_sessions there is only one column -> fol needs to be 1
                        continue;
                    end
                    for pos = 1:numel(positions{tr_pos_con})
                        for stim = 1:num_stim_con
                            if ismember(tr_conditions{tr_pos_con},fieldnames(data_struct)) % if the specified trial conditions are actually analyzed and field of data_struct
                                for col = 1:size([data_struct.(tr_conditions{tr_pos_con}).(positions{tr_pos_con}{pos}).(stim_con{stim}){:,fol_idx}],2); % for all colors (rows in data_struct) in this trial condition
                                    for ch = 1:size(ydata_fname,1)
                                        switch hitrate_choice_or_latency{h_c_l}
                                            case 'choice' % put choice as group level
                                                plot_struct_double{col,tr_pos_con}{p_a,fol}.choice.barvalues{ch,stim} = data_struct.(tr_conditions{tr_pos_con}).(positions{tr_pos_con}{pos}).(stim_con{stim}){col,fol_idx}.(ydata_fname{ch,p_a}); % y value to plot
                                                if strcmp(Settings.bar_plot_choice_mean_or_true_proportion,'true_proportion') % plot the true proportions per session or across sessions
                                                    plot_struct_double{col,tr_pos_con}{p_a,fol}.choice.errors{ch,stim} = 0; % SEM value is 0 because no SEM can be calculated from only one proportion
                                                    plot_struct_double{col,tr_pos_con}{p_a,fol}.choice.bar_plot_mean_or_true_proportion = 'true_proportion';
                                                else % plot the mean and SEM
                                                    plot_struct_double{col,tr_pos_con}{p_a,fol}.choice.errors{ch,stim} = data_struct.(tr_conditions{tr_pos_con}).(positions{tr_pos_con}{pos}).(stim_con{stim}){col,fol_idx}.(SEM_fname{ch}); % SEM values
                                                    plot_struct_double{col,tr_pos_con}{p_a,fol}.choice.bar_plot_mean_or_true_proportion = 'mean_and_SEM';
                                                end
                                                plot_struct_double{col,tr_pos_con}{p_a,fol}.choice.groupnames{ch,1}     = choices{ch};
                                                plot_struct_double{col,tr_pos_con}{p_a,fol}.choice.legend{1,stim}       =legend_entries{stim};
                                                plot_struct_double{col,tr_pos_con}{p_a,fol}.choice.session_s            = strjoin(data_struct.(tr_conditions{tr_pos_con}).(positions{tr_pos_con}{pos}).(stim_con{stim}){col,fol_idx}.session,' - '); % session(s) for this data point
                                                plot_struct_double{col,tr_pos_con}{p_a,fol}.choice.tr_pos_con           = plot_groups{tr_pos_con};
                                                plot_struct_double{col,tr_pos_con}{p_a,fol}.choice.tr_con               = tr_conditions{tr_pos_con};
                                                plot_struct_double{col,tr_pos_con}{p_a,fol}.choice.G_R_ratio = data_struct.(tr_conditions{tr_pos_con}).(positions{tr_pos_con}{pos}).(stim_con{stim}){col,fol_idx}.G_R_ratios(col);
                                                plot_struct_double{col,tr_pos_con}{p_a,fol}.choice.txt{ch,stim} = sprintf('n=%g \nntrials=%g',...
                                                    data_struct.(tr_conditions{tr_pos_con}).(positions{tr_pos_con}{pos}).(stim_con{stim}){col,fol_idx}.(n{1,p_a}),...
                                                    data_struct.(tr_conditions{tr_pos_con}).(positions{tr_pos_con}{pos}).(stim_con{stim}){col,fol_idx}.(ntrials{ch,1}));
                                            case 'hitrate' % put pos as group level (position of the presented target)
                                                plot_struct_single{col,tr_pos_con}{p_a,fol}.choice.barvalues{pos,stim} = data_struct.(tr_conditions{tr_pos_con}).(positions{tr_pos_con}{pos}).(stim_con{stim}){col,fol_idx}.(ydata_fname{ch,p_a}); % y value to plot
                                                if strcmp(Settings.bar_plot_choice_mean_or_true_proportion,'true_proportion') % plot the true proportions per session or across sessions
                                                    plot_struct_single{col,tr_pos_con}{p_a,fol}.choice.errors{pos,stim} = 0; % SEM value is 0 because no SEM can be calculated from only one proportion
                                                    plot_struct_single{col,tr_pos_con}{p_a,fol}.choice.bar_plot_mean_or_true_proportion = 'true_proportion';
                                                else % plot the mean and SEM
                                                    plot_struct_single{col,tr_pos_con}{p_a,fol}.choice.errors{pos,stim} = data_struct.(tr_conditions{tr_pos_con}).(positions{tr_pos_con}{pos}).(stim_con{stim}){col,fol_idx}.(SEM_fname{ch}); % SEM values
                                                    plot_struct_single{col,tr_pos_con}{p_a,fol}.choice.bar_plot_mean_or_true_proportion = 'mean_and_SEM';
                                                end
                                                plot_struct_single{col,tr_pos_con}{p_a,fol}.choice.groupnames{pos,1} = choices{pos};
                                                plot_struct_single{col,tr_pos_con}{p_a,fol}.choice.legend{1,stim} = legend_entries{stim};
                                                plot_struct_single{col,tr_pos_con}{p_a,fol}.choice.session_s = strjoin(data_struct.(tr_conditions{tr_pos_con}).(positions{tr_pos_con}{pos}).(stim_con{stim}){col,fol_idx}.session,' - '); % session(s) for this data point
                                                plot_struct_single{col,tr_pos_con}{p_a,fol}.choice.tr_pos_con = plot_groups{tr_pos_con};
                                                plot_struct_single{col,tr_pos_con}{p_a,fol}.choice.tr_con = tr_conditions{tr_pos_con};
                                                plot_struct_single{col,tr_pos_con}{p_a,fol}.choice.G_R_ratio = data_struct.(tr_conditions{tr_pos_con}).(positions{tr_pos_con}{pos}).(stim_con{stim}){col,fol_idx}.G_R_ratios(col);
                                                plot_struct_single{col,tr_pos_con}{p_a,fol}.choice.txt{pos,stim} = sprintf('n=%g \nntrials=%g',...
                                                    data_struct.(tr_conditions{tr_pos_con}).(positions{tr_pos_con}{pos}).(stim_con{stim}){col,fol_idx}.(n{1,p_a}),...
                                                    data_struct.(tr_conditions{tr_pos_con}).(positions{tr_pos_con}{pos}).(stim_con{stim}){col,fol_idx}.(ntrials{ch,1}));
                                            case 'latency'
                                                switch single_or_double_stimulus{s_d}
                                                    case 'double' % put choice as group level
                                                        plot_struct_double{col,tr_pos_con}{p_a,fol}.latency.barvalues{ch,stim} = data_struct.(tr_conditions{tr_pos_con}).(positions{tr_pos_con}{pos}).(stim_con{stim}){col,fol_idx}.(ydata_fname{ch,p_a}); % y value to plot
                                                        plot_struct_double{col,tr_pos_con}{p_a,fol}.latency.errors{ch,stim} = data_struct.(tr_conditions{tr_pos_con}).(positions{tr_pos_con}{pos}).(stim_con{stim}){col,fol_idx}.(SEM_fname{ch,p_a}); % SEM values
                                                        plot_struct_double{col,tr_pos_con}{p_a,fol}.latency.groupnames{ch,1} = choices{ch};
                                                        plot_struct_double{col,tr_pos_con}{p_a,fol}.latency.legend{1,stim} = legend_entries{stim};
                                                        plot_struct_double{col,tr_pos_con}{p_a,fol}.latency.session_s = strjoin(data_struct.(tr_conditions{tr_pos_con}).(positions{tr_pos_con}{pos}).(stim_con{stim}){col,fol_idx}.session,' - '); % session(s) for this data point
                                                        plot_struct_double{col,tr_pos_con}{p_a,fol}.latency.tr_pos_con = plot_groups{tr_pos_con};
                                                        plot_struct_double{col,tr_pos_con}{p_a,fol}.latency.tr_con = tr_conditions{tr_pos_con};
                                                        plot_struct_double{col,tr_pos_con}{p_a,fol}.latency.G_R_ratio = data_struct.(tr_conditions{tr_pos_con}).(positions{tr_pos_con}{pos}).(stim_con{stim}){col,fol_idx}.G_R_ratios(col);
                                                        % plot_struct_double{col,tr_pos_con}{p_a,fol}.latency.txt{ch,stim} = sprintf(txt_nruns_ses{1,p_a},...
                                                        %     data_struct.(tr_conditions{tr_pos_con}).(positions{tr_pos_con}{pos}).(stim_con{stim}){col,fol_idx}.(nruns_ses{1,p_a}),...
                                                        %     data_struct.(tr_conditions{tr_pos_con}).(positions{tr_pos_con}{pos}).(stim_con{stim}){col,fol_idx}.(ntrials{ch,p_a}));
                                                        plot_struct_double{col,tr_pos_con}{p_a,fol}.latency.txt{ch,stim} = num2str(...
                                                            data_struct.(tr_conditions{tr_pos_con}).(positions{tr_pos_con}{pos}).(stim_con{stim}){col,fol_idx}.(ntrials{ch,p_a}));
                                                        plot_struct_double{col,tr_pos_con}{p_a,fol}.latency.latencies{ch,stim} = data_struct.(tr_conditions{tr_pos_con}).(positions{tr_pos_con}{pos}).(stim_con{stim}){col,fol_idx}.(latencies{ch,p_a});
                                                        if strcmp(Settings.bar_plot_latency_avg_mean_or_true_mean,'mean_of_means')
                                                            plot_struct_double{col,tr_pos_con}{p_a,fol}.latency.bar_plot_avg_mean_or_true_mean = 'mean_of_means';
                                                        elseif strcmp(Settings.bar_plot_latency_avg_mean_or_true_mean,'true_mean_of_concatenated_trials')
                                                            plot_struct_double{col,tr_pos_con}{p_a,fol}.latency.bar_plot_avg_mean_or_true_mean = 'true_mean_of_concatenated_trials';
                                                        end
                                                    case 'single' % put pos as group level
                                                        switch tr_conditions{tr_pos_con} %
                                                            case  'single_targ' % plot latency of hits
                                                                plot_struct_single{col,tr_pos_con}{p_a,fol}.latency.barvalues{pos,stim} = data_struct.(tr_conditions{tr_pos_con}).(positions{tr_pos_con}{pos}).(stim_con{stim}){col,fol_idx}.(ydata_fname{ch,p_a}); % y value to plot
                                                                plot_struct_single{col,tr_pos_con}{p_a,fol}.latency.errors{pos,stim} = data_struct.(tr_conditions{tr_pos_con}).(positions{tr_pos_con}{pos}).(stim_con{stim}){col,fol_idx}.(SEM_fname{ch,p_a}); % SEM values
                                                                plot_struct_single{col,tr_pos_con}{p_a,fol}.latency.groupnames{pos,1} = choices{pos};
                                                                plot_struct_single{col,tr_pos_con}{p_a,fol}.latency.legend{1,stim} = legend_entries{stim};
                                                                plot_struct_single{col,tr_pos_con}{p_a,fol}.latency.session_s = strjoin(data_struct.(tr_conditions{tr_pos_con}).(positions{tr_pos_con}{pos}).(stim_con{stim}){col,fol_idx}.session,' - '); % session(s) for this data point
                                                                plot_struct_single{col,tr_pos_con}{p_a,fol}.latency.tr_pos_con = plot_groups{tr_pos_con};
                                                                plot_struct_single{col,tr_pos_con}{p_a,fol}.latency.tr_con = tr_conditions{tr_pos_con};
                                                                plot_struct_single{col,tr_pos_con}{p_a,fol}.latency.G_R_ratio = data_struct.(tr_conditions{tr_pos_con}).(positions{tr_pos_con}{pos}).(stim_con{stim}){col,fol_idx}.G_R_ratios(col);
                                                                % plot_struct_single{col,tr_pos_con}{p_a,fol}.latency.txt{pos,stim} = sprintf(txt_nruns_ses{1,p_a},...
                                                                %     data_struct.(tr_conditions{tr_pos_con}).(positions{tr_pos_con}{pos}).(stim_con{stim}){col,fol_idx}.(nruns_ses{1,p_a}),...
                                                                %     data_struct.(tr_conditions{tr_pos_con}).(positions{tr_pos_con}{pos}).(stim_con{stim}){col,fol_idx}.(ntrials{ch,p_a}));
                                                                plot_struct_single{col,tr_pos_con}{p_a,fol}.latency.txt{pos,stim} = num2str(...
                                                                    data_struct.(tr_conditions{tr_pos_con}).(positions{tr_pos_con}{pos}).(stim_con{stim}){col,fol_idx}.(ntrials{ch,p_a}));
                                                                plot_struct_single{col,tr_pos_con}{p_a,fol}.latency.latencies{pos,stim} = data_struct.(tr_conditions{tr_pos_con}).(positions{tr_pos_con}{pos}).(stim_con{stim}){col,fol_idx}.(latencies{ch,p_a});
                                                            case  'single_distr' % plot latency of incorrect contra and ipsi choices, respectively
                                                                plot_struct_single{col,tr_pos_con}{p_a,fol}.latency.barvalues{pos,stim} = data_struct.(tr_conditions{tr_pos_con}).(positions{tr_pos_con}{pos}).(stim_con{stim}){col,fol_idx}.(ydata_fname{pos,p_a}); % y value to plot
                                                                plot_struct_single{col,tr_pos_con}{p_a,fol}.latency.errors{pos,stim} = data_struct.(tr_conditions{tr_pos_con}).(positions{tr_pos_con}{pos}).(stim_con{stim}){col,fol_idx}.(SEM_fname{pos,p_a}); % SEM values
                                                                plot_struct_single{col,tr_pos_con}{p_a,fol}.latency.groupnames{pos,1} = choices{pos};
                                                                plot_struct_single{col,tr_pos_con}{p_a,fol}.latency.legend{1,stim} =legend_entries{stim};
                                                                plot_struct_single{col,tr_pos_con}{p_a,fol}.latency.session_s = strjoin(data_struct.(tr_conditions{tr_pos_con}).(positions{tr_pos_con}{pos}).(stim_con{stim}){col,fol_idx}.session,' - '); % session(s) for this data point
                                                                plot_struct_single{col,tr_pos_con}{p_a,fol}.latency.tr_pos_con = plot_groups{tr_pos_con};
                                                                plot_struct_single{col,tr_pos_con}{p_a,fol}.latency.tr_con = tr_conditions{tr_pos_con};
                                                                plot_struct_single{col,tr_pos_con}{p_a,fol}.latency.G_R_ratio = data_struct.(tr_conditions{tr_pos_con}).(positions{tr_pos_con}{pos}).(stim_con{stim}){col,fol_idx}.G_R_ratios(col);
                                                                % plot_struct_single{col,tr_pos_con}{p_a,fol}.latency.txt{pos,stim} = sprintf(txt_nruns_ses{1,p_a},...
                                                                %     data_struct.(tr_conditions{tr_pos_con}).(positions{tr_pos_con}{pos}).(stim_con{stim}){col,fol_idx}.(nruns_ses{1,p_a}),...
                                                                %     data_struct.(tr_conditions{tr_pos_con}).(positions{tr_pos_con}{pos}).(stim_con{stim}){col,fol_idx}.(ntrials{pos,p_a}));
                                                                plot_struct_single{col,tr_pos_con}{p_a,fol}.latency.txt{pos,stim} = num2str(...
                                                                    data_struct.(tr_conditions{tr_pos_con}).(positions{tr_pos_con}{pos}).(stim_con{stim}){col,fol_idx}.(ntrials{pos,p_a}));
                                                                plot_struct_single{col,tr_pos_con}{p_a,fol}.latency.latencies{pos,stim} = data_struct.(tr_conditions{tr_pos_con}).(positions{tr_pos_con}{pos}).(stim_con{stim}){col,fol_idx}.(latencies{pos,p_a});
                                                                
                                                        end
                                                        if strcmp(Settings.bar_plot_latency_avg_mean_or_true_mean,'mean_of_means')
                                                            plot_struct_single{col,tr_pos_con}{p_a,fol}.latency.bar_plot_avg_mean_or_true_mean = 'mean_of_means';
                                                        elseif strcmp(Settings.bar_plot_latency_avg_mean_or_true_mean,'true_mean_of_concatenated_trials')
                                                            plot_struct_single{col,tr_pos_con}{p_a,fol}.latency.bar_plot_avg_mean_or_true_mean = 'true_mean_of_concatenated_trials';
                                                        end
                                                end
                                        end
                                    end
                                end
                            else
                                warning('Some trial conditions may not be plotted. Change ''Settings.analyze_trial_conditions'' to plot data of all trial conditions')
                            end
                        end
                    end
                end
            end
        end
    end
end
end

function plot_bargraph (Settings, plot_struct,stat_struct_ttest)

for col = 1:size(plot_struct,1)
    for tr_pos_con = 1:size(plot_struct,2)
        if ~isempty(plot_struct{col,tr_pos_con})
            per_or_across = size(plot_struct{col,tr_pos_con},1); % how many rows? 1st row = per session(one column per session); 2nd row = across sessions (one column only)
            for p_a = 1:per_or_across
                num_sessions = size([plot_struct{col,tr_pos_con}{p_a,:}],2); % how many columns = sessions
                for fol = 1:num_sessions;
                    fnames_c_l = fieldnames(plot_struct{col,tr_pos_con}{p_a,fol});
                    for c_l = 1:numel(fnames_c_l)
                        fig = figure;
                        tmp_data = plot_struct{col,tr_pos_con}{p_a,fol}.(fnames_c_l{c_l});
                        % nan as double -> cell2mat doesn't work ... KK
                        % added the following loop
                        for ind = 1: numel(tmp_data.barvalues)
                        if ~strcmp(class(tmp_data.barvalues{ind}), 'single')
                            tmp_data.barvalues{ind} = single(tmp_data.barvalues{ind}); 
                        end
                        end
                         for ind = 1: numel(tmp_data.errors)
                        if ~strcmp(class( tmp_data.errors{ind}), 'single')
                            tmp_data.errors{ind} = single(tmp_data.errors{ind}); 
                        end
                        end
                        if any(any(~isnan(cell2mat(tmp_data.barvalues)))) % if not all barvalues are nan, i.e. at least one value is not nan
                            
                            %%% TITLE, LABELS, AXES RANGES %%%
                            title_pses = get_title_session_part(tmp_data.session_s, Settings);
                            title_ptrcon = tmp_data.tr_pos_con;
                            title_pcol = ['G-R-ratio ',num2str(tmp_data.G_R_ratio)];
                            
                            if Settings.statistics_on && Settings.ttests && Settings.FDR % ttests and false discovery rate correction for multiple testing
                                title_pstat = ['ttest-',num2str(Settings.alpha),'-FDR'];
                            elseif Settings.statistics_on && Settings.ttests && ~Settings.FDR
                                title_pstat = ['ttest-',num2str(Settings.alpha)];
                            else
                                title_pstat = '';
                            end
                            
                            if strcmp(fnames_c_l{c_l},'choice') % if choice is to be plotted
                                title_pc_l = 'Choice';
                                if strcmp(tmp_data.tr_con,'single_targ') % single-targ trials
                                    bw_ylabel = 'Hitrate / Choice proportion';
                                    bw_xlabel = 'Target position';
                                elseif strcmp(tmp_data.tr_con,'single_distr') % single-distr trials
                                    bw_ylabel = 'Hitrate / Fixation choice proportion';
                                    bw_xlabel = 'Distractor position';
                                else % double stimuli trials
                                    bw_ylabel = 'Choice proportion';
                                    bw_xlabel = 'Choice position';
                                    title_pc_l = 'Choice';
                                end
                                if strcmp(tmp_data.bar_plot_mean_or_true_proportion,'true_proportion')
                                    title_pprop = 'true_prop';
                                elseif strcmp(tmp_data.bar_plot_mean_or_true_proportion,'mean_and_SEM')
                                    title_pprop = 'mean-SEM';
                                end
                                
                            elseif strcmp(fnames_c_l{c_l},'latency') % if latency is to be plotted
                                bw_ylabel = 'Latency';
                                bw_xlabel = 'Choice position';
                                title_pc_l = bw_ylabel;
                                if strcmp(tmp_data.bar_plot_avg_mean_or_true_mean,'true_mean_of_concatenated_trials')
                                    title_pprop = 'true_mean';
                                elseif strcmp(tmp_data.bar_plot_avg_mean_or_true_mean,'mean_of_means')
                                    title_pprop = 'mean-of-means';
                                end
                            end
                            
                              if strcmp(Settings.hemisphere_of_stimulation,'left')
                                    tmp_data.groupnames = {'ipsi', 'contra'};
                                elseif  strcmp(Settings.hemisphere_of_stimulation,'right')
                                    tmp_data.groupnames = {'contra', 'ipsi'};
                                    
                             end
                            if length(title_pses) > 10
                                if isempty(title_pstat)
                                     bw_title = [title_pc_l,' - ',title_pprop,' - ',title_ptrcon,' - ',title_pcol,sprintf('\n'),title_pses]; % two lines
                                    file_name = [title_pc_l,' - ',title_pprop,' - ',title_ptrcon,' - ',title_pcol]; % one line %title_pses,' - ',

                                else
                                       bw_title = [title_pc_l,' - ',title_pprop,' - ',title_ptrcon,' - ',title_pcol,sprintf('\n'),title_pses,' - ',title_pstat]; %bw_title = []; % two lines
                                    file_name = [title_pc_l,' - ',title_pprop,' - ',title_ptrcon,' - ',title_pcol,' - ',title_pstat]; % one line title_pses,' - ',
                                end
                            else
                                if isempty(title_pstat)
                                    bw_title = [title_pc_l,' - ',title_pprop,' - ',title_pses,' - ',title_ptrcon,' - ',title_pcol]; % one line
                                    file_name = bw_title;
                                else
                                    bw_title = [title_pc_l,' - ',title_pprop,' - ',title_pses,' - ',title_ptrcon,' - ',title_pcol,' - ',title_pstat]; % one line
                                    file_name = bw_title;
                                end
                            end
                            
                            %%% BARWEB CALL & PLOTTING %%%
                            bw_colormap = [100,100,100;96,196,222;195,210,0;239,125,21]./255; gridstatus = [];
                            [handles,xpos,yval,err] = barweb(cell2mat(tmp_data.barvalues), cell2mat(tmp_data.errors), [], tmp_data.groupnames, bw_title, bw_xlabel, bw_ylabel, bw_colormap, gridstatus, tmp_data.legend);
                            set(handles.title,'interpreter','none')
                            set(handles.legend,'interpreter','none')
                            
                            %%% DATAPOINT TEXTS %%%
                            if strcmp(fnames_c_l{c_l},'choice')
                                ylim_min = 0;
                                ylim_max = 1;
                                txt_ypos_max = ylim_max - 0.15;
                                txt_ypos_min = ylim_min;
                            elseif strcmp(fnames_c_l{c_l},'latency')
                                ylim_min = 0; % 0.15
                                ylim_max = 0.4; % 0.6
                                txt_ypos_max = ylim_max - 0.05;
                                txt_ypos_min = ylim_min;
                            end
                            ylim([ylim_min ylim_max]);
                            % set(gca, 'YTick', 0:0.05:0.4); % change for latency plots
                            yval = reshape(yval,1,size(yval,1)*size(yval,2));
                            err = reshape(err,1,size(err,1)*size(err,2));
                            ypos_above = (yval+err)*1.05;
                            ypos_below = (yval-err)*0.90;
                            ypos = ypos_below; % ypos = ypos_above;
                            ypos(isnan(ypos)) = 0; % replaces NaNs with zeros in order to plot the txt for latency categories with no trials; comment this out and no txt will be plotted for these categories
                            % ypos(ypos > txt_ypos_max) = ypos_below(ypos_above > txt_ypos_max); % if text would be placed too high to fit on the plot, plot it below the data point
                            ypos(ypos <= txt_ypos_min) = txt_ypos_min+((ylim_max-ylim_min)*0.05); % if text would be placed too low to fit on the plot, plot it ipsi above ylim_min
                            
                            
                            if  strcmp(tmp_data.bar_plot_avg_mean_or_true_mean,'mean_and_SEM')
                                for ind = 1: size(tmp_data.txt,1)
                                    for in = 1: size(tmp_data.txt,2)
                                a = tmp_data.txt{ind,in};
                                tmp_data.txt{ind,in} = a(1:3); 
                                    end
                                end
                            text(xpos,double(ypos),tmp_data.txt,'HorizontalAlignment','center','VerticalAlignment','bottom');
                            else 
                             text(xpos,double(ypos),tmp_data.txt,'HorizontalAlignment','center','VerticalAlignment','bottom');

                            end
                            
                            %%% SIGNIFICANCE STARS IN PLOTS %%%
                            if Settings.ttests
                                 comparisons = {[1 2]} ;%{[1 2] [1 3] [1 4] [2 3] [2 4] [3 4]}; % combinations of comparisons among one quartett of stimulation conditions
                                all_x_values = reshape(xpos,length(xpos)/length(xpos),length(xpos))';
                               % comparisons = {[1 2] [1 3] [1 4] [2 3] [2 4] [3 4]}; % combinations of comparisons among one quartett of stimulation conditions
                                %all_x_values = reshape(xpos,length(xpos)/4,4)'; % x-values of bars, each column contains one group of 4 stimulation bars
                                if isfield(stat_struct_ttest,tmp_data.tr_con)
                                    if strcmp(tmp_data.tr_pos_con,'contra_targ-ipsi_distr')
                                        position = {'contra'}; % only get tests for the contra target condition
                                    elseif strcmp(tmp_data.tr_pos_con,'contra_distr-ipsi_targ')
                                        position = {'ipsi'}; % only get tests for the ipsi target condition
                                    else
                                        position = fieldnames(stat_struct_ttest.(tmp_data.tr_con));
                                    end
                                    for pos = 1:numel(position)
                                        % choices = fieldnames(stat_struct_ttest.(tmp_data.tr_con).(position{pos}));
                                        c_l_fnames = fieldnames(stat_struct_ttest.(tmp_data.tr_con).(position{pos}));
                                        c_l_idx = cellfun(@(x) ~isempty(strfind(x,fnames_c_l{c_l})),c_l_fnames,'UniformOutput',1); % use only those fields that contain 'choice' or 'latency', respectively
                                        c_l_fnames = c_l_fnames(c_l_idx);
                                        % only get tests of choices towards the stimulus; don't include tests on fixation choices
                                        if strcmp(fnames_c_l{c_l},'choice') && strcmp(tmp_data.tr_pos_con,'single_targ') || strcmp(fnames_c_l{c_l},'choice') && strcmp(tmp_data.tr_pos_con,'single_distr')
                                            for ch = pos % only contra_choices for contra stimuli, only ipsi_choices for ipsi stimuli
                                                %elseif strcmp(fnames_c_l{c_l},'latency')
                                                if Settings.FDR
                                                    sig_comp_idx = logical([stat_struct_ttest.(tmp_data.tr_con).(position{pos}).(c_l_fnames{ch}){col,1}{:,9}]);
                                                else
                                                    sig_comp_idx = logical([stat_struct_ttest.(tmp_data.tr_con).(position{pos}).(c_l_fnames{ch}){col,1}{:,5}]);
                                                end
                                                if sum(sig_comp_idx > 0) % if at least one t-test was significant
                                                    groups = comparisons(sig_comp_idx); % find all significant ttests by column 5
                                                    sig_comp = [comparisons(sig_comp_idx)];
                                                    for comp = 1:numel(sig_comp)
                                                        groups{comp}(1) = all_x_values(sig_comp{comp}(1),ch); % combinations of xpos values that define those pairs of bars
                                                        groups{comp}(2) = all_x_values(sig_comp{comp}(2),ch); % that are significantly different
                                                    end
                                                    if Settings.FDR
                                                        all_p_values = [stat_struct_ttest.(tmp_data.tr_con).(position{pos}).(c_l_fnames{ch}){col,1}{:,10}];
                                                    else
                                                        all_p_values = [stat_struct_ttest.(tmp_data.tr_con).(position{pos}).(c_l_fnames{ch}){col,1}{:,6}];
                                                    end
                                                    stats = all_p_values(sig_comp_idx);
                                                    starlevels = [0.05, 1E-2, 1E-3];
                                                    linewidth = 1.5;
                                                    fontsize = 14;
                                                    
                                                    sigstar_UZ(groups,stats,0,starlevels,linewidth,fontsize); % call sigstar
                                                end
                                            end
                                        else
                                            for ch = 1:numel(c_l_fnames)
                                                
                                                if ch == 1
                                                     comparisons = {[1 3]} ;%{[1 2] [1 3] [1 4] [2 3] [2 4] [3 4]}; % combinations of comparisons among one quartett of stimulation conditions
                                                     all_x_values = reshape(xpos,length(xpos)/length(xpos),length(xpos))';
                                                elseif ch == 2
                                                     comparisons = {[2 4]} ;%{[1 2] [1 3] [1 4] [2 3] [2 4] [3 4]}; % combinations of comparisons among one quartett of stimulation conditions
                                                     all_x_values = reshape(xpos,length(xpos)/length(xpos),length(xpos))';
                                                end
                                                
                                                 
                                                if Settings.FDR
                                                    sig_comp_idx = logical([stat_struct_ttest.(tmp_data.tr_con).(position{pos}).(c_l_fnames{ch}){col,1}{:,9}]);
                                                else
                                                    sig_comp_idx = logical([stat_struct_ttest.(tmp_data.tr_con).(position{pos}).(c_l_fnames{ch}){col,1}{:,5}]);
                                                end
                                                if sum(sig_comp_idx > 0) % if at least one t-test was significant
                                                    groups = comparisons(sig_comp_idx); % find all significant ttests by column 5
                                                    sig_comp = [comparisons(sig_comp_idx)];
                                                    for comp = 1:numel(sig_comp)
                                                        groups{comp}(1) = all_x_values(sig_comp{comp}(1),1); %all_x_values(sig_comp{comp}(1),ch); combinations of xpos values that define those pairs of bars
                                                        groups{comp}(2) = all_x_values(sig_comp{comp}(2),1); %all_x_values(sig_comp{comp}(2),ch); that are significantly different
                                                    end
                                                    if Settings.FDR
                                                        all_p_values = [stat_struct_ttest.(tmp_data.tr_con).(position{pos}).(c_l_fnames{ch}){col,1}{:,10}];
                                                    else
                                                        all_p_values = [stat_struct_ttest.(tmp_data.tr_con).(position{pos}).(c_l_fnames{ch}){col,1}{:,6}];
                                                    end
                                                    stats = all_p_values(sig_comp_idx);
                                                    starlevels = [0.05, 1E-2, 1E-3];
                                                    linewidth = 1.5;
                                                    fontsize = 14;
                                                    
                                                    sigstar_UZ(groups,stats,0,starlevels,linewidth,fontsize); % call sigstar
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                            
                            %%% SAVE PLOT %%%
                            if Settings.save_figures
                                save_dir = get_save_dir(Settings.folders{fol}, Settings);
                                save_figure_as(fig, Settings, file_name, {'.ai','.png'}, save_dir) % save
                                if Settings.close_after_saving
                                    close(gcf)
                                end
                            end
                        else
                            close(gcf)
                        end
                    end
                end
            end
        end
    end
end
end

function [plot_struct_double_II, plot_struct_single_II] = getbarplotstruct_II (Settings,plot_struct_double_II,plot_struct_single_II,data_struct,p_a,stim_con)
num_stim_con = numel(stim_con);
if  strcmp(Settings.Experiment , 'Inactivation') 
    legend_entries = {'Pre','Post'};
elseif  strcmp(Settings.Experiment , 'Microstimulation') 
legend_entries = {'control','-250 ms','-100 ms','+50 ms'};
end
%legend_entries = {'control','-80 ms','at GO','+80 ms'};
single_or_double_stimulus = {'double','single'};
choice_or_latency = {'choice','latency'};
plot_choice_or_latency = [Settings.analyze_choice,Settings.analyze_latency];
nruns_ses = {'num_runs','num_sessions'};
for s_d = 1:numel(single_or_double_stimulus)
    switch single_or_double_stimulus{s_d} % is this going to be a plot of trials with two stimuli or a single stimulus presented
        case 'double' % double stimulus conditions
            
            if strcmp(Settings.hemisphere_of_stimulation,'left')
                plot_groups = {{'targ_targ'},{'ipsi_targ-contra_distr','ipsi_distr-contra_targ','distr_distr'}};
                
            elseif strcmp(Settings.hemisphere_of_stimulation,'right')
                plot_groups = {{'targ_targ'},{'contra_targ-ipsi_distr','contra_distr-ipsi_targ','distr_distr'}};
                
            end
            tr_conditions = {{'targ_targ'},{'targ_distr','targ_distr','distr_distr'}};
            positions = {{'combined'},Settings.analyze_by_position}; % to find the correct field in the data structure
            
        case 'single' % single stimulus conditions
            
            plot_groups = {{'single_targ','single_targ'},{'single_distr','single_distr'}};
            tr_conditions = plot_groups;
            positions = {Settings.analyze_by_position(1:2),Settings.analyze_by_position(1:2)};
            
    end
    for c_l = 1:numel(plot_choice_or_latency)
        if plot_choice_or_latency(c_l)
            for plt = 1:numel(plot_groups) % for each separate plots as defined in plot_groups
                for tr_pos_con = 1:numel(plot_groups{plt}) % for each separate plots as defined in plot_groups
                    switch choice_or_latency{c_l}
                        case 'choice'
                            switch single_or_double_stimulus{s_d}
                                case 'double' % double stimulus conditions
                                    choices = {'contra' 'fixation' 'ipsi'};
                                    if strcmp(Settings.bar_plot_choice_mean_or_true_proportion,'mean_and_SEM') % plot the mean and SEM across runs (plot per session) or across sessions (plot across sessions)
                                        ydata_fname = {'mean_prop_choice_across_runs_contra','mean_prop_choice_across_sessions_contra';...
                                            'mean_prop_choice_across_runs_fixation','mean_prop_choice_across_sessions_fixation';...
                                            'mean_prop_choice_across_runs_ipsi','mean_prop_choice_across_sessions_ipsi'};
                                    elseif strcmp(Settings.bar_plot_choice_mean_or_true_proportion,'true_proportion') % plot the true proportions per session or across sessions
                                        ydata_fname = {'true_prop_choice_per_session_contra','true_prop_choice_across_sessions_contra';...
                                            'true_prop_choice_per_session_fixation','true_prop_choice_across_sessions_fixation';...
                                            'true_prop_choice_per_session_ipsi','true_prop_choice_across_sessions_ipsi'};
                                    end
                                    SEM_fname = {'SEM_choice_contra';...
                                        'SEM_choice_fixation';...
                                        'SEM_choice_ipsi'};
                                    ntrials = {'num_choice_contra';...
                                        'num_choice_fixation';...
                                        'num_choice_ipsi'};
                                case 'single' % single stimulus conditions
                                    if strcmp(Settings.bar_plot_choice_mean_or_true_proportion,'mean_and_SEM') % plot the mean and SEM across runs (plot per session) or across sessions (plot across sessions)
                                        ydata_fname = {'mean_prop_choice_across_runs_','mean_prop_choice_across_sessions_';... % contra or ipsi saccades, respectively to positions
                                            'mean_prop_choice_across_runs_','mean_prop_choice_across_sessions_'}; % fixation choices
                                    elseif strcmp(Settings.bar_plot_choice_mean_or_true_proportion,'true_proportion') % plot the true proportions per session or across sessions
                                        ydata_fname = {'true_prop_choice_per_session_','true_prop_choice_across_sessions_';...
                                            'true_prop_choice_per_session_','true_prop_choice_across_sessions_'};
                                    end
                                    SEM_fname = 'SEM_choice_';
                                    ntrials = 'num_choice_';
                                    choices_single = {'saccade_choice','fixation'};
                            end
                            
                    end
                    % find the sessions to plot
                    if isempty(Settings.load_data)
                        sessions_idx = 1:numel(Settings.folders); % if no data is loaded all sessions in the settings are plotted
                    else
                        tr_condition = fieldnames(data_struct);
                        side = fieldnames(data_struct.(tr_condition{1}));
                        stimulation = fieldnames(data_struct.(tr_condition{1}).(side{1}));
                        array = [data_struct.(tr_condition{1}).(side{1}).(stimulation{1}){1,:}]; %
                        all_loaded_sessions = [array.session];
                        [~,sessions_idx] = ismember(Settings.folders,all_loaded_sessions); % get indices of loaded sessions that are specified in Settings
                    end
                    for fol = 1:numel(sessions_idx)
                        fol_idx = sessions_idx(fol);
                        if p_a == 2 && fol > 1 % in data_struct_across_sessions there is only one column -> fol needs to be 1
                            continue;
                        end
                        %                     pos = tr_pos_con;
                        for stim = 1:num_stim_con
                            if ismember(tr_conditions{plt}{tr_pos_con},fieldnames(data_struct)) % if the specified trial conditions are actually analyzed and field of data_struct
                                for col = 1:size([data_struct.(tr_conditions{plt}{tr_pos_con}).(positions{plt}{tr_pos_con}).(stim_con{stim}){:,fol_idx}],2); % for all colors (rows in data_struct) in this trial condition
                                    for ch = 1:size(ydata_fname,1)
                                        tmp_data = data_struct.(tr_conditions{plt}{tr_pos_con}).(positions{plt}{tr_pos_con}).(stim_con{stim}){col,fol_idx};
                                        switch choice_or_latency{c_l}
                                            case 'choice' % put choice as group level
                                                switch single_or_double_stimulus{s_d}
                                                    case 'double' % double stimulus conditions
                                                        plot_struct_double_II{col,plt}{p_a,fol}.(choices{ch}).barvalues{tr_pos_con,stim} = tmp_data.(ydata_fname{ch,p_a}); % y value to plot
                                                        if strcmp(Settings.bar_plot_choice_mean_or_true_proportion,'true_proportion') % plot the true proportions per session or across sessions
                                                            plot_struct_double_II{col,plt}{p_a,fol}.(choices{ch}).errors{tr_pos_con,stim} = 0; % SEM value is 0 because no SEM can be calculated from only one proportion
                                                            plot_struct_double_II{col,plt}{p_a,fol}.(choices{ch}).bar_plot_mean_or_true_proportion = 'true_proportion';
                                                        else % plot the mean and SEM
                                                            plot_struct_double_II{col,plt}{p_a,fol}.(choices{ch}).errors{tr_pos_con,stim} = tmp_data.(SEM_fname{ch}); % SEM values
                                                            plot_struct_double_II{col,plt}{p_a,fol}.(choices{ch}).bar_plot_mean_or_true_proportion = 'mean_and_SEM';
                                                        end
                                                        plot_struct_double_II{col,plt}{p_a,fol}.(choices{ch}).groupnames{tr_pos_con,1} = plot_groups{plt}{tr_pos_con};
                                                        plot_struct_double_II{col,plt}{p_a,fol}.(choices{ch}).position{tr_pos_con} = positions{plt}{tr_pos_con};
                                                        plot_struct_double_II{col,plt}{p_a,fol}.(choices{ch}).legend{1,stim} = legend_entries{stim};
                                                        plot_struct_double_II{col,plt}{p_a,fol}.(choices{ch}).session_s = strjoin(tmp_data.session,' - '); % session(s) for this data point
                                                        plot_struct_double_II{col,plt}{p_a,fol}.(choices{ch}).tr_pos_con{tr_pos_con} = plot_groups{plt}{tr_pos_con};
                                                        plot_struct_double_II{col,plt}{p_a,fol}.(choices{ch}).tr_con{tr_pos_con} = tr_conditions{plt}{tr_pos_con};
                                                        plot_struct_double_II{col,plt}{p_a,fol}.(choices{ch}).G_R_ratio = tmp_data.G_R_ratios(col);
                                                        %test in the plots
                                                        plot_struct_double_II{col,plt}{p_a,fol}.(choices{ch}).txt{tr_pos_con,stim} = sprintf('n=%g \nntrials=%g',...
                                                            tmp_data.(nruns_ses{p_a}),...
                                                            tmp_data.(ntrials{ch,1}));
                                                        plot_struct_double_II{col,plt}{p_a,fol}.(choices{ch}).num_total_trials(tr_pos_con,stim) = tmp_data.num_total_trials; % get number of total trials
                                                    case 'single' % single stimulus conditions
                                                        choices = {positions{plt}{tr_pos_con} 'fixation'}; % take the saccade option and fixation as two choice options
                                                        ydata_fname_single = [ydata_fname{ch,p_a},choices{ch}];
                                                        SEM_fname_single = [SEM_fname,choices{ch}];
                                                        ntrials_single = [ntrials,choices{ch}];
                                                        plot_struct_single_II{col,plt}{p_a,fol}.(choices_single{ch}).barvalues{tr_pos_con,stim} = tmp_data.(ydata_fname_single); % y value to plot
                                                        if strcmp(Settings.bar_plot_choice_mean_or_true_proportion,'true_proportion') % plot the true proportions per session or across sessions
                                                            plot_struct_single_II{col,plt}{p_a,fol}.(choices_single{ch}).errors{tr_pos_con,stim} = 0; % SEM value is 0 because no SEM can be calculated from only one proportion
                                                            plot_struct_single_II{col,plt}{p_a,fol}.(choices_single{ch}).bar_plot_mean_or_true_proportion = 'true_proportion';
                                                        else % plot the mean and SEM
                                                            plot_struct_single_II{col,plt}{p_a,fol}.(choices_single{ch}).errors{tr_pos_con,stim} = tmp_data.(SEM_fname_single); % SEM values
                                                            plot_struct_single_II{col,plt}{p_a,fol}.(choices_single{ch}).bar_plot_mean_or_true_proportion = 'mean_and_SEM';
                                                        end
                                                        plot_struct_single_II{col,plt}{p_a,fol}.(choices_single{ch}).groupnames{tr_pos_con,1} = positions{plt}{tr_pos_con};
                                                        plot_struct_single_II{col,plt}{p_a,fol}.(choices_single{ch}).position{tr_pos_con} = positions{plt}{tr_pos_con};
                                                        plot_struct_single_II{col,plt}{p_a,fol}.(choices_single{ch}).legend{1,stim} = legend_entries{stim};
                                                        plot_struct_single_II{col,plt}{p_a,fol}.(choices_single{ch}).session_s = strjoin(tmp_data.session,' - '); % session(s) for this data point
                                                        plot_struct_single_II{col,plt}{p_a,fol}.(choices_single{ch}).tr_pos_con{1} = plot_groups{plt}{tr_pos_con};
                                                        plot_struct_single_II{col,plt}{p_a,fol}.(choices_single{ch}).tr_con{1} = tr_conditions{plt}{tr_pos_con};
                                                        plot_struct_single_II{col,plt}{p_a,fol}.(choices_single{ch}).G_R_ratio = tmp_data.G_R_ratios(col);
                                                        % text
                                                        plot_struct_single_II{col,plt}{p_a,fol}.(choices_single{ch}).txt{tr_pos_con,stim} = sprintf('n=%g \nntrials=%g',...
                                                            tmp_data.(nruns_ses{p_a}),...
                                                            tmp_data.(ntrials_single));
                                                        plot_struct_single_II{col,plt}{p_a,fol}.(choices_single{ch}).num_total_trials(tr_pos_con,stim) = tmp_data.num_total_trials; % get number of total trials
                                                end
                                        end
                                    end
                                end
                            else
                                warning('Some trial conditions may not be plotted. Change ''Settings.analyze_trial_conditions'' to plot data of all trial conditions')
                            end
                        end
                    end
                end
            end
        end
    end
end
end

function plot_bargraph_II (Settings, plot_struct,stat_struct_ttest)
for col = 1:size(plot_struct,1)
    for plt = 1:size(plot_struct,2)
        if ~isempty(plot_struct{col,plt})
            per_or_across = size(plot_struct{col,plt},1); % how many rows? 1st row = per session(one column per session); 2nd row = across sessions (one column only)
            for p_a = 1:per_or_across
                num_sessions = size([plot_struct{col,plt}{p_a,:}],2); % how many columns = sessions
                for fol = 1:num_sessions;
                    choices = fieldnames(plot_struct{col,plt}{p_a,fol});
                    for ch = 1:numel(choices)
                        fig = figure;
                        tmp_data = plot_struct{col,plt}{p_a,fol}.(choices{ch});
                        if any(any(~isnan(cell2mat(tmp_data.barvalues)))) % if not all barvalues are nan, i.e. at least one value is not nan
                            
                            %%% TITLE, LABELS, AXES RANGES %%%
                            % Change default axes fonts.
                            set(0,'DefaultAxesFontName', 'Arial')
                            set(0,'DefaultAxesFontSize', 12)
                            
                            % Change default text fonts.
                            set(0,'DefaultTextFontname', 'Arial')
                            set(0,'DefaultTextFontSize', 12)
                            
                            title_pses = get_title_session_part(tmp_data.session_s, Settings);
                            if isequal(tmp_data.tr_pos_con,{'contra_targ-ipsi_distr','contra_distr-ipsi_targ','distr_distr'})|| isequal(tmp_data.tr_pos_con,{'ipsi_targ-contra_distr','ipsi_distr-contra_targ','distr_distr'})  % double stimulus with distractors trials
                                title_ptrcon = 'targ_distr-distr_distr';
                            else
                                title_ptrcon = strjoin(tmp_data.tr_pos_con,' - ');
                            end
                            
                            if isequal(tmp_data.tr_pos_con,{'single_targ'})|| isequal(tmp_data.tr_pos_con,{'single_distr'})
                                if strcmp(Settings.hemisphere_of_stimulation,'left')
                                    tmp_data.groupnames = {'ipsi', 'contra'};
                                elseif  strcmp(Settings.hemisphere_of_stimulation,'ipsi')
                                    tmp_data.groupnames = {'contra', 'right'};
                                    
                                end
                            end
                            
                            
                            title_pcol = ['G-R-ratio ',num2str(tmp_data.G_R_ratio)];
                            
                            if Settings.statistics_on && Settings.ttests && Settings.FDR % ttests and false discovery rate correction for multiple testing
                                title_pstat = ['ttest-',num2str(Settings.alpha),'-FDR'];
                            elseif Settings.statistics_on && Settings.ttests && ~Settings.FDR
                                title_pstat = ['ttest-',num2str(Settings.alpha)];
                            else
                                title_pstat = '';
                            end
                            
                            if strcmp(choices{ch},'contra') % double stimulus trials, if contra choice is to be plotted
                                title_ch  = Settings.GraphProperties.contra.title_ch;
                                bw_ylabel = Settings.GraphProperties.contra.bw_ylabel ;
                            elseif strcmp(choices{ch},'ipsi') % double stimulus trials, ipsi choice
                                title_ch = Settings.GraphProperties.ipsi.title_ch;
                                bw_ylabel = Settings.GraphProperties.ipsi.bw_ylabel ;
                            elseif strcmp(choices{ch},'fixation') % double & single stimulus trials, fixation choice
                                title_ch = 'Fixation_choices';
                                if isequal(tmp_data.tr_con,{'single_distr'}) % single_distr trials
                                    bw_ylabel = 'Fixation Selection (correct)';
                                elseif isequal(tmp_data.tr_con,{'single_targ'}) % single_targ trials
                                    bw_ylabel = 'Fixation Selection (error)';
                                else % double stimulus trials
                                    bw_ylabel = 'Fixation Selection';
                                end
                            elseif strcmp(choices{ch},'saccade_choice') % if single stimulus trials saccade choices is to be plotted
                                title_ch = 'Saccade_choices';
                                if isequal(tmp_data.tr_con,{'single_targ'}) % single_targ trials
                                    bw_ylabel = 'Target Selection (correct)';
                                else
                                    bw_ylabel = 'Distractor Selection (error)';
                                end
                            end
                            bw_xlabel = 'Trial condition';
                            if strcmp(tmp_data.bar_plot_mean_or_true_proportion,'true_proportion')
                                title_pprop = 'true_prop';
                            elseif strcmp(tmp_data.bar_plot_mean_or_true_proportion,'mean_and_SEM')
                                title_pprop = 'mean-SEM';
                            end
                            
                            if length(title_pses) > 10 %KK
                                if isempty(title_pstat)
                                    bw_title = []; %[title_ch,' - ',title_pprop,' - ',title_ptrcon,' - ',title_pcol,sprintf('\n'), title_pses]; % two lines 
                                    %                                     file_name = ['Choice - no_txt - ',title_pprop,' - ',title_pses,' - ',title_ptrcon,' - ',title_ch,' - ',title_pcol]; % one line
                                    file_name = ['Choice-',title_pprop,'-',title_ptrcon,'-',title_ch,'-',title_pcol]; % one line %title_pses,'-'
                                    
                                else
                                   bw_title = []; % bw_title = [title_ch,' - ',title_pprop,' - ',title_ptrcon,' - ',title_pcol,sprintf('\n'),title_pses,' - ',title_pstat]; % two lines
                                    %                                     file_name = ['Choice - no_txt - ',title_pprop,' - ',title_pses,' - ',title_ptrcon,' - ',title_ch,' - ',title_pcol,' - ',title_pstat]; % one line
                                    file_name = ['Choice-',title_pprop,'-',title_ptrcon,'-',title_ch,'-',title_pcol,'-',title_pstat]; % one line % title_pses,'-',
                                end
                            else
                                if isempty(title_pstat)
                                    bw_title = [title_ch,' - ',title_pprop,' - ',title_pses,' - ',title_ptrcon,' - ',title_pcol]; % one line
                                    %                                     filename = ['Choice - no_txt - ',title_pprop,' - ',title_pses,' - ',title_ptrcon,' - ',title_ch,' - ',title_pcol]; % one line
                                    file_name = ['Choice-',title_pprop,'-',title_pses,'-',title_ptrcon,'-',title_ch,'-',title_pcol]; % one line
                                else
                                    bw_title = [title_ch,' - ',title_pprop,' - ',title_pses,' - ',title_ptrcon,' - ',title_pcol,' - ',title_pstat]; % one line
                                    %                                     file_name = ['Choice - no_txt - ',title_pprop,' - ',title_pses,' - ',title_ptrcon,' - ',title_ch,' - ',title_pcol,' - ',title_pstat]; % one line
                                    file_name = ['Choice-',title_pprop,'-',title_pses,'-',title_ptrcon,'-',title_ch,'-',title_pcol,'-',title_pstat]; % one line
                                end
                            end
                            
                            %%% BARWEB CALL & PLOTTING %%%
                            bw_colormap = [100,100,100;96,196,222;195,210,0;239,125,21]./255; gridstatus = [];
                            [handles,xpos,yval,err] = barweb(cell2mat(tmp_data.barvalues), cell2mat(tmp_data.errors), [], tmp_data.groupnames, bw_title, bw_xlabel, bw_ylabel, bw_colormap, gridstatus, tmp_data.legend);
                            set(handles.title,'interpreter','none')
                            set(handles.legend,'interpreter','none')
                            set(handles.ca,'TicklabelInterpreter','none') %% KK: changed it from comment to code
                            
                            %%% DATAPOINT TEXTS %%%
                            ylim_min                = 0;
                            ylim_max                = 1;
                            txt_ypos_max            = ylim_max - 0.15;
                            txt_ypos_min            = ylim_min;
                            ylim([ylim_min ylim_max]);
                            yval                    = reshape(yval,1,size(yval,1)*size(yval,2));
                            err                     = reshape(err,1,size(err,1)*size(err,2));
                            ypos_above              = (yval+err)*1.05;
                            ypos_below              = (yval-err)*0.90;
                            ypos                    = ypos_above;
                            ypos(isnan(ypos))       = 0; % replaces NaNs with zeros in order to plot the txt for latency categories with no trials; comment this out and no txt will be plotted for these categories
                            ypos(ypos > txt_ypos_max) = ypos_below(ypos_above > txt_ypos_max); % if text would be placed too high to fit on the plot, plot it below the data point
                            ypos(ypos <= txt_ypos_min) = txt_ypos_min+((ylim_max-ylim_min)*0.05); % if text would be placed too low to fit on the plot, plot it ipsi above ylim_min
                            % adds the
                            if  strcmp(tmp_data.bar_plot_mean_or_true_proportion,'mean_and_SEM')
                                for ind = 1: size(tmp_data.txt,1)
                                    for in = 1: size(tmp_data.txt,2)
                                a = tmp_data.txt{ind,in};
                                tmp_data.txt{ind,in} = a(1:3); 
                                    end
                                end
                            text(xpos,double(ypos),tmp_data.txt,'HorizontalAlignment','center','VerticalAlignment','bottom');
                            else 
                             text(xpos,double(ypos),tmp_data.txt,'HorizontalAlignment','center','VerticalAlignment','bottom');

                            end
                            %add the text for number of trials
                            text(1,0.5,[num2str(min(min(tmp_data.num_total_trials))),' - ',num2str(max(max(tmp_data.num_total_trials)))]);
                            %%% SIGNIFICANCE STARS IN PLOTS %%%
                            if Settings.ttests
                                comparisons = {[1 2]} ;%{[1 2] [1 3] [1 4] [2 3] [2 4] [3 4]}; % combinations of comparisons among one quartett of stimulation conditions
                                all_x_values = reshape(xpos,length(xpos)/length(xpos),length(xpos))'; % x-values of bars, each column contains one group of 4 stimulation bars
                                for tr_con = 1:numel(tmp_data.tr_con)
                                    if isfield(stat_struct_ttest,tmp_data.tr_con{tr_con})
                                        if isequal(tmp_data.tr_pos_con,{'single_targ'}) || isequal(tmp_data.tr_pos_con,{'single_distr'})
                                            num_pos = numel(tmp_data.position);
                                            i = 1;
                                        else
                                            num_pos = tr_con;
                                            i = tr_con;
                                        end
                                        for pos = i:num_pos
                                            % choices = fieldnames(stat_struct_ttest.(tmp_data.tr_con).(position{pos}));
                                            c_l_fnames = fieldnames(stat_struct_ttest.(tmp_data.tr_con{tr_con}).(tmp_data.position{pos}));
                                            if isequal(tmp_data.tr_pos_con,{'single_targ'}) || isequal(tmp_data.tr_pos_con,{'single_distr'})
                                                c_l_idx = cellfun(@(x) ~isempty(strfind(x,[tmp_data.position{pos},'_choice'])),c_l_fnames,'UniformOutput',1); % use only those fields that contain 'choice' or 'latency', respectively
                                                c_l_fnames = c_l_fnames(c_l_idx); % get only fields containing 'contra_', 'ipsi_' or 'fixation_' -'choice'
                                                if pos == 2
                                                    comparisons = {[2 4]} ;
                                                else
                                                    comparisons = {[1 3]} ;
                                                end
                                            else
                                                c_l_idx = cellfun(@(x) ~isempty(strfind(x,[choices{ch},'_choice'])),c_l_fnames,'UniformOutput',1); % use only those fields that contain 'choice' or 'latency', respectively
                                                c_l_fnames = c_l_fnames(c_l_idx); % get only fields containing 'contra_', 'ipsi_' or 'fixation_' -'choice'
                                            end
                                            % only get tests of choices towards the stimulus; don't include tests on fixation choices
                                            if Settings.FDR
                                                sig_comp_idx = logical([stat_struct_ttest.(tmp_data.tr_con{tr_con}).(tmp_data.position{pos}).(c_l_fnames{1}){col,1}{:,9}]);
                                            else
                                                sig_comp_idx = logical([stat_struct_ttest.(tmp_data.tr_con{tr_con}).(tmp_data.position{pos}).(c_l_fnames{1}){col,1}{:,5}]);
                                            end
                                            %% Plot the stars
                                            if sum(sig_comp_idx > 0) % if at least one t-test was significant
                                                groups = comparisons(sig_comp_idx); % find all significant ttests by column 5
                                                sig_comp = [comparisons(sig_comp_idx)];
                                                for comp = 1:numel(sig_comp)
                                                    groups{comp}(1) = all_x_values(sig_comp{comp}(1),1); % pos KK changed% combinations of xpos values that define those pairs of bars
                                                    groups{comp}(2) = all_x_values(sig_comp{comp}(2),1); % pos KK changed %that are significantly different
                                                end
                                                if Settings.FDR
                                                    all_p_values = [stat_struct_ttest.(tmp_data.tr_con{tr_con}).(tmp_data.position{pos}).(c_l_fnames{1}){col,1}{:,10}];
                                                else
                                                    all_p_values = [stat_struct_ttest.(tmp_data.tr_con{tr_con}).(tmp_data.position{pos}).(c_l_fnames{1}){col,1}{:,6}];
                                                end
                                                stats = all_p_values(sig_comp_idx);
                                                starlevels = [0.05, 1E-2, 1E-3];
                                                linewidth = 1.5;
                                                fontsize = 14;
                                                
                                                sigstar_UZ(groups,stats,0,starlevels,linewidth,fontsize); % call sigstar
                                            end
                                        end
                                    end
                                end
                            end
                            
                            %%% SAVE PLOT %%%
                            if Settings.save_figures
                                save_dir = get_save_dir(Settings.folders, Settings);
                                save_figure_as(fig, Settings, file_name, {'.ai','.png'}, save_dir) % save
                                if Settings.close_after_saving
                                    close(gcf)
                                end
                            end
                        else
                            close(gcf)
                        end
                    end
                end
            end
        end
    end
end
end

function hist_RTs(Settings,plot_struct)
bw_colormap = [100,100,100;96,196,222;195,210,0;239,125,21]./255;

if  strcmp(Settings.Experiment , 'Inactivation') 
    legend_entries = {'Pre','Post'};
    legend_names = {'Pre','Post'};

elseif  strcmp(Settings.Experiment , 'Microstimulation') 
legend_entries = {'control','-250 ms','-100 ms','+50 ms'};
legend_names = {'control','-250 ms','-100 ms','+50 ms'}; %kk

end
%legend_names = {'control','-80 ms','at GO','+80 ms'};
for col = 1:size(plot_struct,1)
    for tr_pos_con = 1:size(plot_struct,2)
        if ~isempty(plot_struct{col,tr_pos_con})
            per_or_across = size(plot_struct{col,tr_pos_con},1); % how many rows? 1st row = per session(one column per session); 2nd row = across sessions (one column only)
            for p_a = 1:per_or_across
                num_sessions = size([plot_struct{col,tr_pos_con}{p_a,:}],2); % how many columns = sessions
                for fol = 1:num_sessions;
                    if ismember('latency',fieldnames(plot_struct{col,tr_pos_con}{p_a,fol}))
                        tmp_data = plot_struct{col,tr_pos_con}{p_a,fol}.latency;
                        num_choices = size(tmp_data.latencies,1);
                        fig = figure;
                        figure_size = get(fig,'Position');
                        figure_size(3) = figure_size(3)*num_choices;
                        set(fig, 'Position', figure_size)
                        for ch = 1:num_choices
                            han(ch) = subplot(1,num_choices,ch); hold on;
                            for stim = 1:size(tmp_data.latencies,2)
                                if ~isempty(tmp_data.latencies{ch,stim}) % if there is data for the histogram
                                    
                                    data2hist = double(tmp_data.latencies{ch,stim});
                                    % title
                                    title_pc_l = 'Latency_Hist';
                                    title_pses = get_title_session_part(tmp_data.session_s, Settings);
                                    title_ptrcon = tmp_data.tr_pos_con;
                                    title_pcol = ['G-R-ratio ',num2str(tmp_data.G_R_ratio)];
                                    title_pch = [tmp_data.groupnames{ch,1},'_saccades'];
                                    hist_title = [title_pc_l,' - ',title_ptrcon,' - ',title_pcol,' - ',title_pch]; %title_pses,' - ',
                                    nbins = round(numel(data2hist) / 5);
                                    if nbins < 10
                                        nbins = 10;
                                    elseif nbins > 50 % never have more than 50 bins
                                        nbins = 50;
                                    end
                                    %%% PLOT THE HISTOGRAM %%%
                                    data = histc(data2hist,(0.12:0.02:0.5));
                                    data = hist2per(data);
                                    h= plot((0.12:0.02:0.5),data,'Color',bw_colormap(stim,:),'LineWidth',2);
                                    % histogram(data2hist,nbins)
                                    t = title(hist_title);
                                    set(t,'interpreter','none')
                                    set(gcf,'color','w'); % sets background color to white
                                    % set(handles.legend,'interpreter','none')
                                    % add text
                                    xlim([0.12 0.5]);
                                    ylim([0 50]);
                                    ylimits = get(gca,'YLim');
                                    yrange = ylimits(2) - ylimits(1);
                                    xpos = 0.4;
                                    hh = legend(legend_names);
                                    set(hh,'interpreter','none');
                                    text(xpos,(yrange-1/4*yrange)- yrange/6*stim,tmp_data.txt{ch,stim},'HorizontalAlignment','left','VerticalAlignment','bottom','Color',h.Color);
                                end
                            end
                            
                            xlabel('Reaction Time [ms]')
                            ylabel('No. saccades [%]');
                            % subplot_position = get(han(ch), 'position');
                            % le = 0.13*ch +0.775*(ch-1); bo = 0.11; wi = 0.775; he = 0.815;
                            % set(han(ch), 'position', [le, bo, wi, he]);
                            
                        end
                        % figure('units','normalized','position',[.1 .1 .4 .4])
                        
                        %%% SAVE PLOT %%%
                        if Settings.save_figures
                            session = strsplit(tmp_data.session_s ,' - ');
                            file_name = [title_pc_l,' - ',title_ptrcon,' - ',title_pcol]; %,title_pses,' - '
                            save_dir = get_save_dir(session, Settings);
                            save_figure_as(fig, Settings, file_name, {'.ai','.png'}, save_dir) % save
                            if Settings.close_after_saving
                                close(gcf)
                            end
                        end
                    end
                end
            end
        end
    end
end

end

function save_dir = get_save_dir(session, Settings)
if iscell(session)
    for ses = 1:numel(session)
        for dtset = 1:numel(fieldnames(Settings.Datasets))
            Dtsets = fieldnames(Settings.Datasets);
            if ismember(session{ses},Settings.Datasets.(Dtsets{dtset}))
                save_dir{dtset,ses} = Settings.save_dir.(Dtsets{dtset}); % get save directory respective to the dataset
            end
        end
    end
elseif ischar(session)
    for dtset = 1:numel(fieldnames(Settings.Datasets))
        Dtsets = fieldnames(Settings.Datasets);
        if ismember(session,Settings.Datasets.(Dtsets{dtset}))
            save_dir{dtset,1} = Settings.save_dir.(Dtsets{dtset}); % get save directory respective to the dataset
        end
    end
end
if exist('save_dir','var')
    save_dir =  save_dir(~cellfun('isempty',save_dir)); % delete empty cells so that 'unique' can be applied
    save_dir = char(unique(save_dir));
else % if no save directory is specified, i.e. if this session is not included in any of the Datasets
    save_dir = Settings.save_dir.No_save_dir_specified;
end
if ~isequal(size(save_dir,1),1) % if data from several datasets is included
    save_dir = Settings.save_dir.No_save_dir_specified;
end
end

function file = df_to_txt(file)
column = strcmp(file(1,:),'df');
df = file(2:end,column);
if iscell(df)
    df_mat = double(cell2mat(df));
    for row = 1:size(df_mat,1)
        df_str{row,1} = num2str(df_mat(row,:));
        df_cell{row,1} = strrep(df_str{row,1},'  ',',');
    end
    file(2:end,column) = df_cell;
end
end

function save_figure_as (fig, Settings, file_name, file_format, save_dir)
% for dtset = 1:numel(fieldnames(Settings.Datasets))
%     Dtsets = fieldnames(Settings.Datasets);
%     if ismember(session,Settings.Datasets.(Dtsets{dtset}))
%         save_dir = Settings.save_dir.(Dtsets{dtset}); % get save directory respective to the dataset
%     end
% end
% if nargin < 5 || ~exist('save_dir','var') % if no save directory is specified, i.e. if this session is not included in any of the Datasets
%     save_dir = Settings.save_dir.No_save_dir_specified; %
% end
for frmt = 1:numel(file_format)
    compl_save_dir = [save_dir filesep 'Figures_' file_format{frmt}(2:end)];
    if ~isequal(exist(compl_save_dir, 'dir'),7)
        mkdir(compl_save_dir)
    end
    switch file_format{frmt}
        case '.ai'
            set(fig,'Renderer','Painters');
            set(fig,'PaperPositionMode','auto')
            compl_filename = char(fullfile(compl_save_dir,strcat(file_name,file_format{frmt})));
            print(fig,'-depsc',compl_filename);
        case '.png'
            compl_filename = char(fullfile(compl_save_dir,strcat(file_name,file_format{frmt})));
            
            % Saves image larger but turned by 90 degrees
            %             set(fig,'PaperOrientation','landscape', 'PaperUnits','normalized', 'PaperPosition', [0 0 1 1]);
            %             print(fig,'-dpng','-r0',compl_filename);
            
            % Turns image to correct orientation but size is smaller
            %             set(fig_handle,'PaperSize',fliplr(get(fig_handle,'PaperSize')),'PaperUnits','inches'); % this flips the image to the ipsi, because by default it will be printed tilted to the contra
            %             compl_filename = char(fullfile(compl_save_dir,strcat(file_name,file_format{frmt})));
            %             print(fig_handle,'-dpng','-r0',compl_filename);
            
            % Save in the size displayed on the screen
            set(fig,'PaperPositionMode','auto')
            print('-dpng','-r0',compl_filename)
            
            % Specify size of the saved file
            %             fig.PaperUnits = 'inches';
            %             fig.PaperPosition = [0 0 6 3];
            %             print('-dpng','-r0',compl_filename)
        case '.pdf'
            set(fig,'PaperOrientation','landscape', 'PaperUnits','normalized', 'PaperPosition', [0 0 1 1]);
            compl_filename = char(fullfile(compl_save_dir,strcat(file_name,file_format{frmt})));
            print(fig, '-dpdf', '-r300',compl_filename)
    end
end
disp(['Plot saved: ',compl_filename])
end

function title_pses = get_title_session_part(session_s, Settings)
    title_pses = session_s;
end

