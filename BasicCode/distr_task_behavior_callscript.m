%% Call analyze_hitrate
% This script calls the function 'distr_task_behavior_analysis'.
% For details on the function type 'help distr_task_behavior_analysis'.
function distr_task_behavior_callscript(Monkey)
clear Settings
%clear all
Settings.Monkey = Monkey;
Computer = 'DPZ_KristinPC'; % 'DPZ','mine'
switch Computer
    case 'mine'
        Settings.directory                              = 'D:\MASTERS\Server_02\Data\Curius\setup1_eye';
        path = 'D:\MASTERS\';
    case 'DPZ'
        Settings.directory                              = 'Z:\dag02\Data\Curius\setup1_eye';
        path = 'Z:\dag02\Projects\';
    case 'DPZ_KristinPC'
       % path = ['Y:\Projects\Pulv_distractor_spatial_choice\Data\',Settings.Monkey, filesep, 'behavior',filesep];
        path = ['Y:\Projects\Pulv_Inac_ECG_respiration\Data\',Settings.Monkey, filesep, 'behavior', filesep ];
end

%% Data settings

switch Settings.Monkey
    
   %% probblem: 20211207
    case 'Bacchus'
        Settings.hemisphere_of_stimulation  = 'left';          % hemisphere of the brain where the stimulation happens. Input: 'right' or 'left'
        Settings.directory                  = 'Y:\Data\Bacchus';
        Settings.Experiment = 'Electrophysiology';
        Settings.Datasets.EpyhsExpDPul_20210706_until_20210930       = {'20210706','20210714','20210715','20210826',...
            '20210722','20210723','20210729','20210730', '20210803','20210805','20210806' ,'20210826', '20210827',...
            '20210829','20210905' ,'20210930','20210906'}; % Add sessions to a dataset to make sure it is saved in its respective folder
         Settings.Datasets.EpyhsExpDPul_202110_until_202112       = {'20211001','20211005','20211007','20211012',...
            '20211013','20211014','20211019','20211025', '20211027','20211028','20211102' ,'20211103','20211108', '20211111','20211116','20211117' ,...
            '20211214','20211222' };
        
%         Settings.Datasets.EpyhsExpDPul_202201_until_202203       = {'20220105','20220106','20220125','20220126',...
%             '20220203','20220211','20220221','20220222', '20220224','20220225','20220309' ,'20220310'};
%        % Settings.Datasets.Training_20210616_20210617        = {'20210616', '20210617' }; % Add sessions to a dataset to make sure it is saved in its respective folder
        Settings.folders                               =  Settings.Datasets.EpyhsExpDPul_202110_until_202112  ; %Settings.Datasets.TRAIN_pre_stim; Settings.Datasets.PSF;
        
        Settings.exclude_runs_of_sessions              =  [{'none'}] ;%[{'20190725','2','8','10','4','5','11','13','12','14'}];  %[{'none'}];
        Settings.load_data                             = ''; %[path, 'PSF\PSF_20171018-1340.mat']; %[path, 'PSF\PSF_20171018-1340.mat']; %''; %
        
        Settings.save_dataset                          = '\Behavior_Select_RT'; % ; %'Behavior_Select_RT' Behavior_Baseline 'Behavior_BaselineVPul';'Behavior_InactivationVPul';
        Settings.save_figures                                       = 1;
        Settings.close_after_saving                                 = 1;
        Settings.save_dir.EpyhsExpDPul_202110_until_202112            = [path,'EpyhsExpDPul_202110_until_202112\'];
        
        Settings.pathExcel          = 'Y:\Logs\Inactivation\Bacchus';
        Settings.Excelfilenames     = 'Bacchus_inactivation_log.xlsx';
        table                       = readtable([Settings.pathExcel,filesep, 'Bacchus_inactivation_log.xlsx' ]);
        
        
    case 'Curius'
        Settings.Experiment = 'Microstimulation'; %Electrophysiology
        Settings.hemisphere_of_stimulation              = 'right';          % hemisphere of the brain where the stimulation happens. Input: 'right' or 'left'
        Settings.directory                              = 'Y:\Data\Curius';
        Settings.exclude_runs_of_sessions               =  [{'none'}] ;%[{'20190725','2','8','10','4','5','11','13','12','14'}];  %[{'none'}];
        Settings.load_data                              = ''; %[path, 'PSF\PSF_20171018-1340.mat']; %[path, 'PSF\PSF_20171018-1340.mat']; %''; %
        
        Settings.save_dataset                           = '\Behavior_Select_RT'; % ; %'Behavior_Inactivation' Behavior_Baseline 'Behavior_BaselineVPul';'Behavior_InactivationVPul';
        Settings.save_figures                                       = 1;
        Settings.close_after_saving                                 = 1;
        Settings.pathExcel          = 'Y:\Logs\Inactivation\Curius';
        Settings.Excelfilenames     = 'Curius_Inactivation_log_since201905.xlsx';
        table                       = readtable([Settings.pathExcel,filesep, 'Curius_Inactivation_log_since201905.xlsx' ]);
        
        if  strcmp(Settings.Experiment , 'Inactivation')
             Settings.directory                              = 'Y:\Data\Curius';            
            Settings.Datasets.Inactivation_11Sess_Behav_ECG     = {'20190729','20190801','20190809','20190807','20190814','20190820','20190826','20190905', '20190913'}; % ,Add sessions to a dataset to make sure it is saved in its respective folder
           Settings.Datasets.Baseline_11Sess_Behav_ECG         = {'20190802','20190804','20190806','20190808','20190811','20190813', '20190815','20190903','20190910'}; %  , '20190912',Add sessions to a dataset to make sure it is saved in its respective folder

            
            % Settings.Datasets.Inactivation_20200422_20200424_20200429_20200505_20200513_20200518_20200625 = {'20200422','20200424','20200429','20200505', '20200513', '20200518', '20200625'}; % Add sessions to a dataset to make sure it is saved in its respective folder
            %Settings.Datasets.Baseline_20200427_20200504_20200507_20200508_20200512_20200623_20200626 = {'20200427', '20200504', '20200507', '20200508', '20200512', '20200623' , '20200626'}; % Add sessions to a dataset to make sure it is saved in its respective folder
            %Settings.folders                               =  Settings.Datasets.Inactivation_20200424_20200429_20200505_20200513_20200518  ; %Settings.Datasets.TRAIN_pre_stim; Settings.Datasets.PSF;
            %Settings.folders                               =  Settings.Datasets.Baseline_20200427_20200504_20200507_20200508_20200512_20200623_20200626  ; %Settings.Datasets.TRAIN_pre_stim; Settings.Datasets.PSF;
           
            Settings.folders                               =  Settings.Datasets.Inactivation_11Sess_Behav_ECG; 
           Settings.folders                               =  Settings.Datasets.Baseline_11Sess_Behav_ECG; 
            
            %  Settings.folders                               =  Settings.Datasets.Inactivation_20190729_20190801_20190809_20190814_20190820_20190905_20190913  ; %Settings.Datasets.TRAIN_pre_stim; Settings.Datasets.PSF;
            %   Settings.folders                              =  Settings.Datasets.Baseline_20190912_20190910_20190903_20190815_20190808_20190806_20190802  ; %Settings.Datasets.TRAIN_pre_stim; Settings.Datasets.PSF;
            % Settings.save_dir.Inactivation_20200422_20200424_20200429_20200505_20200513_20200518_20200625           = [path,'Inactivation_20200422_20200424_20200429_20200505_20200513_20200518_20200625\'];
            
            % Settings.save_dir.Inactivation_20200424_20200429_20200505_20200513_20200518            = [path,'Inactivation_20200424_20200429_20200505_20200513_20200518\'];
            %Settings.save_dir.Baseline_20200427_20200504_20200507_20200508_20200512_20200623_20200626            = [path,'Baseline_20200427_20200504_20200507_20200508_20200512_20200623_20200626\'];
            
            % Settings.save_dir.Baseline_20190912_20190910_20190903_20190815_20190808_20190806_20190802            = [path,'Baseline_20190912_20190910_20190903_20190815_20190808_20190806_20190802\'];
            % Settings.save_dir.Inactivation_20190729_20190801_20190809_20190814_20190820_20190905_20190913            = [path,'Inactivation_20190729_20190801_20190809_20190814_20190820_20190905_20190913\'];
            Settings.save_dir.Inactivation_11Sess_Behav_ECG            = [path,'Inactivation_11Sess_Behav_ECG\'];
           Settings.save_dir.Baseline_11Sess_Behav_ECG            = [path,'Baseline_11Sess_Behav_ECG\'];


        elseif  strcmp(Settings.Experiment , 'Microstimulation')
            Settings.directory                              = 'Y:\Data\Curius\setup1_eye';
            Settings.save_dataset                           = '\Behavior_Select_RT'; % ; %'Behavior_Inactivation' Behavior_Baseline 'Behavior_BaselineVPul';'Behavior_InactivationVPul';

            path = ['Y:\Projects\Pulv_distractor_spatial_choice\', Settings.Experiment, filesep, 'Data\',Settings.Monkey, filesep];
            Settings.Datasets.STIM_early_20161020_until_20161209 = {'20161020', '20161021', '20161026', '20161027', '20161028', '20161102' , '20161103', '20161104',...
                 '20161108', '20161109', '20161110', '20161111', '20161206', '20161207', '20161208', '20161209'}; % Add sessions to a dataset to make sure it is saved in its respective folder
            Settings.folders                               =  Settings.Datasets.STIM_early_20161020_until_20161209  ; %Settings.Datasets.TRAIN_pre_stim; Settings.Datasets.PSF;
            Settings.save_dir.STIM_early_20161020_until_20161209            = [path,'STIM_early_20161020_until_20161209\'];
        elseif strcmp(Settings.Experiment , 'Electrophysiology')
            
        Settings.directory                             = 'Y:\Data\Curius';
        Settings.Datasets.Exp20210311_20210318            = {'20210311' '20210318'}; % Add sessions to a dataset to make sure it is saved in its respective folder
        Settings.folders                               =  Settings.Datasets.Exp20210311_20210318  ; %Settings.Datasets.TRAIN_pre_stim; Settings.Datasets.PSF;
        Settings.save_dir.Exp20210311_20210318            = [path,'Exp20210311_20210318\'];
        Settings.save_dataset                          = 'Training_Select_RT'; % ; %'Behavior_Select_RT' Behavior_Baseline 'Behavior_BaselineVPul';'Behavior_InactivationVPul';

            
        end
        
            %%For Curius defined by UWE
            case 'Cornelius'
                Settings.Experiment = 'Microstimulation';
                Settings.directory                              = 'Y:\Data\Cornelius';
                Settings.hemisphere_of_stimulation              = 'right';          % hemisphere of the brain where the stimulation happens. Input: 'right' or 'left'
                Settings.exclude_runs_of_sessions               =  [{'none'}];
                Settings.load_data                              = ''; %[path, 'PSF\PSF_20171018-1340.mat']; %[path, 'PSF\PSF_20171018-1340.mat']; %''; %
                Settings.save_dataset                           = '\Behavior_Select_RT'; % ; %'Behavior_Inactivation' Behavior_Baseline 'Behavior_BaselineVPul';'Behavior_InactivationVPul';
                Settings.save_figures                                       = 0;
                Settings.close_after_saving                                 = 0;
                Settings.pathExcel = 'Y:\Logs\Inactivation\Cornelius';
                Settings.Excelfilenames = 'Cornelius_Inactivation_log_since201901.xlsx';
                table = readtable([Settings.pathExcel,filesep, 'Cornelius_Inactivation_log_since201901.xlsx' ]);
                
                
                     if  strcmp(Settings.Experiment , 'Inactivation')
                %Settings.Datasets.Inactivation_20190124_20190129_20190201_20190207_20190214_20190228_20190314         = {'20190124','20190129','20190201','20190207','20190214','20190228', '20190314'}; % ,Add sessions to a dataset to make sure it is saved in its respective folder
                Settings.Datasets.Inactivation_11Sess_Behav_ECG         = {'20190124','20190129','20190201','20190207','20190214','20190228', '20190314','20190828','20190904','20190910', '20191011'}; % ,Add sessions to a dataset to make sure it is saved in its respective folder
                Settings.Datasets.Baseline_11Sess_Behav_ECG         = {'20190131','20190213','20190216', '20190227', '20190304', '20190313', '20190403','20190913', '20191007', '20191010', '20191014'}; % ,Add sessions to a dataset to make sure it is saved in its respective folder

                % Settings.Datasets.Baseline_20190131_20190213_20190216_20190227_20190304_20190313                      = {'20190131','20190213','20190216', '20190227', '20190304', '20190313'}; % Add sessions to a dataset to make sure it is saved in its respective folder
                % Settings.Datasets.Baseline_20190121_20190131_20190213_20190216_20190227_20190304_20190313             = {'20190121','20190131','20190213','20190216', '20190227', '20190304', '20190313'}; % Add sessions to a dataset to make sure it is saved in its respective folder
                % Settings.Datasets.InactivationVPul_20190404_20190408_20190430_20190509     = {'20190404','20190408','20190430','20190509'};
                % Settings.Datasets.BaselineVPul_20190403_20190424_20190429_20190508     = {'20190403','20190424','20190429', '20190508'};
               % Settings.folders                                = Settings.Datasets.Inactivation_20190124_20190129_20190201_20190207_20190214_20190228_20190314 ;
                Settings.folders                                = Settings.Datasets.Inactivation_11Sess_Behav_ECG ;
                Settings.folders                                = Settings.Datasets.Baseline_11Sess_Behav_ECG ;

                % Settings.folders                                =  Settings.Datasets.Baseline_20190131_20190213_20190216_20190227_20190304_20190313 ; %Settings.Datasets.TRAIN_pre_stim; Settings.Datasets.PSF;
                % Settings.folders                                =  Settings.Datasets.Baseline_20190121_20190131_20190213_20190216_20190227_20190304_20190313 ; %Settings.Datasets.TRAIN_pre_stim; Settings.Datasets.PSF;
                % Settings.folders                                = Settings.Datasets.InactivationVPul_20190404_20190408_20190430_20190509;
                %Settings.folders                                = Settings.Datasets.BaselineVPul_20190403_20190424_20190429_20190508;
                %Settings.folders                                = Settings.Datasets.Baseline ; %Settings.Datasets.TRAIN_pre_stim; Settings.Datasets.PSF;
                %% Why did we exclude the runs?
                %Settings.exclude_runs_of_sessions               =  [{'20190129','10','11','14'}]; % [{'20171103','11'}]; %[{'20171026','9'}];...
                %Settings.exclude_runs_of_sessions               =  [{'20190124','1','5'}]; % [{'20171103','11'}]; %[{'20171026','9'}];...
                   %Settings.save_dataset                              = 'Behavior_Baseline_20190121_20190131_20190206'; % ; %'PSF'  'Training_DifficultyDistractors';'';
                % Settings.save_dir.Inactivation_20190215           = strcat(fullfile( path,Settings.Monkey,filesep,Settings.Experiment ,filesep) ,Settings.Experiment, '_',  Settings.Datasets.Inactivation_20190215);
                %Settings.save_dir.Baseline_20190912           = [path,'Baseline_20190912\'];
                %Settings.save_dir.Inactivation_20190124_20190129_20190201_20190207_20190214_20190228_20190314                 = [path,'Inactivation_20190124_20190129_20190201_20190207_20190214_20190228_20190314\'];
                Settings.save_dir.Inactivation_11Sess_Behav_ECG                 = [path,'Inactivation_11Sess_Behav_ECG\'];
                Settings.save_dir.Baseline_11Sess_Behav_ECG                 = [path,'Baseline_11Sess_Behav_ECG\'];

                % Settings.save_dir.Baseline_20190131_20190213_20190216_20190227_20190304_20190313               = [path,'Baseline_20190131_20190213_20190216_20190227_20190304_20190313\'];
                % Settings.save_dir.Baseline_20190121_20190131_20190213_20190216_20190227_20190304_20190313               = [path,'Baseline_20190121_20190131_20190213_20190216_20190227_20190304_20190313\'];
                % Settings.save_dir.InactivationVPul_20190404_20190408_20190430_20190509       = [path,'InactivationVPul_20190404_20190408_20190430_20190509\'];
                % Settings.save_dir.BaselineVPul_20190403_20190424_20190429_20190508           = [path,'BaselineVPul_20190403_20190424_20190429_20190508\'];
                     elseif  strcmp(Settings.Experiment , 'Microstimulation')
                         
            path = ['Y:\Projects\Pulv_distractor_spatial_choice\', Settings.Experiment, filesep, 'Data\',Settings.Monkey, filesep];
           % Settings.Datasets.PSF_4session_DifficultyDiffToCur = { '20171012', '20171013','20171017', '20171018'}; % '20171010', '20171011'
           % Settings.Datasets.PSF_5session_DifficultyDiffToCur = { '20171010', '20171011'}; %
            % Settings.Datasets.PSF_5session_DifficultyDiffToCur = { '20171010', '20171019'}; % 
          
            
            Settings.Datasets.STIM_late_20180109_until_20180123 = { '20180109','20180111', '20180113', '20180114','20180115', '20180119', '20180121', '20180122', '20180123'}; % '20171126', '20171127', '20171201', '20171202', '20171203', '20171204' , '20171205', '20171206',...
             %    '20171210', '20171211', '20171212', '20171213',
            Settings.folders                               =  Settings.Datasets.STIM_late_20180109_until_20180123  ; %Settings.Datasets.TRAIN_pre_stim; Settings.Datasets.PSF;
            Settings.save_dir.STIM_late_20180109_until_20180123            = [path,'STIM_late_20180109_until_20180123\'];


                     end
                
           
        end %monkey specific
        
        
        %% load excel-file
        
        
        Settings.PreInjection_Runs = [];
        
        for indSession = 1: length(Settings.folders)
            if sum(table.date == str2num([Settings.folders{indSession}])) == 0
                Settings.SeesionInfo{indSession} = 'Training';
                Run =  NaN ;
                Settings.PreInjection_Runs{indSession} = NaN;
            else
                Settings.PreInjection{indSession} =   table.injection(table.date == str2num([Settings.folders{indSession}]))';
                Settings.SessionExp{indSession} =   table.experiment(table.date == str2num([Settings.folders{indSession}]))';
                Settings.SessionInfo = Settings.SessionExp{1,1}{1,1};
                %Which run is it for Pre?
                Run = table.run(table.date == str2num([Settings.folders{indSession}]))';
                Settings.PreInjection_Runs{indSession} =  Run(strcmp(Settings.PreInjection{indSession}, 'Pre'));
            end
        end
        
        
        % % find all dates which are related to a specific experiment
        % ExperimentNames = unique(table.experiment);
        %
        % s = [ExperimentNames(2)];
        % Settings.Datasets.(char(ExperimentNames(2))) = num2cell(unique(table.date(strcmp(table.experiment,'Injection'))))';
        % % readout the runs for pre vs post-injection for each Date
        % Injection.PrePost = num2cell(table.injection(strcmp(table.experiment,'Injection')))';
        
        %How to add the variables to the SETTING structure
        
        %% Task settings
        Settings.target_color_dim           = [128,0,0];
        Settings.fixation_color_dim         = [60,60,60];
        
        %% Analysis settings
        Settings.task_type                  = 2; % 1 - calibration;
        % 2 - direct;
        % 3 - memory;
        
        Settings.necessary_trial_conditions = {}; %'single_targ', exclude runs  'targ_targ_2HF'that do not contain all necessary trial conditions; % {'targ_targ','distr_distr','targ_distr','single_targ','single_distr'}{'distr_distr','targ_distr'}
        % PSF
        %Settings.analyze_trial_conditions   = {{'targ_distr'}, {'targ_targ' 'distr_distr'}}; %{{'targ_targ' 'targ_distr'}, {'targ_targ' 'distr_distr'}};
        % Baseline & Stimulation Sessions
        %% Why does the order of these combinations matter?
        %Settings.analyze_trial_conditions   = {{'single_targ' 'single_distr'}, {'single_targ' 'targ_distr'}, {'targ_targ' 'targ_distr'}, {'targ_targ' 'distr_distr'}}; %{{'targ_targ' 'targ_distr'}, {'targ_targ' 'distr_distr'}};
        % trial conditions that shall be analyzed; senseful combinations: {{'single_targ' 'single_distr'}, {'single_targ' 'targ_distr'}, {'targ_targ' 'targ_distr'}, {'targ_targ' 'distr_distr'}}
        Settings.analyze_trial_conditions   = {{'single_targ' 'single_distr'}, {'targ_targ_2HF' 'distr_distr_2HF'},{'targ_distr_1HF' 'targ_distr_2HF'}, {'targ_targ_1HF' 'distr_distr_1HF'}}; %,{{'targ_targ' 'targ_distr'}, {'targ_targ' 'distr_distr'}};
        %{'targ_targ_2HF' 'distr_distr_2HF'}
        Settings.analyze_distr_colors       = ['any']; % e.g. 'any'        % 'any'
        % OR Not functioning: [60 60 0]
        
        %% Plot settings
        Settings.min_num_trials_per_run     = 30; %100
        
        Settings.analyze_hitrate            = 1; %?
        Settings.hitrate_mode               = 'correct target selected';   % 'success' - only successful trials are counted as hits
        % 'correct target selected' - all trials where the correct target was selected are counted as hits (incl. target hold aborts)
        Settings.analyze_latency            = 1;
        Settings.analyze_choice             = 1;
        Settings.analyze_error_rates        = 0;                           % analyzes error rates of fixation break, target acquisition and target hold aborts
        Settings.analyze_eye_traces         = 0;                           % analyzes eye traces of correct trials
        %Settings.analyze_by_position = {'left','right','combined'};       % {'contra','ipsi','all'} - plot a separate line for each specified entry of position of target (in target-distr and single target trials)
        Settings.analyze_by_position = {'contra','ipsi','combined'};                                                                    % or distractor (in single distr trials) and, or left and right targets combined
        %Settings.analyze_by_position = {'contra','contra_up','contra_down','ipsi','ipsi_up','ipsi_down','combined'};   % {[]} = {'all'} -  Don't analyze position. Overall analysis of above specified trials
        Settings.analyze_across_sessions    = 1;                           % automatically set on if plot_across_sessions is on; analysis per session is always done
        
        
        %%
        Settings.plot_per_session                           = 0;
        Settings.plot_across_sessions                       = 0;                % if only 1 session is analyzed this is set false by default
        Settings.plot_error_traces                          = 0;                % plots error traces of tacq errors in Ipsi-targ - contra easy distractor trials; only across sessions
        Settings.line_plot                                  = 1;
        Settings.fit_psychometric_curve                     = 1;
        Settings.calc_goodness_of_fit                       = 1;
        Settings.bar_plot                                   = 0; % for latency     -> x-axis: same plot fixation, contra, ipsi .. different task condition
        Settings.bar_plot_II                                = 0; % only for choice -> seprate plots for fixation, contra, ipsi selection
        Settings.bar_plot_choice_mean_or_true_proportion    = 'mean_and_SEM';   % 'mean_and_SEM' - plots the mean and SEM of choice/hitrate proportions across runs (sessions) in a plot per session (across sessions)
        % 'true_proportion'  - plots the true proportion per session (per all sessions), i.e. calculates only one proportion across all trials of a session (of all sessions)
        Settings.bar_plot_latency_avg_mean_or_true_mean     =  'mean_of_means' ;  % 'mean_of_means' - plots the mean of the mean latencies per run (session) in a plot per session (across sessions); SEM is the SEM across these means
        % 'true_mean_of_concatenated_trials' - plots the true mean of all trials per session (per all sessions), i.e. calculates only one mean of all trials of a session (of all sessions)
        Settings.hist_RTs                                   = 0;
        
        
        
        if  strcmp(Settings.Experiment , 'Inactivation') &&  strcmp(Settings.hemisphere_of_stimulation, 'left')
            Settings.GraphProperties.contra.title_ch    = 'Contra_choices';
            Settings.GraphProperties.ipsi.title_ch     = 'Ipsi_choices';
            Settings.GraphProperties.Fixation.title_ch = 'Fixation_choices';
            Settings.GraphProperties.contra.bw_ylabel   = 'Contralesional Selection';
            Settings.GraphProperties.ipsi.bw_ylabel    ='Ipsilesional  Selection';
            
        elseif  strcmp(Settings.Experiment , 'Inactivation') &&  strcmp(Settings.hemisphere_of_stimulation, 'right')
            Settings.GraphProperties.ipsi.title_ch     = 'Ipsi_choices';
            Settings.GraphProperties.contra.title_ch      = 'Contra_choices';
            Settings.GraphProperties.Fixation.title_ch  = 'Fixation_choices';
            Settings.GraphProperties.ipsi.bw_ylabel   = 'Ipsilesional  Selection';
            Settings.GraphProperties.contra.bw_ylabel    = 'Contralesional Selection';
            
        elseif  (strcmp(Settings.Experiment , 'Microstimulation')  || strcmp(Settings.Experiment ,'Electrophysiology'))&&   strcmp(Settings.hemisphere_of_stimulation, 'right')
            Settings.GraphProperties.ipsi.title_ch    = 'Ipsi_choices';
            Settings.GraphProperties.contra.title_ch     = 'Contra_choices';
            Settings.GraphProperties.Fixation.title_ch = 'Fixation_choices';
            
        elseif  strcmp(Settings.Experiment , 'Microstimulation') &&  strcmp(Settings.hemisphere_of_stimulation, 'left')
            Settings.GraphProperties.contra.title_ch    = 'Contra_choices';
            Settings.GraphProperties.ipsi.title_ch     = 'Ipsi_choices';
            Settings.GraphProperties.Fixation.title_ch = 'Fixation_choices';
        end
        
        %% Statistics Settings
        % ERROR RATE: difference in fixation breaks between stimulation and no
        % stimulation
        
        
        Settings.statistics_on                                                          = 0;
        Settings.save_stat_output                                                       = 0; % save statistical results as mat and txt files to pre-defined save_dirs, may return an error
        Settings.alpha                                                                  = 0.05; % for all analyses incl. ANOVAs
        Settings.ttests                                                                 = 1; % t-tests for bar plots; results of these t-tests may interfere (be overwritten) by posthoc t-tests, so turn on only one or the other
        Settings.FDR                                                                    = 1; % % compare to control, correct for multiple testing using the False Discovery Rate procedure (Benjamini & Hochberg), for Settings.ttests and Settings.posthoc_ttest
        Settings.choice_rmanova3way_DoubleStimDisplay_Difficulty_Stimulation            = 0; % this design was dismissed and not used in final analysis
        Settings.choice_rmanova2way_DistrDistrTargTarg_Stimulation                      = 0; %1
        Settings.choice_rmanova3way_TargDistr_Difficulty_Stimulation                    = 0; %1
        Settings.choice_rmanova3way_SingleStimDisplay_Side_Stimulation                  = 0; %1
        Settings.error_rate_rmanova2way_DistrDistrTargTarg_Stimulation                  = 0; %1
        Settings.error_rate_rmanova3way_TargDistr_Difficulty_Stimulation                = 0; %1
        Settings.error_rate_rmanova3way_SingleStimDisplay_Side_Stimulation              = 0; %1
        Settings.latency_rmanova2way_SingleStimDisplay_Side_Stimulation                 = 0; %1
        Settings.latency_diff_rmanova2way_SingleStimDisplay_Side_Stimulation            = 0; %1
        Settings.latency_rmanova3way_DoubleStimDisplay_Choice_Stimulation               = 0; %1
        Settings.latency_diff_rmanova3way_DoubleStimDisplay_Choice_Stimulation          = 0; %1 % define display in line 1056 in distr_task_behavior_statistics
        Settings.latency_anova2way_SingleStimDisplay_Side_Stimulation_ses_by_ses        = 0; % calculates the latency_diff ANOVAS session by session instead of across sessions
        Settings.latency_anova3way_DoubleStimDisplay_Choice_Stimulation_ses_by_ses      = 0; % not working; calculates the latency_diff ANOVAS session by session instead of across sessions
        Settings.latency_kruskal_wallis_Session_by_session                              = 0; % not used in final analysis
        Settings.latency_posthoc_Mann_Whitney_U                                         = 0; % not used in final analysis
        Settings.latency_differences_wilcoxon                                           = 0; % not used in final analysis
        Settings.choice_anova2way_Baseline_Control                                      = 0;
        Settings.posthoc_ttest                                                          = 0; % post-hoc t tests for ANOVAs, not plotted in bar plots
        
        
        %% Call analyze_saccades
        
        [data_struct_per_run, data_struct_per_session, data_struct_across_sessions, plot_struct_double, plot_struct_single, errors, stat_results_anova, stat_struct_ttest, stat_struct_ttest_anova, stat_struct_chi2] = ...
            distr_task_behavior_analysis(Settings);
        
        %%%%
        % Settings.bar_plot_latency_avg_mean_or_true_mean     = 'true_mean_of_concatenated_trials';
        % Settings.plot_per_session = true;
        % Settings.plot_across_sessions = false;
        % Settings.bar_plot_mean_or_true_proportion = 'true_proportion';
        % [data_struct_per_run, data_struct_per_session, data_struct_across_sessions, plot_struct_double, plot_struct_single, errors, stat_results] = distr_task_behavior_analysis(Settings);
        %
        % a=2;
        
