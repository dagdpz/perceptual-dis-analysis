keys.project_versions={''};
keys.project_version='dPul_inj_LIP_Bac_8s'; %Is this still used ?
keys.filelist_formatted={};

%% what to plot
keys.plot.single_cells =0;
keys.plot.waveforms=0;
keys.plot.population_PSTH_legends=1;  
%% to check carefully
% keys.position_and_plotting_arrangements             ={'hands'};
keys.position_and_plotting_arrangements             ={'hands_inactivation'};

%% computation settings
keys.cal.datasets                   =[10];
keys.cal.effectors                  =[0];
keys.cal.reach_hand                 =[0];
keys.cal.types                      =[4];
keys.cal.units_from_sorting_table   =1;

%% batching
keys.batching.combine_monkeys           =1;
keys.batching.monkeys                   ={'Bacchus'};
% keys.Tesla.date                      ='[20160217 20180101]';
%keys.Linus.date                       ='[20161103 20180101]';
keys.Bacchus.date                       ='[20201112 20210401]';
keys.batching.targets                  = {'LIP_L','LIP_R'};
keys.contra_ipsi_relative_to           ='perturbation_site';
keys.plot.polars_on_extra_figure        =0;

%% criterions to exclude trials and units

keys.cal.stablity                       =[1];                     % not assigning sorting table information if criterion is not met. Therefore only excludes when taking only units in the table
keys.cal.single_rating                  =[1,2,3];                   % not assigning sorting table information if criterion is not met. Therefore only excludes when taking only units in the tabl
keys.cal.SNR_rating                     =[1,2,3];                 % not assigning sorting table information if criterion is not met. Therefore only excludes when taking only units in the tabl
keys.cal.min_trials_per_condition       =5;                         % minimum trials per conditon (look at ph_arrange_positions to see how conditions are defined)
keys.cal.min_spikes_per_unit            =10;                        % excluding units that have in total less spikes (workaround for sortcode assignment bug) - to be removed
keys.cal.perturbation_groups            ={12,2};       % which perturbation values from excel table will be assigned to control and perturbation for comparisons and population analysis

%% epochs

  keys.EPOCHS_PER_TYPE{4}={...
     'INI',      2,	-0.4,	-0.1,   'INI';...
    'Facq',     3,	-0.4,	-0.1,   'INI';...
     'Fhol',     6,	-0.425,	-0.125,      'INI';...
    'Cue',      6,	0.05,   0.35,   'INI';...
%     'EDel',     8, 	0.3,    0.6,   'INI';...
    'Del',      4, 	-0.5,   -0.2,      'INI';...
     'PreS',     60,	-0.3, 	0,  'INI';...
    'PeriS',	60,	-0.15, 	0.15,   'INI';...
     'PostS',	61,	0.05,   0.35,    'INI';...
%     'PreR',     62,	-0.3, 	-0.05,  'INI';...
%     'PeriR',	62,	-0.05, 	0.15,   'INI';...
%     'PostR',	63,	0.05,   0.2,    'INI';...
%      'Thol',     20,	-0.1,   0,      'INI';...     
    };

keys.WINDOWS_PER_TYPE{4}={...
    'Delay Period', 6,	-0.33,  0.8;... %1.35
%     'Reach',        62,	-0.35,  0.7;... %-0.35
     'Saccade',        60,	-0.8,  0.3;... %-0.35
    };  
    
    keys.ANOVAS_PER_TYPE(4).epoch={'INI' 'Facq';...
    'INI' 'Fhol';...
    'INI' 'Cue';...
    'INI' 'Del';...
    'INI' 'PreS';...
    'INI' 'PeriS';...
    'INI' 'PostS';...
    %     'INI' 'PreR';...
    %     'INI' 'PeriR';...
    'INI' 'Thol'};


% keys.ANOVAS_PER_TYPE(4).spaceLR            ={'INI','Facq','Fhol','Cue','Del','PreR','PeriR','PostR','Thol'}';
% keys.ANOVAS_PER_TYPE(4).positions          ={'Facq','Fhol','Cue','Del','PreR','PeriR','PostR','Thol'}';
% keys.ANOVAS_PER_TYPE(4).hands              ={'INI','Facq','Fhol','Cue','Del','PreR','PeriR','PostR','Thol'}';
% keys.ANOVAS_PER_TYPE(4).SxH                ={'INI','Facq','Fhol','Cue','Del','PreR','PeriR','PostR','Thol'}';
% keys.ANOVAS_PER_TYPE(4).main               ={'INI','Facq','Fhol','Cue','Del','PreR','PeriR','PostR','Thol'}';

keys.ANOVAS_PER_TYPE(4).spaceLR            ={'INI','Facq','Fhol','Cue','Del','PreS','PeriS','PostS','Thol'}';
keys.ANOVAS_PER_TYPE(4).positions          ={'Facq','Fhol','Cue','Del','PreS','PeriS','PostS','Thol'}';
keys.ANOVAS_PER_TYPE(4).hands              ={'INI','Facq','Fhol','Cue','Del','PreS','PeriS','PostS','Thol'}';
keys.ANOVAS_PER_TYPE(4).SxH                ={'INI','Facq','Fhol','Cue','Del','PreS','PeriS','PostS','Thol'}';
keys.ANOVAS_PER_TYPE(4).main               ={'INI','Facq','Fhol','Cue','Del','PreS','PeriS','PostS','Thol'}';


%% tuning table selection
% Most important parameter to understand how cell seleciton works in the population analysis, is the trial criterion per condition:
% units with less than this amount of trials for one of the conditions will be excluded. 
keys.cal.min_trials_per_condition       =5;
% This criterion is applied every time the tuning table is loaded via ph_load_extended_tuning_table and ph_reduce_tuning_table. 
% The tricky part here is to define the conditions in which we ask for a specific number of trials. 
% Each condition is defined by tasktype(contains type,effector,and
% position_arrangement),position (depends on arrangement),hand used and choice/instruced trial (MISSING HERE: PERTURBATION)
keys.tt.tasktypes                   ={'Dsa_han'}; % typically only one tasktype defines
keys.tt.type_effectors                   ={'Dsa'}; % typically only one tasktype defines
keys.tt.hands                       =[];
keys.tt.choices                     =[0]; %IMPORTANT: and also not really perfect, for choice trials trial criterion is applied by hemifield, not by position.
% Each unique combination of the above parameters has to contain at least keys.cal.min_trials_per_condition trials, if not the cell is excluded in ph_reduce_tuning_table
keys.tt.selection                   ={};                         % easy to use if there is a parameter in the tuning table for which you want your cells to have the same value
%                                       'in_NH_TIhol_position_Msac_opt','true'};  % each row in the cell arryáy will be used to exclude cells that don't have the specifie characteristic
% keys.tt.unselect                    ={'SNR_rating',3,'Single_rating',3}; % easy to use if there is a parameter in the tuning table for which your cells shouldn't have a specific value
keys.tt.unselect                    ={};
keys.tt.combine_tuning_properties   ={'hand_tuning','in_IH_Facq_epoch_Ddsa_han','in_CH_Facq_epoch_Ddsa_han'}; % results in enen/ensu/suen/en-/-en/su-/-su/--
% ph_load_extended_tuning_table will create an additional column from combining existing columns: 
% 1st entry is the name of the new column (to refer to it), and the following entries specify the columns which should be combined.
% this additional column can be used for grouping, specific selection and/or unselection ph_load_extended_tuning_table also creates new columns which can be used
% for multicomparison correction (sort of) by excluding cells that did not show certain combination of ANOVA effects (this is ONLY used in cell_counts so far
keys.tt.epoch_criterion             ='none'; % only relevant for cell counts
keys.tt.space_criterion             ='none';
keys.tt.hands_criterion             ='none';
keys.tt.SXH_criterion               ='none';
keys.tt.trial_criterion_in          = 'per_hemifield_and_perturbation';
keys.tt.trial_criterion_ch          = 'per_hemifield_and_perturbation';

%% population PSTH settings
cc=0;
% 1
cc=cc+1;
keys.pop(cc).tt.hands                 	= [0];
keys.pop(cc).tt.perturbations          	= [0 1];
keys.pop(cc).tt.choices                	= 0;
keys.pop(cc).tt.selection             	= {};
% keys.pop(cc).tt.unselect                ={'SNR_rating',3,'Single_rating',3}; 
keys.pop(cc).tt.unselect                ={}; 
keys.pop(cc).tt.tasktypes               = {'Dsa_han'};
keys.pop(cc).group_parameter            = 'ungrouped'; %hand_tuning
keys.pop(cc).conditions_to_plot         = {'Dsac'}; 
keys.pop(cc).epoch_PF                   = 'Cue';               % epoch in which preference defines target location for "pref" plots
keys.pop(cc).epoch_RF                   = 'Cue';               % epoch for which gaussian response fields will be plotted (if plot_RF ~ 0)
keys.pop(cc).epoch_BL                   = 'Fhol';                % Epoch to subtract trial by trial (if FR_subtract_baseline ~ 0)
keys.pop(cc).FR_subtract_baseline       = 0;   
keys.pop(cc).group_excluded             = {};%{'susu','ensu','suen','-su','su-','--'};
keys.pop(cc).epoch_for_normalization    = 'Fhol';               % epoch used for (divisive) normalization
keys.pop(cc).normalization              = 'by_perturbation';        % separate (divisive) normalization factor for trials grouped by effector; other options:
            