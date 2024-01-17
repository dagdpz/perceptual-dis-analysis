function SDT_doubleSameStimuli_2H(monkey, folder,dataset,folder_baseline,dataset_baseline,experiment,path_SaveFig,path_SaveData,path_getData, ResponseBias_DisplayColor)

%Input
close all;

SaveGraph =1; 
SaveTable =1; 

roundValue = 3; NonParametri = 0; 
%% Colors For Plotting
Color.Ipsi = [0 0 1 ]; 
Color.Contra = [1 0 1 ]; 
% Size of Writing
% Size of Writing
fs = 28; % font size
MarkSize_GraphFAR_HR = 12;
MarkSize = 20;
if strcmp(experiment, 'Microstimulation')
    xlim_SDT_easy = [-1 4];
    xlim_SDT_diff = [-1 4];

else
xlim_SDT_easy = [1 5];
xlim_SDT_diff = [-1 4];
end

MarkSize_GraphFAR_HR_PerSession = 7; 
MarkSize_GraphFAR_HR = 17; 
MarkSize_CritDpr_small = 8; 
MarkSize_CritDpr = 20;
LineWidthSize = 3;


if strcmp(experiment, 'Inactivation') && isempty(folder_baseline) == 0
    ControlFolder = '_BaselinePost' ;
    xLegend = {'Ctr pst' 'Ina pst'};
elseif strcmp(experiment, 'Inactivation')&& isempty(folder_baseline) == 1
    ControlFolder = '_InactivationPre' ;
    xLegend = {'Ina pre' 'Ina pst'};
elseif strcmp(experiment, 'Microstimulation')
    ControlFolder = '' ;
    xLegend = {' Ctr' ' Sti'};
else
    ControlFolder = '' ;
    xLegend = {' pre' ' pst'};
end


%% load the data
if strcmp(experiment, 'Ephys')
    Files = dir([path_getData, folder , filesep, '*.mat' ]); 
    load([path_getData, folder, filesep, Files(end).name ])
    disp( ['analyzed dataset:  ',  filesep,folder, filesep, Files(end).name ])
    TimeLine        = {'stim_off', 'go_signal'};
    
    
elseif strcmp(experiment, 'Inactivation')
   Files = dir([path_getData, monkey ,'\behavior\' folder , filesep, '*.mat' ]); 
    load([path_getData, monkey ,'\behavior\' folder, filesep, Files(end).name ])
    disp( ['analyzed dataset:  ',  filesep,folder, filesep, Files(end).name ])

    TimeLine        = {'stim_off', 'go_signal'};
    if isempty(folder_baseline)== 0
        Bline =load([path_getData, monkey ,'\behavior\' folder_baseline, filesep, dataset_baseline ]);
    end
elseif strcmp(experiment, 'Microstimulation')
    load(['Y:\Projects\Pulv_distractor_spatial_choice\Microstimulation\Data',filesep, monkey , filesep, folder, filesep, dataset ])
    if strcmp(monkey, 'Curius')
        if  strcmp(folder, 'STIM_early_20161020_until_20161209')
            TimeLine        = {'stim_off', 'minus_80ms'};
            Stimulation = 'Early'; 
        elseif  strcmp(folder, 'STIM_late_20161020_until_20161209')
            TimeLine        = {'stim_off', 'plus_80ms'};
            Stimulation = 'Late'; 
        end
    elseif strcmp(monkey, 'Cornelius')
        if  strcmp(folder, 'STIM_early_20171126_until_20180123')
            TimeLine        = {'stim_off', 'minus_250ms'};
            Stimulation = 'Early'; 
        elseif  strcmp(folder, 'STIM_late_20180109_until_20180123') %STIM_late_20171126_until_20180123
            TimeLine        = {'stim_off', 'minus_50ms'};
            Stimulation = 'Late'; 
        end
    end
end

%% extract the data from the datafile
% single trials
%targets : successfull.. as saccaded (Hit) & unsuccessful as staying(Miss)
%contra = ipsi
% go_signal = pre
maxSess =  length(data_struct_per_session.targ_targ_2HF.combined.go_signal);
positionExp     = {'ipsi', 'contra'};
TimeLineExp         = {'pre', 'post'};
Difficulty          = {'easy', 'difficult'};
Stimuli = {'targ_targ', 'distr_distr'};
%data_struct_per_session.single_distr.(positionExp{indPos}).(TimeLine{indTime}){2,indSess}.G_R_ratios
for indSess = 1:   maxSess
    for indPos = 1:length(positionExp)
        for indTime = 1:length(TimeLine)
            if   strcmp(positionExp{indPos}, 'contra')
                Session.(positionExp{indPos}).(TimeLineExp{indTime})(indSess)               = data_struct_per_session.targ_targ_2HF.combined.(TimeLine{indTime}){1,indSess}.session;

                NumTotalTargets.(positionExp{indPos}).(TimeLineExp{indTime})(indSess)           = data_struct_per_session.targ_targ_2HF.combined.(TimeLine{indTime}){1,indSess}.num_total_trials;
                Hit.(positionExp{indPos}).(TimeLineExp{indTime})(indSess)                       = data_struct_per_session.targ_targ_2HF.combined.(TimeLine{indTime}){1,indSess}.num_choice_contra;
                Miss.(positionExp{indPos}).(TimeLineExp{indTime})(indSess)                      = data_struct_per_session.targ_targ_2HF.combined.(TimeLine{indTime}){1,indSess}.num_choice_fixation;
                
                NumTotalDistractor.easy.(positionExp{indPos}).(TimeLineExp{indTime})(indSess)   = data_struct_per_session.distr_distr_2HF.combined.(TimeLine{indTime}){2,indSess}.num_total_trials;
                NumTotalDistractor.difficult.(positionExp{indPos}).(TimeLineExp{indTime})(indSess) = data_struct_per_session.distr_distr_2HF.combined.(TimeLine{indTime}){1,indSess}.num_total_trials;
                CR.easy.(TimeLineExp{indTime})(indSess)                                         = data_struct_per_session.distr_distr_2HF.combined.(TimeLine{indTime}){2,indSess}.num_choice_fixation;
                CR.difficult.(TimeLineExp{indTime})(indSess)                                    = data_struct_per_session.distr_distr_2HF.combined.(TimeLine{indTime}){1,indSess}.num_choice_fixation;
                FA.easy.(positionExp{indPos}).(TimeLineExp{indTime})(indSess)                   = data_struct_per_session.distr_distr_2HF.combined.(TimeLine{indTime}){2,indSess}.num_choice_contra;
                FA.difficult.(positionExp{indPos}).(TimeLineExp{indTime})(indSess)              = data_struct_per_session.distr_distr_2HF.combined.(TimeLine{indTime}){1,indSess}.num_choice_contra;
            elseif     strcmp(positionExp{indPos}, 'ipsi')
                Session.(positionExp{indPos}).(TimeLineExp{indTime})(indSess)               = data_struct_per_session.targ_targ_2HF.combined.(TimeLine{indTime}){1,indSess}.session;

                NumTotalTargets.(positionExp{indPos}).(TimeLineExp{indTime})(indSess)               = data_struct_per_session.targ_targ_2HF.combined.(TimeLine{indTime}){1,indSess}.num_total_trials;
                Hit.(positionExp{indPos}).(TimeLineExp{indTime})(indSess)                           = data_struct_per_session.targ_targ_2HF.combined.(TimeLine{indTime}){1,indSess}.num_choice_ipsi;
                Miss.(positionExp{indPos}).(TimeLineExp{indTime})(indSess)                          = data_struct_per_session.targ_targ_2HF.combined.(TimeLine{indTime}){1,indSess}.num_choice_fixation;            
                NumTotalDistractor.easy.(positionExp{indPos}).(TimeLineExp{indTime})(indSess)       = data_struct_per_session.distr_distr_2HF.combined.(TimeLine{indTime}){2,indSess}.num_total_trials;
                NumTotalDistractor.difficult.(positionExp{indPos}).(TimeLineExp{indTime})(indSess)  = data_struct_per_session.distr_distr_2HF.combined.(TimeLine{indTime}){1,indSess}.num_total_trials;
                CR.easy.(TimeLineExp{indTime})(indSess)                                             = data_struct_per_session.distr_distr_2HF.combined.(TimeLine{indTime}){2,indSess}.num_choice_fixation;
                CR.difficult.(TimeLineExp{indTime})(indSess)                                         = data_struct_per_session.distr_distr_2HF.combined.(TimeLine{indTime}){1,indSess}.num_choice_fixation;
                FA.easy.(positionExp{indPos}).(TimeLineExp{indTime})(indSess)                       = data_struct_per_session.distr_distr_2HF.combined.(TimeLine{indTime}){2,indSess}.num_choice_ipsi;
                FA.difficult.(positionExp{indPos}).(TimeLineExp{indTime})(indSess)                  = data_struct_per_session.distr_distr_2HF.combined.(TimeLine{indTime}){1,indSess}.num_choice_ipsi;
            end
        end
    end
end


%% ONLY FOR INACTIVATION: extract the Data for the Baseline-session
if strcmp(experiment, 'Inactivation') && isempty(folder_baseline) == 0
    maxSess =  length(Bline.data_struct_per_session.targ_targ_2HF.combined.stim_off);
    positionExp     = {'ipsi', 'contra'};
    TimeLineExp         = {'pre', 'post'};
    Difficulty          = {'easy', 'difficult'};
    

    
    for indSess = 1:   maxSess
        for indPos = 1:length(positionExp)
            for indTime = 2 % use the Post-Control session to overwrite the pre
                if   strcmp(positionExp{indPos}, 'contra')
                    Ctr_Session.(positionExp{indPos}).(TimeLineExp{1})(indSess)               = Bline.data_struct_per_session.targ_targ_2HF.combined.(TimeLine{indTime}){1,indSess}.session;

                    NumTotalTargets.(positionExp{indPos}).(TimeLineExp{1})(indSess)               = Bline.data_struct_per_session.targ_targ_2HF.combined.(TimeLine{indTime}){1,indSess}.num_total_trials;
                    Hit.(positionExp{indPos}).(TimeLineExp{1})(indSess)                           = Bline.data_struct_per_session.targ_targ_2HF.combined.(TimeLine{indTime}){1,indSess}.num_choice_contra;
                    Miss.(positionExp{indPos}).(TimeLineExp{1})(indSess)                          = Bline.data_struct_per_session.targ_targ_2HF.combined.(TimeLine{indTime}){1,indSess}.num_choice_fixation;                  
                    NumTotalDistractor.easy.(positionExp{indPos}).(TimeLineExp{1})(indSess)       = Bline.data_struct_per_session.distr_distr_2HF.combined.(TimeLine{indTime}){2,indSess}.num_total_trials;
                    NumTotalDistractor.difficult.(positionExp{indPos}).(TimeLineExp{1})(indSess)  = Bline.data_struct_per_session.distr_distr_2HF.combined.(TimeLine{indTime}){1,indSess}.num_total_trials;
                    CR.easy.(TimeLineExp{1})(indSess)                                             = Bline.data_struct_per_session.distr_distr_2HF.combined.(TimeLine{indTime}){2,indSess}.num_choice_fixation;
                    CR.difficult.(TimeLineExp{1})(indSess)                                        = Bline.data_struct_per_session.distr_distr_2HF.combined.(TimeLine{indTime}){1,indSess}.num_choice_fixation;
                    FA.easy.(positionExp{indPos}).(TimeLineExp{1})(indSess)                       = Bline.data_struct_per_session.distr_distr_2HF.combined.(TimeLine{indTime}){2,indSess}.num_choice_contra;
                    FA.difficult.(positionExp{indPos}).(TimeLineExp{1})(indSess)                  = Bline.data_struct_per_session.distr_distr_2HF.combined.(TimeLine{indTime}){1,indSess}.num_choice_contra;
                    
                elseif     strcmp(positionExp{indPos}, 'ipsi')
                    Ctr_Session.(positionExp{indPos}).(TimeLineExp{1})(indSess)               = Bline.data_struct_per_session.targ_targ_2HF.combined.(TimeLine{indTime}){1,indSess}.session;

                    NumTotalTargets.(positionExp{indPos}).(TimeLineExp{1})(indSess)               =  Bline.data_struct_per_session.targ_targ_2HF.combined.(TimeLine{indTime}){1,indSess}.num_total_trials;
                    Hit.(positionExp{indPos}).(TimeLineExp{1})(indSess)                           =  Bline.data_struct_per_session.targ_targ_2HF.combined.(TimeLine{indTime}){1,indSess}.num_choice_ipsi;
                    Miss.(positionExp{indPos}).(TimeLineExp{1})(indSess)                          =  Bline.data_struct_per_session.targ_targ_2HF.combined.(TimeLine{indTime}){1,indSess}.num_choice_fixation;
                    
                    NumTotalDistractor.easy.(positionExp{indPos}).(TimeLineExp{1})(indSess)       =  Bline.data_struct_per_session.distr_distr_2HF.combined.(TimeLine{indTime}){2,indSess}.num_total_trials;
                    NumTotalDistractor.difficult.(positionExp{indPos}).(TimeLineExp{1})(indSess)  =  Bline.data_struct_per_session.distr_distr_2HF.combined.(TimeLine{indTime}){1,indSess}.num_total_trials;
                    CR.easy.(TimeLineExp{1})(indSess)                                             =  Bline.data_struct_per_session.distr_distr_2HF.combined.(TimeLine{indTime}){2,indSess}.num_choice_fixation;
                    CR.difficult.(TimeLineExp{1})(indSess)                                        =  Bline.data_struct_per_session.distr_distr_2HF.combined.(TimeLine{indTime}){1,indSess}.num_choice_fixation;
                    FA.easy.(positionExp{indPos}).(TimeLineExp{1})(indSess)                       =  Bline.data_struct_per_session.distr_distr_2HF.combined.(TimeLine{indTime}){2,indSess}.num_choice_ipsi;
                    FA.difficult.(positionExp{indPos}).(TimeLineExp{1})(indSess)                  =  Bline.data_struct_per_session.distr_distr_2HF.combined.(TimeLine{indTime}){1,indSess}.num_choice_ipsi;
%                     
                end
            end
        end
    end
    
    
end



%%  avoid 0 or Inf probabilities
%use the approach regardless wheather or not extreme values are obtained
for indPos = 1:length(positionExp) %contra vs ipsi selection
    for indTime = 1:length(TimeLine) %pre vs post
        Hit.(positionExp{indPos}).(TimeLineExp{indTime})            = Hit.(positionExp{indPos}).(TimeLineExp{indTime}) + 0.5 ;
        Miss.(positionExp{indPos}).(TimeLineExp{indTime})           = Miss.(positionExp{indPos}).(TimeLineExp{indTime}) + 0.5 ;
        FA.easy.(positionExp{indPos}).(TimeLineExp{indTime})        = FA.easy.(positionExp{indPos}).(TimeLineExp{indTime}) + 0.5 ;
        FA.difficult.(positionExp{indPos}).(TimeLineExp{indTime})   = FA.difficult.(positionExp{indPos}).(TimeLineExp{indTime}) + 0.5;
    end
end
        CR.easy.(TimeLineExp{indTime})                              = CR.easy.(TimeLineExp{indTime})  + 0.5;
        CR.difficult.(TimeLineExp{indTime})                         = CR.difficult.(TimeLineExp{indTime})  + 0.5;
%% calculate d'prime & criterion
% for ipsi vs contra
% control vs inactivation
% difficult vs easy

    for indPos = 1:length(positionExp) %contra vs ipsi selection
        for indTime = 1:length(TimeLine) %pre vs post
            %target is presented
            if indPos == 1 ; indPos2 = 2 ; else indPos2 = 1 ; end
            pHit.(positionExp{indPos}).(TimeLineExp{indTime})=    Hit.(positionExp{indPos}).(TimeLineExp{indTime}) ./ ...
                sum(Hit.(positionExp{indPos}).(TimeLineExp{indTime})  + Miss.(positionExp{indPos}).(TimeLineExp{indTime}) + Hit.(positionExp{indPos2}).(TimeLineExp{indTime}),1);
            %distractor is presented
            pFA.easy.(positionExp{indPos}).(TimeLineExp{indTime}) = FA.easy.(positionExp{indPos}).(TimeLineExp{indTime}) ./ ...
                sum(FA.easy.(positionExp{indPos}).(TimeLineExp{indTime})+  CR.easy.(TimeLineExp{indTime}) + FA.easy.(positionExp{indPos2}).(TimeLineExp{indTime}),1);
            
            pFA.difficult.(positionExp{indPos}).(TimeLineExp{indTime}) = FA.difficult.(positionExp{indPos}).(TimeLineExp{indTime}) ./ ...
                sum(FA.difficult.(positionExp{indPos}).(TimeLineExp{indTime})+  CR.difficult.(TimeLineExp{indTime})+FA.difficult.(positionExp{indPos2}).(TimeLineExp{indTime}),1);
            
              
          %distractor is presented
        pMiss.(TimeLineExp{indTime}) = Miss.ipsi.(TimeLineExp{indTime}) ./ ...
                sum(Hit.(positionExp{indPos}).(TimeLineExp{indTime})  + Miss.(positionExp{indPos}).(TimeLineExp{indTime}) + Hit.(positionExp{indPos2}).(TimeLineExp{indTime}),1);
        
        pCR.easy.(TimeLineExp{indTime}) = CR.easy.(TimeLineExp{indTime}) ./ ...
                sum(FA.easy.(positionExp{indPos}).(TimeLineExp{indTime})+  CR.easy.(TimeLineExp{indTime}) + FA.easy.(positionExp{indPos2}).(TimeLineExp{indTime}),1);
        pCR.difficult.(TimeLineExp{indTime}) = CR.difficult.(TimeLineExp{indTime}) ./ ...
                sum(FA.difficult.(positionExp{indPos}).(TimeLineExp{indTime})+  CR.difficult.(TimeLineExp{indTime}) + FA.difficult.(positionExp{indPos2}).(TimeLineExp{indTime}),1);
       
            
            [dprime.easy.(positionExp{indPos}).(TimeLineExp{indTime}),beta.easy.(positionExp{indPos}).(TimeLineExp{indTime}),criterion.easy.(positionExp{indPos}).(TimeLineExp{indTime})] = ...
                testsim_dprime( pHit.(positionExp{indPos}).(TimeLineExp{indTime}), pFA.easy.(positionExp{indPos}).(TimeLineExp{indTime}));
            [dprime.difficult.(positionExp{indPos}).(TimeLineExp{indTime}),beta.difficult.(positionExp{indPos}).(TimeLineExp{indTime}),criterion.difficult.(positionExp{indPos}).(TimeLineExp{indTime})] = ...
                testsim_dprime( pHit.(positionExp{indPos}).(TimeLineExp{indTime}), pFA.difficult.(positionExp{indPos}).(TimeLineExp{indTime}));
            
        end
    end
    
%
%% %%%%%%%%  SAVE THE DATASET %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T = table([nan;nan], [nan;nan], [nan;nan]) ; 
Table = [];
            Data.dprime = dprime ; 
            Data.criterion = criterion ; 
            Data.Session = Session;
            if strcmp(experiment, 'Inactivation'); Data.Ctr_Session = Ctr_Session; end

        all_DV = {'criterion', 'dprime'};
        for indVar = 1: length(all_DV)
            for indPos = 1: length(positionExp)
                for indTime = 1: length(TimeLineExp)
                    for  indDiff = 1: length(Difficulty)
                        Values = []; Hemifield = []; TimeLineExperiment = []; StimulusTyp = []; DisDifficulty = []; Monkey = [];  DependentVariable = []; Date = []; 
                        Experiment = [];
                        
                        Values = Data.(all_DV{indVar}).(Difficulty {indDiff}).(positionExp {indPos}).(TimeLineExp {indTime})';
                        Hemifield       = repmat(positionExp(indPos), length(Values),1);
                        TimeLineExperiment     =  repmat(TimeLineExp(indTime), length(Values),1);
                        StimulusTyp     = repmat({'DoubleSameStim'}, length(Values),1);
                        DisDifficulty   = repmat(Difficulty (indDiff), length(Values),1);
                        Monkey          = repmat(monkey, length(Values),1);

                        if strcmp(experiment, 'Inactivation')&& strcmp((TimeLineExp {indTime}), 'post')
                            Date        = Data.Session.(positionExp{indPos}).(TimeLineExp{indTime})';
                            Experiment    = repmat({'Ina'}, length(Date),1);
                            
                        elseif strcmp(experiment, 'Inactivation')&& strcmp((TimeLineExp {indTime}), 'pre')
                            TimeLineExperiment = repmat({'post'}, length(Values),1);
                            Date        = Data.Ctr_Session.(positionExp {indPos}).(TimeLineExp {indTime})';
                            Experiment    = repmat({'Ctr'}, length(Date),1);

                        elseif strcmp(experiment, 'Microstimulation') 
                            Date        = Data.Session.(positionExp {indPos}).(TimeLineExp {indTime})';
                            if strcmp(Stimulation, 'Late')
                            Experiment    = repmat({'Late'}, length(Date),1);
                            elseif strcmp(Stimulation, 'Early')
                            Experiment    = repmat({'Early'}, length(Date),1);  
                            end
                            if  strcmp((TimeLineExp {indTime}), 'pre')
                            TimeLineExperiment = repmat({'off'}, length(Values),1);
                            else
                            TimeLineExperiment = repmat(TimeLine(2), length(Values),1);  
                            end
                        end
                        
                        DependentVariable   = repmat((all_DV(indVar)), length(Values),1);
                        
                        Tab = table(Monkey, DisDifficulty, Values, Hemifield ,TimeLineExperiment,  StimulusTyp,Experiment,Date,DependentVariable );
                        Table = [ Table; Tab ];
                        Tab = [];
                    end
                end
            end
        end

if strcmp(experiment, 'Inactivation') && SaveTable == 1;        
%addtoDropbox = 'C:\Users\kkaduk\Dropbox\PhD\Projects\InaDPul_SDT\AGit_DPul_SDT\AGit_InaDPul_SDT\Data\'; 
 filename =[path_SaveFig,['Table', ControlFolder], filesep, monkey, '_SDTvar_DoubleSameStim' ] ;  

writetable(Table, filename , 'Delimiter', ' ')
disp(['SAVED   ', filename])

elseif SaveTable == 1
addtoDropbox =  'C:\Users\kkaduk\Dropbox\PhD\Projects\InaDPul_SDT\AGit_DPul_SDT\AGit_InaDPul_SDT\Data\'; 
%save([addtoDropbox,  monkey, filesep,'SDTvariables_Microstimulation', filesep ,monkey, '_SDTvar_DoubleSameStim' ],'Table');
writetable(Table, [addtoDropbox, monkey, filesep,'SDTvariables_Microstimulation', filesep, experiment, filesep,monkey,'_', Stimulation, '_SDTvar_DoubleSameStim' ], 'Delimiter', ' ')
disp(['SAVED   ', addtoDropbox, monkey, filesep,monkey, '_SDTvar_DoubleSameStim' ])

end



%% Statistic

TablePwerte = [];

if strcmp(experiment, 'Inactivation')&& isempty(folder_baseline) == 0 && NonParametri == 0
    %statistic - parametric: independent sessions
    for indPos = 1:length(positionExp) %contra vs ipsi selection
            [h,p,ci,stats]  = ttest2(dprime.easy.(positionExp{indPos}).pre, dprime.easy.(positionExp{indPos}).post);
            dprime.easy.(positionExp{indPos}).pvalue = p;
            dprime.easy.(positionExp{indPos}).tstat = stats.tstat;   
            disp(['dprime ', 'easy ' , (positionExp{indPos}) '_',num2str(round(p,2)),'_t=', num2str(round(dprime.easy.(positionExp{indPos}).tstat,2))] )

            TabPwerte = [ ];
            TabPwerte = table({monkey}, {'DoubleSameStim'}, {'dprime'},{ 'easy'} , (positionExp(indPos)) ,round(p,roundValue),{'t'} , round(dprime.easy.(positionExp{indPos}).tstat,2), {'indepTtest '} );
            TablePwerte = [ TablePwerte; TabPwerte ];
            
            [h,p,ci,stats]  = ttest2(criterion.easy.(positionExp{indPos}).pre, criterion.easy.(positionExp{indPos}).post);
            criterion.easy.(positionExp{indPos}).pvalue = p;
            criterion.easy.(positionExp{indPos}).tstat = stats.tstat;   
            disp(['criterion ', 'easy ' , (positionExp{indPos}) '_',num2str(round(p,2)),'_t=', num2str(round(criterion.easy.(positionExp{indPos}).tstat,2))] )
            
            TabPwerte = [ ];
            TabPwerte = table({monkey }, {'DoubleSameStim'}, {'criterion'},{ 'easy'} , (positionExp(indPos)) ,round(p,roundValue),{'t'} , round(criterion.easy.(positionExp{indPos}).tstat,2), {'indepTtest '} );
            TablePwerte = [ TablePwerte; TabPwerte ];
            
            [h,p,ci,stats]  = ttest2(dprime.difficult.(positionExp{indPos}).pre, dprime.difficult.(positionExp{indPos}).post);
            dprime.difficult.(positionExp{indPos}).pvalue = p;
            dprime.difficult.(positionExp{indPos}).tstat = stats.tstat;   
            disp(['dprime ', 'difficult ' , (positionExp{indPos}) '_',num2str(round(p,2)),'_t=', num2str(round(dprime.difficult.(positionExp{indPos}).tstat,2))] )

            TabPwerte = [ ];
            TabPwerte = table({monkey}, {'DoubleSameStim'}, {'dprime'},{ 'difficult'} , (positionExp(indPos)) ,round(p,roundValue),{'t'} , round(dprime.difficult.(positionExp{indPos}).tstat,2), {'indepTtest '} );
            TablePwerte = [ TablePwerte; TabPwerte ];

            [h,p,ci,stats]  = ttest2(criterion.difficult.(positionExp{indPos}).pre, criterion.difficult.(positionExp{indPos}).post);
            criterion.difficult.(positionExp{indPos}).pvalue = p;
            criterion.difficult.(positionExp{indPos}).tstat = stats.tstat;   
            disp(['criterion ', 'difficult ' , (positionExp{indPos}) '_',num2str(round(p,2)),'_t=', num2str(round(criterion.difficult.(positionExp{indPos}).tstat,2))] )

            TabPwerte = [ ];
            TabPwerte = table({monkey}, {'DoubleSameStim'}, {'criterion'},{ 'difficult'} , (positionExp(indPos)) ,round(p,roundValue),{'t'} , round(criterion.difficult.(positionExp{indPos}).tstat,2), {'indepTtest '} );
            TablePwerte = [ TablePwerte; TabPwerte ];
           
      
          [h,p,ci,stat] = ttest2(pFA.easy.(positionExp{indPos}).pre, pFA.easy.(positionExp{indPos}).post);
        pFA.easy.(positionExp{indPos}).pvalue = p;
        pFA.easy.(positionExp{indPos}).tstat = stat.tstat;
        median_Pre  = median(pFA.easy.(positionExp{indPos}).pre) ; 
        median_Post = median(pFA.easy.(positionExp{indPos}).post) ;
        pFA.easy.(positionExp{indPos}).r = norminv(p)/sqrt(length(pFA.easy.(positionExp{indPos}).post) + length(pFA.easy.(positionExp{indPos}).pre)); 
        disp([median_Pre ,' ',  median_Post, 'pFA ', 'easy ' , (positionExp{indPos}) '_',num2str(round(p,roundValue)),'_t=', num2str(round( pFA.easy.(positionExp{indPos}).tstat,2))] )

            TabPwerte = [ ];
            TabPwerte = table({monkey}, {'DoubleSameStim'}, {'pFA'},{ 'easy'} , (positionExp(indPos)) ,round(p,roundValue),{'t'} , round(pFA.easy.(positionExp{indPos}).tstat,3), {'indepTtest2'} );
            TablePwerte = [ TablePwerte; TabPwerte ];
       
        
 
        [h,p,ci,stat] = ttest2(pFA.difficult.(positionExp{indPos}).pre, pFA.difficult.(positionExp{indPos}).post);
        pFA.difficult.(positionExp{indPos}).pvalue = p;
        pFA.difficult.(positionExp{indPos}).tstat = stat.tstat;
        pFA.difficult.(positionExp{indPos}).r = norminv(p)/sqrt(length(pFA.difficult.(positionExp{indPos}).post) + length(pFA.difficult.(positionExp{indPos}).pre)); 
        median_Pre  = median(pFA.difficult.(positionExp{indPos}).pre) ; 
        median_Post = median(pFA.difficult.(positionExp{indPos}).post) ;
        disp([median_Pre ,' ',  median_Post, 'pFA ', 'difficult ' , (positionExp{indPos}) '_',num2str(round(p,roundValue)),'_t=', num2str(round( pFA.difficult.(positionExp{indPos}).tstat,2))] )
            TabPwerte = [ ];
            TabPwerte = table({monkey}, {'DoubleSameStim'}, {'pFA'},{ 'difficult'} , (positionExp(indPos)) ,round(p,roundValue),{'t'} , round(pFA.difficult.(positionExp{indPos}).tstat,3), {'indepTtest2'} );
            TablePwerte = [ TablePwerte; TabPwerte ];
       

        
        [h,p,ci,stat] = ttest2(pHit.(positionExp{indPos}).pre, pHit.(positionExp{indPos}).post);
        pHit.(positionExp{indPos}).pvalue = p;
        pHit.(positionExp{indPos}).tstat = stat.tstat;
        pHit.(positionExp{indPos}).r = norminv(p)/sqrt(length(pHit.(positionExp{indPos}).post) + length(pHit.(positionExp{indPos}).pre)); 
        median_Pre  = median(pHit.(positionExp{indPos}).pre) ; 
        median_Post = median(pHit.(positionExp{indPos}).post) ;
        disp([median_Pre ,' ',  median_Post, 'pHi '  (positionExp{indPos}) '_',num2str(round(p,2)),'_t=', num2str(round(pHit.(positionExp{indPos}).tstat,2))] )
            TabPwerte = [ ];
            TabPwerte = table({monkey}, {'DoubleSameStim'}, {'pHi'},{ 'nan'} , (positionExp(indPos)) ,round(p,roundValue),{'t'} , round(pHit.(positionExp{indPos}).tstat,3), {'indepTtest2'} );
            TablePwerte = [ TablePwerte; TabPwerte ];
       

        
        [h,p,ci,stat] = ttest2(pMiss.pre, pMiss.post);
        pMiss.pvalue = p;
        [h,p,ci,stat] = ttest2(pCR.easy.pre,pCR.easy.post);
        pCR.easy.pvalue = p;

        [h,p,ci,stat] = ttest2(pCR.difficult.pre,pCR.difficult.post);
        pCR.difficult.pvalue = p;
    
    end   
    

else
     %statistic - parametric: dependent sessions
    for indPos = 1:length(positionExp) %contra vs ipsi selection
            [p,h,stat] = ranksum(dprime.easy.(positionExp{indPos}).pre, dprime.easy.(positionExp{indPos}).post);
            dprime.easy.(positionExp{indPos}).pvalue = p;
           
            TabPwerte = [ ];
            TabPwerte = table({monkey }, {'DoubleSameStim'}, {'dprime'},{ 'easy'} , (positionExp(indPos)) ,round(p,roundValue),{'U'} , stat.ranksum, {'ranksum '} );
            TablePwerte = [ TablePwerte; TabPwerte ];
            disp(p)
            p = [ ];
         
            [p,h,stat] = ranksum( criterion.easy.(positionExp{indPos}).post, criterion.easy.(positionExp{indPos}).pre);
            criterion.easy.(positionExp{indPos}).pvalue = p;
            
             TabPwerte = [ ];
            TabPwerte = table({monkey }, {'DoubleSameStim'}, {'criterion'},{ 'easy'} , (positionExp(indPos)) ,round(p,roundValue),{'U'} ,stat.ranksum, {'ranksum '} );
            TablePwerte = [ TablePwerte; TabPwerte ];
                      disp(p); p = [ ];
            [p,h,stat] = ranksum(dprime.difficult.(positionExp{indPos}).pre, dprime.difficult.(positionExp{indPos}).post);
            dprime.difficult.(positionExp{indPos}).pvalue = p;
            
             TabPwerte = [ ];
            TabPwerte = table({monkey }, {'DoubleSameStim'}, {'dprime'},{ 'difficult'} , (positionExp(indPos)) ,round(p,roundValue),{'U'} ,stat.ranksum, {'ranksum '} );
            TablePwerte = [ TablePwerte; TabPwerte ];
                      disp(p); p = [ ];
            [p,h,stat] = ranksum( criterion.difficult.(positionExp{indPos}).post, criterion.difficult.(positionExp{indPos}).pre);
            criterion.difficult.(positionExp{indPos}).pvalue = p;
             TabPwerte = [ ];
            TabPwerte = table({monkey }, {'DoubleSameStim'}, {'criterion'},{ 'difficult'} , (positionExp(indPos)) ,round(p,roundValue),{'U'} , stat.ranksum, {'ranksum '} );
            TablePwerte = [ TablePwerte; TabPwerte ];            
            
        [H,P,CI,STATS] = ttest(pFA.easy.(positionExp{indPos}).pre, pFA.easy.(positionExp{indPos}).post);
        pFA.easy.(positionExp{indPos}).pvalue = P;
        median_Pre  = median(pFA.easy.(positionExp{indPos}).pre) ; 
        median_Post = median(pFA.easy.(positionExp{indPos}).post) ;
       pFA.easy.(positionExp{indPos}).r = norminv(p)/sqrt(length(pFA.easy.(positionExp{indPos}).post) + length(pFA.easy.(positionExp{indPos}).pre)); 
        disp([median_Pre ,' ',  median_Post, 'pFA ', 'easy ' , (positionExp{indPos}) '_',num2str(round(P,roundValue)),'_r=', num2str(round(pFA.easy.(positionExp{indPos}).r,2))] )

            TabPwerte = [ ];
            TabPwerte = table({monkey}, {'DoubleSameStim'}, {'pFA'},{ 'easy'} , (positionExp(indPos)) ,round(P,roundValue),{'t'} , round(STATS.tstat,3), {'depTtest'} );
            TablePwerte = [ TablePwerte; TabPwerte ];
 
        [H,P,CI,STATS] = ttest(pFA.difficult.(positionExp{indPos}).pre, pFA.difficult.(positionExp{indPos}).post);
        pFA.difficult.(positionExp{indPos}).pvalue = P;
        median_Pre  = median(pFA.difficult.(positionExp{indPos}).pre) ; 
        median_Post = median(pFA.difficult.(positionExp{indPos}).post) ;
         pFA.difficult.(positionExp{indPos}).r = norminv(p)/sqrt(length(pFA.difficult.(positionExp{indPos}).post) + length(pFA.difficult.(positionExp{indPos}).pre)); 

        disp([median_Pre ,' ',  median_Post, 'pFA ', 'difficult ' , (positionExp{indPos}) '_',num2str(round(P,roundValue)),'_r=', num2str(round(pFA.difficult.(positionExp{indPos}).r,3))] )
            TabPwerte = [ ];
            TabPwerte = table({monkey}, {'DoubleSameStim'}, {'pFA'},{ 'difficult'} , (positionExp(indPos)) ,round(P,roundValue),{'t'} , round(STATS.tstat,3), {'depTtest'} );
            TablePwerte = [ TablePwerte; TabPwerte ];
 
        
       [H,P,CI,STATS] = ttest(pHit.(positionExp{indPos}).pre, pHit.(positionExp{indPos}).post);
        pHit.(positionExp{indPos}).pvalue = P;
        median_Pre  = median(pHit.(positionExp{indPos}).pre) ; 
        median_Post = median(pHit.(positionExp{indPos}).post) ;
        pHit.(positionExp{indPos}).r = norminv(p)/sqrt(length(pHit.(positionExp{indPos}).post) + length(pHit.(positionExp{indPos}).pre)); 

        disp([median_Pre ,' ',  median_Post, 'pHi '  (positionExp{indPos}) '_',num2str(round(P,roundValue)),'_r=', num2str(round(pHit.(positionExp{indPos}).r,2))] )
        TabPwerte = [ ];
            TabPwerte = table({monkey}, {'DoubleSameStim'}, {'pHi'},{ 'no'} , (positionExp(indPos)) ,round(P,roundValue),{'t'} , round(STATS.tstat,3), {'depTtest'} );
            TablePwerte = [ TablePwerte; TabPwerte ];
 
        
        [H,P,CI,STATS] = ttest(pMiss.pre, pMiss.post);
        pMiss.pvalue = P;
       [H,P,CI,STATS] = ttest2(pCR.easy.pre,pCR.easy.post);
        pCR.easy.pvalue = P;

        [H,P,CI,STATS] = ttest(pCR.difficult.pre,pCR.difficult.post);
        pCR.difficult.pvalue =P;
       

    end
end

%addtoDropbox = 'C:\Users\kkaduk\Dropbox\PhD\Projects\InaDPul_SDT\AGit_DPul_SDT\AGit_InaDPul_SDT\Data\'; 
if strcmp(experiment, 'Inactivation') && SaveTable
 filename =[path_SaveFig,['Table', ControlFolder], filesep, monkey, '_SDTvarPwerte_DoubleSameStim.xlsx' ] ;  
 %filename =[addtoDropbox,filesep, monkey,filesep,   'SDTvariables' , filesep, experiment, filesep, monkey, '_SDTvarPwerte_DoubleSameStim.xlsx' ] ;  
writetable(TablePwerte,filename,'Sheet',1,  'Range' ,'A1' )

disp(['SAVED  ' filename ])
   
elseif strcmp(experiment, 'Microstimulation') && SaveTable
 filename =[path_SaveFig,['Table', ControlFolder], filesep, monkey, '_SDTvarPwerte_DoubleSameStim.xlsx' ] ;  
%    filename = [addtoDropbox ,  monkey, filesep, 'SDTvariables' , filesep, experiment, filesep, 'Stat' ,filesep, monkey,'_', Stimulation, '_SDTvarPwerte_DoubleSameStim.xlsx' ] ;  
writetable(TablePwerte,filename,'Sheet',1,  'Range' ,'A1' )

disp(['SAVED  ' , path_SaveFig ,filesep,  'Stat' ,filesep, monkey,'_' , Stimulation,  '_SDTvarPwerte_DoubleSameStim' ])

end

%% Determine the if the dprime differs against the H0?
for indPos = 1:length(positionExp) %contra vs ipsi selection
    for indTime = 1:length(TimeLine) %pre vs post
        [p,h,stat] = signrank(dprime.easy.(positionExp{indPos}).pre);
        dprime.easy.(positionExp{indPos}).P_Bias_pre = p;
        [p,h,stat] = signrank(dprime.easy.(positionExp{indPos}).post);
        dprime.easy.(positionExp{indPos}).P_Bias_post = p;
        [p,h,stat] = signrank(dprime.difficult.(positionExp{indPos}).pre);
        dprime.difficult.(positionExp{indPos}).P_Bias_pre = p;
        [p,h,stat] = signrank( dprime.difficult.(positionExp{indPos}).post);
        dprime.difficult.(positionExp{indPos}).P_Bias_post = p;
    end
end
%% Determine the response Bias by testing the criterion again the H0?
for indPos = 1:length(positionExp) %contra vs ipsi selection
    for indTime = 1:length(TimeLine) %pre vs post
        [p,h,stat] = signrank(criterion.easy.(positionExp{indPos}).pre);
        criterion.easy.(positionExp{indPos}).P_Bias_pre = p;
        [p,h,stat] = signrank(criterion.easy.(positionExp{indPos}).post);
        criterion.easy.(positionExp{indPos}).P_Bias_post = p;
        [p,h,stat] = signrank(criterion.difficult.(positionExp{indPos}).pre);
        criterion.difficult.(positionExp{indPos}).P_Bias_pre = p;
        [p,h,stat] = signrank( criterion.difficult.(positionExp{indPos}).post);
        criterion.difficult.(positionExp{indPos}).P_Bias_post = p;
    end
end
%% CRITERION - WHICH RESPONSE BIAS?


for indPos = 1:length(positionExp) %contra vs ipsi selection
    if ResponseBias_DisplayColor == 1
        Color_MoreContra    = [1 0 1]; %cyan
        Color_LessContra    = [0 0 1]; %magenta
        Color_Neutral       = [0.5 0.5 0.5];
    elseif ResponseBias_DisplayColor == 0
        Color_MoreContra    = [0 0 0];
        Color_LessContra    = [0 0 0];
        Color_Neutral       = [0 0 0];
    end
    for indDiff = 1: length(Difficulty)
        if  criterion.(Difficulty{indDiff}).(positionExp{indPos}).P_Bias_pre < 0.05
            if strcmp(positionExp{indPos}, 'contra')
                %
                if nanmean(criterion.(Difficulty{indDiff}).(positionExp{indPos}).pre) < 0 %negative Criterion
                    criterion.(Difficulty{indDiff}).(positionExp{indPos}).BiasDirection_pre = 'moreContra' ;
                    criterion.(Difficulty{indDiff}).(positionExp{indPos}).Bias_pre = 'Go' ;
                    criterion.(Difficulty{indDiff}).(positionExp{indPos}).Bias_pre_Color = Color_MoreContra ;
                    
                elseif nanmean(criterion.(Difficulty{indDiff}).(positionExp{indPos}).pre) > 0
                    criterion.(Difficulty{indDiff}).(positionExp{indPos}).BiasDirection_pre = 'lessContra' ;
                    criterion.(Difficulty{indDiff}).(positionExp{indPos}).Bias_pre = 'NoGo' ;
                    criterion.(Difficulty{indDiff}).(positionExp{indPos}).Bias_pre_Color = Color_LessContra ;
                    
                end
                
            elseif  strcmp(positionExp{indPos}, 'ipsi')
                % the ipsilateral Criterion is reversed .. positive
                % Criterion means after reversion 
                if nanmean(-criterion.(Difficulty{indDiff}).(positionExp{indPos}).pre) < 0 %negative Criterion
                    criterion.(Difficulty{indDiff}).(positionExp{indPos}).BiasDirection_pre = 'moreContra' ;
                    criterion.(Difficulty{indDiff}).(positionExp{indPos}).Bias_pre = 'NoGo' ;
                    criterion.(Difficulty{indDiff}).(positionExp{indPos}).Bias_pre_Color = Color_MoreContra ;
                    
                elseif nanmean(-criterion.(Difficulty{indDiff}).(positionExp{indPos}).pre) > 0
                    criterion.(Difficulty{indDiff}).(positionExp{indPos}).BiasDirection_pre = 'lessContra' ;
                    criterion.(Difficulty{indDiff}).(positionExp{indPos}).Bias_pre = 'Go' ;
                    criterion.(Difficulty{indDiff}).(positionExp{indPos}).Bias_pre_Color = Color_LessContra ;
                    
                end
                
            end
        else
            criterion.(Difficulty{indDiff}).(positionExp{indPos}).BiasDirection_pre = 'neutral' ;
            criterion.(Difficulty{indDiff}).(positionExp{indPos}).Bias_pre = 'neutral' ;
            criterion.(Difficulty{indDiff}).(positionExp{indPos}).Bias_pre_Color = Color_Neutral ;
            
        end
    end
end
%%% POST
for indPos = 1:length(positionExp) %contra vs ipsi selection
    if ResponseBias_DisplayColor == 1
        Color_MoreContra    = [1 0 1]; %cyan
        Color_LessContra    = [0 0 1]; %magenta
        Color_Neutral       = [0.5 0.5 0.5];
    elseif ResponseBias_DisplayColor == 0
        Color_MoreContra    = Color.Ipsi;
        Color_LessContra    = Color.Ipsi;
        Color_Neutral       = Color.Ipsi;
    end
    for indDiff = 1: length(Difficulty)
        if  criterion.(Difficulty{indDiff}).(positionExp{indPos}).P_Bias_post < 0.05
            if strcmp(positionExp{indPos}, 'contra')
                if nanmean(criterion.(Difficulty{indDiff}).(positionExp{indPos}).post) < 0 %negative Criterion
                    criterion.(Difficulty{indDiff}).(positionExp{indPos}).BiasDirection_post = 'moreContra' ;
                    criterion.(Difficulty{indDiff}).(positionExp{indPos}).Bias_post = 'Go' ;
                    criterion.(Difficulty{indDiff}).(positionExp{indPos}).Bias_post_Color = Color_MoreContra ;
                    
                elseif nanmean(criterion.(Difficulty{indDiff}).(positionExp{indPos}).post) > 0
                    criterion.(Difficulty{indDiff}).(positionExp{indPos}).BiasDirection_post = 'lessContra' ;
                    criterion.(Difficulty{indDiff}).(positionExp{indPos}).Bias_post = 'NoGo' ;
                    criterion.(Difficulty{indDiff}).(positionExp{indPos}).Bias_post_Color = Color_LessContra ;
                    
                end
            elseif  strcmp(positionExp{indPos}, 'ipsi')
                if nanmean(-criterion.(Difficulty{indDiff}).(positionExp{indPos}).post) < 0 %negative Criterion
                    criterion.(Difficulty{indDiff}).(positionExp{indPos}).BiasDirection_post = 'moreContra' ;
                    criterion.(Difficulty{indDiff}).(positionExp{indPos}).Bias_post = 'NoGo' ;
                    criterion.(Difficulty{indDiff}).(positionExp{indPos}).Bias_post_Color = Color_MoreContra ;
                    
                elseif nanmean(-criterion.(Difficulty{indDiff}).(positionExp{indPos}).post) > 0
                    criterion.(Difficulty{indDiff}).(positionExp{indPos}).BiasDirection_post = 'lessContra' ;
                    criterion.(Difficulty{indDiff}).(positionExp{indPos}).Bias_post = 'Go' ;
                    criterion.(Difficulty{indDiff}).(positionExp{indPos}).Bias_post_Color = Color_LessContra ;
                    
                end
                
            end
        else
            criterion.(Difficulty{indDiff}).(positionExp{indPos}).BiasDirection_post = 'neutral' ;
            criterion.(Difficulty{indDiff}).(positionExp{indPos}).Bias_post = 'neutral' ;
            criterion.(Difficulty{indDiff}).(positionExp{indPos}).Bias_post_Color = Color_Neutral  ;
        end
    end
end

%% GRAPH


%% Criterion vs Dprime
figure('Position',[200 200 1200 900],'PaperPositionMode','auto'); % ,'PaperOrientation','landscape'
set(gcf,'Name','Dprime vs Criterion');
set(gcf,'Color',[1 1 1]);


   
ha(indPos) = subplot(1,2,1);
for i = 1: length(dprime.easy.contra.pre)
plot([dprime.easy.contra.pre(i),dprime.easy.contra.post(i)],[criterion.easy.contra.pre(i),criterion.easy.contra.post(i)], 'o','color',[Color.Contra, 0.4] , 'MarkerSize',MarkSize-8,'markerfacecolor',[1 1 1 ]); hold on;
plot([dprime.easy.ipsi.pre(i),dprime.easy.ipsi.post(i)], [-criterion.easy.ipsi.pre(i),-criterion.easy.ipsi.post(i)], 'o','color',[Color.Ipsi, 0.4] , 'MarkerSize',MarkSize-8,'markerfacecolor',[1 1 1 ]); hold on;% reverse direction of criterion for ipsi
plot(dprime.easy.contra.post(i),criterion.easy.contra.post(i), 'o','color',Color.Contra , 'MarkerSize',MarkSize-10,'markerfacecolor',Color.Contra); hold on;
plot(dprime.easy.ipsi.post(i),-criterion.easy.ipsi.post(i), 'o','color',Color.Ipsi , 'MarkerSize',MarkSize-10,'markerfacecolor',Color.Ipsi); hold on;% reverse direction of criterion for ipsi
end
plot([nanmean(dprime.easy.contra.pre),nanmean(dprime.easy.contra.post) ],[nanmean(criterion.easy.contra.pre),nanmean(criterion.easy.contra.post)], 'o-','color',Color.Contra , 'MarkerSize',MarkSize,'markerfacecolor',[1 1 1 ]); hold on;
plot([nanmean(dprime.easy.ipsi.pre),nanmean(dprime.easy.ipsi.post)],[-nanmean(criterion.easy.ipsi.pre) ,-nanmean(criterion.easy.ipsi.post)], 'o-','color',Color.Ipsi , 'MarkerSize',MarkSize,'markerfacecolor',[1 1 1 ]); hold on;% reverse direction of criterion for ipsi

plot(nanmean(dprime.easy.contra.post),nanmean(criterion.easy.contra.post), 'o','color',Color.Contra , 'MarkerSize',MarkSize,'markerfacecolor',Color.Contra); hold on;
plot(nanmean(dprime.easy.ipsi.post),-nanmean(criterion.easy.ipsi.post), 'o','color',Color.Ipsi , 'MarkerSize',MarkSize,'markerfacecolor',Color.Ipsi); hold on;% reverse direction of criterion for ipsi

if dprime.easy.ipsi.pvalue < 0.05
    ymax = min(-nanmean(criterion.easy.ipsi.pre) ,-nanmean([criterion.easy.ipsi.post])) -0.3;
    ext_sigline([nanmean([dprime.easy.ipsi.pre]),nanmean([dprime.easy.ipsi.post])],dprime.easy.ipsi.pvalue,[ ],(ymax ),'x', Color.Ipsi); hold on;
end
if criterion.easy.ipsi.pvalue < 0.05
    y1 = -1*nanmean(criterion.easy.ipsi.pre) ;
    y2 = -1*nanmean([criterion.easy.ipsi.post]);
    ymax = max(nanmean(dprime.easy.ipsi.pre) ,nanmean([dprime.easy.ipsi.post])) ;
    ext_sigline([y1,y2],criterion.easy.ipsi.pvalue,[],ymax+ 0.4,'y', Color.Ipsi); hold on;
end


if dprime.easy.contra.pvalue < 0.05
    ymax = max(nanmean(criterion.easy.contra.pre) ,nanmean([criterion.easy.contra.post]));
    ext_sigline([nanmean([dprime.easy.contra.pre]),nanmean([dprime.easy.contra.post])],dprime.easy.contra.pvalue,[ ],(ymax +0.2),'x', Color.Contra); hold on;
end
if criterion.easy.contra.pvalue < 0.05
    y1 = nanmean(criterion.easy.contra.pre) ;
    y2 = nanmean([criterion.easy.contra.post]);
    ymax = max(nanmean(dprime.easy.contra.pre) ,nanmean([dprime.easy.contra.post])) ;
    ext_sigline([y1,y2],criterion.easy.contra.pvalue,[],ymax +0.4,'y', Color.Contra); hold on;
end
axis square
xlabel('sensitivity','fontsize',fs,'fontweight','b', 'Interpreter', 'none')
ylabel('criterion','fontsize',fs,'fontweight','b', 'Interpreter', 'none')
set(gca,'ylim',[-3 3],'xlim',[-1 5],'fontsize',fs)
title(' easy distractor')
text(0.5,-2.8, 'more Contra & ipsi:NoGo, contra:Go', 'Color', 'k')
text(0.5,2.8, 'less Contra & ipsi:Go, contra:NoGo', 'Color', 'k')

ha(indPos) = subplot(1,2,2);
for i = 1: length(dprime.difficult.contra.pre)
plot([dprime.difficult.contra.pre(i),dprime.difficult.contra.post(i)],[criterion.difficult.contra.pre(i),criterion.difficult.contra.post(i)], 'o','color',[Color.Contra, 0.4] , 'MarkerSize',MarkSize-8,'markerfacecolor',[1 1 1 ]); hold on;
plot([dprime.difficult.ipsi.pre(i),dprime.difficult.ipsi.post(i)], [-criterion.difficult.ipsi.pre(i),-criterion.difficult.ipsi.post(i)], 'o','color',[Color.Ipsi, 0.4] , 'MarkerSize',MarkSize-8,'markerfacecolor',[1 1 1 ]); hold on;% reverse direction of criterion for ipsi
plot(dprime.difficult.contra.post(i),criterion.difficult.contra.post(i), 'o','color',Color.Contra , 'MarkerSize',MarkSize-10,'markerfacecolor',Color.Contra); hold on;
plot(dprime.difficult.ipsi.post(i),-criterion.difficult.ipsi.post(i), 'o','color',Color.Ipsi , 'MarkerSize',MarkSize-10,'markerfacecolor',Color.Ipsi); hold on;% reverse direction of criterion for ipsi
end
plot([nanmean(dprime.difficult.contra.pre),nanmean(dprime.difficult.contra.post)],[nanmean(criterion.difficult.contra.pre) ,nanmean(criterion.difficult.contra.post)], 'o-','color',Color.Contra , 'MarkerSize',MarkSize,'markerfacecolor',[1 1 1 ]); hold on;
plot([nanmean(dprime.difficult.ipsi.pre),nanmean(dprime.difficult.ipsi.post)],[ -nanmean(criterion.difficult.ipsi.pre) ,-nanmean(criterion.difficult.ipsi.post)], 'o-','color',Color.Ipsi , 'MarkerSize',MarkSize,'markerfacecolor',[1 1 1 ]); hold on;% reverse direction of criterion for ipsi

plot(nanmean(dprime.difficult.contra.post),nanmean(criterion.difficult.contra.post), 'o','color',Color.Contra , 'MarkerSize',MarkSize,'markerfacecolor',[Color.Contra]); hold on;
plot(nanmean(dprime.difficult.ipsi.post),-nanmean(criterion.difficult.ipsi.post), 'o','color',Color.Ipsi , 'MarkerSize',MarkSize,'markerfacecolor',[Color.Ipsi ]); hold on;% reverse direction of criterion for ipsi

if dprime.difficult.ipsi.pvalue < 0.05
    ymax = max(nanmean(criterion.difficult.ipsi.pre) ,nanmean([criterion.difficult.ipsi.post]));
    ext_sigline([nanmean([dprime.difficult.ipsi.pre]),nanmean([dprime.difficult.ipsi.post])],dprime.difficult.ipsi.pvalue,[ ],(ymax -1 ),'x', Color.Ipsi); hold on;
end
if criterion.difficult.ipsi.pvalue < 0.05
    y1 = -1*nanmean(criterion.difficult.ipsi.pre) ;
    y2 = -1*nanmean([criterion.difficult.ipsi.post]);
    ymax = max(nanmean(dprime.difficult.ipsi.pre) ,nanmean([dprime.difficult.ipsi.post])) ;
    ext_sigline([y1,y2],criterion.difficult.ipsi.pvalue,[],ymax+ 0.5,'y', Color.Ipsi); hold on;
end


if dprime.difficult.contra.pvalue < 0.05
    ymax = max(nanmean(criterion.difficult.contra.pre) ,nanmean([criterion.difficult.contra.post]));
    ext_sigline([nanmean([dprime.difficult.contra.pre]),nanmean([dprime.difficult.contra.post])],dprime.difficult.contra.pvalue,[ ],(ymax+ 0.8 ),'x', Color.Contra ); hold on;
end
if criterion.difficult.contra.pvalue < 0.05
    y1 = nanmean(criterion.difficult.contra.pre) ;
    y2 = nanmean([criterion.difficult.contra.post]);
    ymax = max(nanmean(dprime.difficult.contra.pre) ,nanmean([dprime.difficult.contra.post])) ;
    ext_sigline([y1,y2],criterion.difficult.contra.pvalue,[],ymax+ 0.5,'y', Color.Contra ); hold on;
end
axis square
xlabel('sensitivity','fontsize',fs,'fontweight','b', 'Interpreter', 'none')
ylabel('criterion','fontsize',fs,'fontweight','b', 'Interpreter', 'none')
set(gca,'ylim',[-3 3],'xlim',[-1 4],'fontsize',fs)
text(-0.5,-2.8, 'more Contra & ipsi:NoGo, contra:Go', 'Color', 'k')
text(-0.5,2.8, 'less Contra & ipsi:Go, contra:NoGo', 'Color', 'k')


%% save the graphs
if SaveGraph
h= figure(1);
%savefig(h, [path_SaveFig ,filesep,'fig',filesep,  'DoubleStimuli_Crit_Dprime_' monkey,'_' experiment ,'_' folder '.fig'])
print(h,[path_SaveFig ,filesep,  ['png', ControlFolder]  ,filesep,  'DoubleSameStimuli_Crit_Dprime' monkey,'_' experiment,'_' folder '.png'], '-dpng')
set(h,'Renderer','Painters');
set(h,'PaperPositionMode','auto')
compl_filename =  [path_SaveFig , ['ai', ControlFolder]  ,filesep,  'DoubleSameStimuli_Crit_Dprime' monkey,'_' experiment ,'_', folder ,'.ai'] ;
print(h,'-depsc',compl_filename);
close all;
end



%%
figure('Position',[200 200 1200 900],'PaperPositionMode','auto'); % ,'PaperOrientation','landscape'
set(gcf,'Name',['Fixation - Miss',ControlFolder ]);
set(gcf,'Color',[1 1 1]);

    ha(1) = subplot(2,2,1); %1 % 5
        
    plot([1;2],[pMiss.pre; pMiss.post], 'o','color',[0 0 0] ,'MarkerSize',10,'markerfacecolor','k'); hold on;
    line(1:2,[pMiss.pre;pMiss.post],'Color',[0 0 0],'LineWidth', 2)
    for i = 1: length([pMiss.post])
        text(0.9,pMiss.pre(i),num2str(i),'fontsize',10)
        text(2.1,pMiss.post(i),num2str(i),'fontsize',10)
    end
    plot([1;2],[nanmean(pMiss.pre); nanmean(pMiss.post)], 'o','color',[1 0 0] ,'MarkerSize',10,'markerfacecolor','r'); hold on;
    line(1:2,[nanmean(pMiss.pre); nanmean(pMiss.post)],'Color',[1 0 0],'LineWidth', 2)
   
     if pMiss.pvalue < 0.05
        ext_sigline([1,2], pMiss.pvalue,[],0.8,'x' )
    end
    
    ylabel( 'Fixation as Miss','fontsize',14,'fontweight','b', 'Interpreter', 'none' );
    set(gca,'xlim',[0 3],'Xtick',1:2,'XTickLabel',{'pre' 'post'},'fontsize',20);
    set(gca,'ylim',[0 1])
    
    
      
     ha(2) = subplot(2,2,2); %1 % 5
        
    plot([1;2],[pCR.easy.pre;    pCR.easy.post], 'o','color',[0 0 0] ,'MarkerSize',10,'markerfacecolor','k'); hold on;
    line(1:2,[pCR.easy.pre;      pCR.easy.post],'Color',[0 0 0],'LineWidth', 2)
    for i = 1: length([pCR.easy.post])
        text(0.9,pCR.easy.pre(i),num2str(i),'fontsize',10)
        text(2.1,pCR.easy.post(i),num2str(i),'fontsize',10)
    end
    plot([1;2],[nanmean(pCR.easy.pre); nanmean(pCR.easy.post)], 'o','color',[1 0 0] ,'MarkerSize',10,'markerfacecolor','r'); hold on;
    line(1:2,[nanmean(pCR.easy.pre); nanmean(pCR.easy.post)],'Color',[1 0 0],'LineWidth', 2)
   
     if pCR.easy.pvalue < 0.05
        ext_sigline([1,2], pCR.easy.pvalue,[],0.2,'x' )
    end
    
    ylabel( 'Fixation as CR (easy)','fontsize',fs,'fontweight','b', 'Interpreter', 'none' );
    set(gca,'xlim',[0 3],'Xtick',1:2,'XTickLabel',{'pre' 'post'},'fontsize',20);
    set(gca,'ylim',[0 1])
    
    
     ha(3) = subplot(2,2,4); %1 % 5
        
    plot([1;2],[pCR.difficult.pre; pCR.difficult.post], 'o','color',[0 0 0] ,'MarkerSize',10,'markerfacecolor','k'); hold on;
    line(1:2,[pCR.difficult.pre;pCR.difficult.post],'Color',[0 0 0],'LineWidth', 2)
    for i = 1: length([pCR.difficult.post])
        text(0.9,pCR.difficult.pre(i),num2str(i),'fontsize',10)
        text(2.1,pCR.difficult.post(i),num2str(i),'fontsize',10)
    end
    plot([1;2],[nanmean(pCR.difficult.pre); nanmean(pCR.difficult.post)], 'o','color',[1 0 0] ,'MarkerSize',10,'markerfacecolor','r'); hold on;
    line(1:2,[nanmean(pCR.difficult.pre); nanmean(pCR.difficult.post)],'Color',[1 0 0],'LineWidth', 2)
   
     if pCR.difficult.pvalue < 0.05
        ext_sigline([1,2], pCR.difficult.pvalue,[],0.8,'x' )
    end
    
    ylabel( 'Fixation as CR (difficult)','fontsize',fs,'fontweight','b', 'Interpreter', 'none' );
    set(gca,'xlim',[0 3],'Xtick',1:2,'XTickLabel',{'pre' 'post'},'fontsize',20);
    set(gca,'ylim',[0 1])
    

if SaveGraph
    h = figure(1);
    print(h,[path_SaveFig ,filesep,  ['png', ControlFolder] ,filesep,  'DoubleSameStimuli_Fixation' monkey,'_' experiment ,'_' folder, '.png'], '-dpng')
    set(h,'Renderer','Painters');
    set(h,'PaperPositionMode','auto')
    compl_filename =  [path_SaveFig , ['ai', ControlFolder]  ,filesep,  'DoubleSameStimuli_Fixation.' monkey,'_' experiment ,'_', folder ,'.ai'] ;
    print(h,'-depsc',compl_filename);
    close all;
end

%% graph - From the paper Luo & Maunsell (2018)
figure('Position',[200 200 1200 900],'PaperPositionMode','auto'); % ,'PaperOrientation','landscape'
set(gcf,'Name',['Hitrate and FalseAlarmRate',ControlFolder ]);
set(gcf,'Color',[1 1 1]);

ha(1) = subplot(1,2,1);
for i = 1: length(pFA.easy.ipsi.pre)
%Fill the circle which is post
plot(pFA.easy.ipsi.post(i), pHit.ipsi.post(i), 'o' ,'color',Color.Ipsi, 'MarkerSize',MarkSize_GraphFAR_HR-5, 'MarkerFaceColor',Color.Ipsi,'LineWidth', 2); hold on;
plot(pFA.easy.contra.post(i), pHit.contra.post(i), 'o','color',Color.Contra , 'MarkerSize',MarkSize_GraphFAR_HR-5,'markerfacecolor',Color.Contra,'LineWidth', 2); hold on;

line([pFA.easy.ipsi.pre(i),pFA.easy.ipsi.post(i)], [pHit.ipsi.pre(i),pHit.ipsi.post(i)],'Color',[Color.Ipsi,0.3] , 'MarkerSize',MarkSize_GraphFAR_HR-5,'markerfacecolor',[1 1 1],'LineWidth', 2); hold on;
line([pFA.easy.contra.pre(i),pFA.easy.contra.post(i)], [pHit.contra.pre(i),pHit.contra.post(i)],'Color',[Color.Contra, 0.3] , 'MarkerSize',MarkSize_GraphFAR_HR-5,'markerfacecolor',[1 1 1],'LineWidth', 2); hold on;

end
%ipsi Post
plot([nanmean(pFA.easy.ipsi.pre),nanmean(pFA.easy.ipsi.post)], [nanmean(pHit.ipsi.pre),nanmean(pHit.ipsi.post)], 'o-','color',Color.Ipsi , 'MarkerSize',MarkSize_GraphFAR_HR,'markerfacecolor',[1 1 1],'LineWidth', 2); hold on;
plot([nanmean(pFA.easy.contra.pre),nanmean(pFA.easy.contra.post)], [nanmean(pHit.contra.pre),nanmean(pHit.contra.post)], 'o-','color',Color.Contra , 'MarkerSize',MarkSize_GraphFAR_HR,'markerfacecolor',[1 1 1],'LineWidth', 2); hold on;
%Fill the circle which is post
plot(nanmean(pFA.easy.ipsi.post), nanmean(pHit.ipsi.post), 'o','color',Color.Ipsi , 'MarkerSize',MarkSize_GraphFAR_HR-1,'markerfacecolor',Color.Ipsi,'LineWidth', 2); hold on;
plot(nanmean(pFA.easy.contra.post), nanmean(pHit.contra.post), 'o','color',Color.Contra , 'MarkerSize',MarkSize_GraphFAR_HR-1,'markerfacecolor',Color.Contra,'LineWidth', 2); hold on;

title('easy')

line([0 1],[1 0],'Color',[0 0 0],'LineStyle',':');
set(gca,'ylim',[0 1],'xlim',[0 1],'fontsize',fs)
xlabel( 'FA rate','fontsize',fs,'fontweight','b', 'Interpreter', 'none' );
ylabel( 'Hitrate','fontsize',fs,'fontweight','b', 'Interpreter', 'none' );
axis square
   if pFA.easy.ipsi.pvalue < 0.05
        y1 = nanmean(pFA.easy.ipsi.pre) ;
        y2 = nanmean(pFA.easy.ipsi.post);
        ymax = min(nanmean(pHit.ipsi.pre) ,nanmean(pHit.ipsi.post));
        ext_sigline([y1,y2],pFA.easy.ipsi.pvalue,[],ymax -0.2,'x', Color.Ipsi); hold on;
    end
    if pHit.ipsi.pvalue < 0.05
        y1 = nanmean(pHit.ipsi.pre) ;
        y2 = nanmean(pHit.ipsi.post);
        ymax = max(nanmean(pFA.easy.ipsi.pre) ,nanmean(pFA.easy.ipsi.post));
        ext_sigline([y1,y2],pHit.ipsi.pvalue,[],ymax +0.2,'y',Color.Ipsi); hold on;
    end   
    if pFA.easy.contra.pvalue < 0.05
        y1 = nanmean(pFA.easy.contra.pre) ;
        y2 = nanmean(pFA.easy.contra.post);
        ymax = min(nanmean(pHit.contra.pre) ,nanmean(pHit.contra.post));
        ext_sigline([y1,y2],pFA.easy.contra.pvalue,[],ymax -0.2,'x', Color.Contra); hold on;
    end
    if pHit.contra.pvalue < 0.05
        y1 = nanmean(pHit.contra.pre) ;
        y2 = nanmean(pHit.contra.post);
        ymax = max(nanmean(pFA.easy.contra.pre) ,nanmean(pFA.easy.contra.post));
        ext_sigline([y1,y2],pHit.contra.pvalue,[],ymax +0.2,'y', Color.Contra); hold on;
    end  




%%

ha(2) = subplot(1,2,2);

for i = 1: length(pFA.easy.ipsi.pre)
plot(pFA.difficult.ipsi.post(i), pHit.ipsi.post(i), 'o' ,'color',Color.Ipsi, 'MarkerSize',MarkSize_GraphFAR_HR-5, 'MarkerFaceColor',Color.Ipsi,'LineWidth', 2); hold on;
plot(pFA.difficult.contra.post(i), pHit.contra.post(i), 'o','color',Color.Contra , 'MarkerSize',MarkSize_GraphFAR_HR-5,'markerfacecolor',Color.Contra,'LineWidth', 2); hold on;

line([pFA.difficult.ipsi.pre(i),pFA.difficult.ipsi.post(i)], [pHit.ipsi.pre(i),pHit.ipsi.post(i)],'Color',[Color.Ipsi, 0.3] , 'MarkerSize',MarkSize_GraphFAR_HR-5,'markerfacecolor',[1 1 1],'LineWidth', 2); hold on;
line([pFA.difficult.contra.pre(i),pFA.difficult.contra.post(i)], [pHit.contra.pre(i),pHit.contra.post(i)],'Color',[Color.Contra, 0.3] , 'MarkerSize',MarkSize_GraphFAR_HR-5,'markerfacecolor',[1 1 1],'LineWidth', 2); hold on;
end
%ipsi Post
plot([nanmean(pFA.difficult.ipsi.pre),nanmean(pFA.difficult.ipsi.post)], [nanmean(pHit.ipsi.pre),nanmean(pHit.ipsi.post)], 'o-','color',Color.Ipsi , 'MarkerSize',MarkSize_GraphFAR_HR,'markerfacecolor',[1 1 1],'LineWidth', 2); hold on;
plot([nanmean(pFA.difficult.contra.pre),nanmean(pFA.difficult.contra.post)], [nanmean(pHit.contra.pre),nanmean(pHit.contra.post)], 'o-','color',Color.Contra , 'MarkerSize',MarkSize_GraphFAR_HR,'markerfacecolor',[1 1 1],'LineWidth', 2); hold on;
plot(nanmean(pFA.difficult.ipsi.post), nanmean(pHit.ipsi.post), 'o','color',Color.Ipsi , 'MarkerSize',MarkSize_GraphFAR_HR-1,'markerfacecolor',Color.Ipsi,'LineWidth', 2); hold on;
plot(nanmean(pFA.difficult.contra.post), nanmean(pHit.contra.post), 'o','color',Color.Contra , 'MarkerSize',MarkSize_GraphFAR_HR-1,'markerfacecolor',Color.Contra,'LineWidth', 2); hold on;
legend('ipsi', 'contra', 'Location', 'SouthEast')
title('difficult')

line([0 1],[1 0],'Color',[0 0 0],'LineStyle',':');
set(gca,'ylim',[0 1],'xlim',[0 1],'fontsize',fs)
xlabel( 'FA rate','fontsize',fs,'fontweight','b', 'Interpreter', 'none' );
ylabel( 'Hitrate','fontsize',fs,'fontweight','b', 'Interpreter', 'none' );
axis square

   if pFA.difficult.ipsi.pvalue < 0.05
        y1 = nanmean(pFA.difficult.ipsi.pre) ;
        y2 = nanmean(pFA.difficult.ipsi.post);
        ymax = min(nanmean(pHit.ipsi.pre) ,nanmean(pHit.ipsi.post));
        ext_sigline([y1,y2],pFA.difficult.ipsi.pvalue,[],ymax -0.3,'x', Color.Ipsi); hold on;
    end
    if pHit.ipsi.pvalue < 0.05
        y1 = nanmean(pHit.ipsi.pre) ;
        y2 = nanmean(pHit.ipsi.post);
        ymax = max(nanmean(pFA.difficult.ipsi.pre) ,nanmean(pFA.difficult.ipsi.post));
        ext_sigline([y1,y2],pHit.ipsi.pvalue,[],ymax +0.2,'y',Color.Ipsi); hold on;
    end   
    if pFA.difficult.contra.pvalue < 0.05
        y1 = nanmean(pFA.difficult.contra.pre) ;
        y2 = nanmean(pFA.difficult.contra.post);
        ymax = min(nanmean(pHit.contra.pre) ,nanmean(pHit.contra.post));
        ext_sigline([y1,y2],pFA.difficult.contra.pvalue,[],ymax -0.2,'x', Color.Contra); hold on;
    end
    if pHit.contra.pvalue < 0.05
        y1 = nanmean(pHit.contra.pre) ;
        y2 = nanmean(pHit.contra.post);
        ymax = max(nanmean(pFA.difficult.contra.pre) ,nanmean(pFA.difficult.contra.post));
        ext_sigline([y1,y2],pHit.contra.pvalue,[],ymax +0.2,'y',Color.Contra); hold on;
    end  


if SaveGraph
    h = figure(1);
    % savefig(h, [path_SaveFig ,'fig', filesep,  'TarDistStimuli_HR_FAR_' monkey,'_' experiment ,'_' folder ,'.fig'])
    print(h,[path_SaveFig ,filesep,  ['png', ControlFolder] ,filesep,  'DoubleSameStimuli_HR_FAR_AfterLuoetal' monkey,'_' experiment ,'_' folder, '.png'], '-dpng')
    set(h,'Renderer','Painters');
    set(h,'PaperPositionMode','auto')
    compl_filename =  [path_SaveFig , ['ai', ControlFolder]  ,filesep,  'DoubleSameStimuli_HR_FAR_AfterLuoetal.' monkey,'_' experiment ,'_', folder ,'.ai'] ;
    print(h,'-depsc',compl_filename);
    close all;
end
%% graph - Easy_Graph1HR_FAR_Graph2Sensitivity_Criterion
figure('Position',[200 200 1200 900],'PaperPositionMode','auto'); % ,'PaperOrientation','landscape'
set(gcf,'Name',['Easy_Graph1HR_FAR_Graph2Sensitivity_Criterion',ControlFolder ]);
set(gcf,'Color',[1 1 1]); 

ha(2) = subplot(1,2,1);
%ipsi Post

for i = 1: length(pFA.easy.ipsi.pre)
plot(pFA.easy.ipsi.post(i), pHit.ipsi.post(i), 'o' ,'color',Color.Ipsi, 'MarkerSize',MarkSize_GraphFAR_HR_PerSession, 'MarkerFaceColor',Color.Ipsi,'LineWidth', 2); hold on;
plot(pFA.easy.contra.post(i), pHit.contra.post(i), 'o','color',Color.Contra , 'MarkerSize',MarkSize_GraphFAR_HR_PerSession,'markerfacecolor',Color.Contra,'LineWidth', 2); hold on;

line([pFA.easy.ipsi.pre(i),pFA.easy.ipsi.post(i)], [pHit.ipsi.pre(i),pHit.ipsi.post(i)],'Color',[Color.Ipsi, 0.3] , 'MarkerSize',MarkSize_GraphFAR_HR_PerSession,'markerfacecolor',[1 1 1],'LineWidth', 2); hold on;
line([pFA.easy.contra.pre(i),pFA.easy.contra.post(i)], [pHit.contra.pre(i),pHit.contra.post(i)],'Color',[Color.Contra, 0.3] , 'MarkerSize',MarkSize_GraphFAR_HR_PerSession,'markerfacecolor',[1 1 1],'LineWidth', 2); hold on;

%Fill the circle which is post
end
plot([nanmean(pFA.easy.ipsi.pre),nanmean(pFA.easy.ipsi.post)], [nanmean(pHit.ipsi.pre),nanmean(pHit.ipsi.post)], 'o-','color',Color.Ipsi , 'MarkerSize',MarkSize_GraphFAR_HR,'markerfacecolor',[1 1 1],'LineWidth', LineWidthSize); hold on;
plot([nanmean(pFA.easy.contra.pre),nanmean(pFA.easy.contra.post)], [nanmean(pHit.contra.pre),nanmean(pHit.contra.post)], 'o-','color',Color.Contra , 'MarkerSize',MarkSize_GraphFAR_HR,'markerfacecolor',[1 1 1],'LineWidth', LineWidthSize); hold on;
plot(nanmean(pFA.easy.ipsi.post), nanmean(pHit.ipsi.post), 'o','color',Color.Ipsi , 'MarkerSize',MarkSize_GraphFAR_HR-1,'markerfacecolor',Color.Ipsi,'LineWidth', LineWidthSize); hold on;
plot(nanmean(pFA.easy.contra.post), nanmean(pHit.contra.post), 'o','color',Color.Contra , 'MarkerSize',MarkSize_GraphFAR_HR-1,'markerfacecolor',Color.Contra,'LineWidth', LineWidthSize); hold on;
legend('ipsi', 'contra', 'Location', 'South')

%title('easy')

line([0 1],[1 0],'Color',[0 0 0],'LineStyle',':');
set(gca,'ylim',[0 1],'xlim',[0 1],'fontsize',fs)
xlabel( 'FA rate','fontsize',fs,'fontweight','b', 'Interpreter', 'none' );
ylabel( 'Hitrate','fontsize',fs,'fontweight','b', 'Interpreter', 'none' );
axis square


   if pFA.easy.ipsi.pvalue < 0.05
        y1 = nanmean(pFA.easy.ipsi.pre) ;
        y2 = nanmean(pFA.easy.ipsi.post);
        ymax = min(nanmean(pHit.ipsi.pre) ,nanmean(pHit.ipsi.post));
        ext_sigline([y1,y2],pFA.easy.ipsi.pvalue,[],ymax -0.3,'x', Color.Ipsi); hold on;
    end
    if pHit.ipsi.pvalue < 0.05
        y1 = nanmean(pHit.ipsi.pre) ;
        y2 = nanmean(pHit.ipsi.post);
        ymax = max(nanmean(pFA.easy.ipsi.pre) ,nanmean(pFA.easy.ipsi.post));
        ext_sigline([y1,y2],pHit.ipsi.pvalue,[],ymax +0.2,'y',Color.Ipsi); hold on;
    end   
    if pFA.easy.contra.pvalue < 0.05
        y1 = nanmean(pFA.easy.contra.pre) ;
        y2 = nanmean(pFA.easy.contra.post);
        ymax = min(nanmean(pHit.contra.pre) ,nanmean(pHit.contra.post));
        ext_sigline([y1,y2],pFA.easy.contra.pvalue,[],ymax -0.2,'x', Color.Contra); hold on;
    end
    if pHit.contra.pvalue < 0.05
        y1 = nanmean(pHit.contra.pre) ;
        y2 = nanmean(pHit.contra.post);
        ymax = max(nanmean(pFA.easy.contra.pre) ,nanmean(pFA.easy.contra.post));
        ext_sigline([y1,y2],pHit.contra.pvalue,[],ymax +0.2,'y', Color.Contra); hold on;
    end 
    
    
ha(indPos) = subplot(1,2,2);
for i = 1: length(dprime.easy.contra.pre)
plot([dprime.easy.contra.pre(i),dprime.easy.contra.post(i)],[criterion.easy.contra.pre(i),criterion.easy.contra.post(i)], 'o','color',[Color.Contra, 0.4] , 'MarkerSize',MarkSize_CritDpr_small,'markerfacecolor',[1 1 1 ]); hold on;
plot([dprime.easy.ipsi.pre(i),dprime.easy.ipsi.post(i)], [-criterion.easy.ipsi.pre(i),-criterion.easy.ipsi.post(i)], 'o','color',[Color.Ipsi, 0.4] , 'MarkerSize',MarkSize_CritDpr_small,'markerfacecolor',[1 1 1 ]); hold on;% reverse direction of criterion for ipsi
plot(dprime.easy.contra.post(i),criterion.easy.contra.post(i), 'o','color',Color.Contra , 'MarkerSize',MarkSize_CritDpr_small,'markerfacecolor',Color.Contra); hold on;
plot(dprime.easy.ipsi.post(i),-criterion.easy.ipsi.post(i), 'o','color',Color.Ipsi , 'MarkerSize',MarkSize_CritDpr_small,'markerfacecolor',Color.Ipsi); hold on;% reverse direction of criterion for ipsi
end
plot([nanmean(dprime.easy.contra.pre),nanmean(dprime.easy.contra.post)],[nanmean(criterion.easy.contra.pre) ,nanmean(criterion.easy.contra.post)], 'o-','color',Color.Contra , 'MarkerSize',MarkSize_CritDpr,'markerfacecolor',[1 1 1 ], 'LineWidth', LineWidthSize); hold on;
plot([nanmean(dprime.easy.ipsi.pre),nanmean(dprime.easy.ipsi.post)],[ -nanmean(criterion.easy.ipsi.pre) ,-nanmean(criterion.easy.ipsi.post)], 'o-','color',Color.Ipsi , 'MarkerSize',MarkSize_CritDpr,'markerfacecolor',[1 1 1 ], 'LineWidth', LineWidthSize); hold on;% reverse direction of criterion for ipsi

plot(nanmean(dprime.easy.contra.post),nanmean(criterion.easy.contra.post), 'o','color',Color.Contra , 'MarkerSize',MarkSize_CritDpr,'markerfacecolor',Color.Contra); hold on;
plot(nanmean(dprime.easy.ipsi.post),-nanmean(criterion.easy.ipsi.post), 'o','color',Color.Ipsi , 'MarkerSize',MarkSize_CritDpr,'markerfacecolor',Color.Ipsi); hold on;% reverse direction of criterion for ipsi

if dprime.easy.ipsi.pvalue < 0.05
    ymax = max(nanmean(criterion.easy.ipsi.pre) ,nanmean([criterion.easy.ipsi.post]));
    ext_sigline([nanmean([dprime.easy.ipsi.pre]),nanmean([dprime.easy.ipsi.post])],dprime.easy.ipsi.pvalue,[ ],(ymax -1 ),'x', Color.Ipsi); hold on;
end
if criterion.easy.ipsi.pvalue < 0.05
    y1 = -1*nanmean(criterion.easy.ipsi.pre) ;
    y2 = -1*nanmean([criterion.easy.ipsi.post]);
    ymax = max(nanmean(dprime.easy.ipsi.pre) ,nanmean([dprime.easy.ipsi.post])) ;
    ext_sigline([y1,y2],criterion.easy.ipsi.pvalue,[],ymax+ 0.2,'y', Color.Ipsi); hold on;
end


if dprime.easy.contra.pvalue < 0.05
    ymax = max(nanmean(criterion.easy.contra.pre) ,nanmean([criterion.easy.contra.post]));
    ext_sigline([nanmean([dprime.easy.contra.pre]),nanmean([dprime.easy.contra.post])],dprime.easy.contra.pvalue,[ ],(ymax+ 0.8 ),'x', Color.Contra); hold on;
end
if criterion.easy.contra.pvalue < 0.05
    y1 = nanmean(criterion.easy.contra.pre) ;
    y2 = nanmean([criterion.easy.contra.post]);
    ymax = max(nanmean(dprime.easy.contra.pre) ,nanmean([dprime.easy.contra.post])) ;
    ext_sigline([y1,y2],criterion.easy.contra.pvalue,[],ymax+ 0.1,'y',Color.Contra); hold on;
end
axis square
xlabel('sensitivity','fontsize',fs,'fontweight','b', 'Interpreter', 'none')
ylabel('criterion','fontsize',fs,'fontweight','b', 'Interpreter', 'none')
set(gca,'ylim',[-3 3],'xlim',xlim_SDT_easy,'fontsize',fs)
        
text(xlim_SDT_easy(1),-2.8, 'more Contra (ipsi:NoGo, contra:Go)', 'Color', 'k','fontsize',14)
text(xlim_SDT_easy(1),2.8, 'less Contra (ipsi:Go, contra:NoGo)', 'Color', 'k','fontsize',14)

%% save the graphs
if SaveGraph
h= figure(1);
%savefig(h, [path_SaveFig ,filesep,'fig',filesep,  'DoubleStimuli_Crit_Dprime_' monkey,'_' experiment ,'_' folder '.fig'])
print(h,[path_SaveFig ,filesep,  ['png', ControlFolder]  ,filesep,  'DoubleSameStimuli_easy_HR_FAR_Dprime_Criterion' monkey,'_','.png'], '-dpng') % experiment,'_' folder 
set(h,'Renderer','Painters');
set(h,'PaperPositionMode','auto')
compl_filename =  [path_SaveFig , ['ai', ControlFolder]  ,filesep,  'DoubleSameStimuli_easy_HR_FAR_Dprime_Criterion' monkey,'_' ,'.ai'] ;% experiment ,'_', folder
print(h,'-depsc',compl_filename); 
close all;
end

%% GRAPH - TALK
figure('Position',[200 200 1200 900],'PaperPositionMode','auto'); % ,'PaperOrientation','landscape'
set(gcf,'Name',['Difficult_Graph1HR_FAR_Graph2Sensitivity_Criterion',ControlFolder ]);
set(gcf,'Color',[1 1 1]);
 
ha(2) = subplot(1,2,1);
for i = 1: length(pFA.easy.ipsi.pre)
plot(pFA.difficult.ipsi.post(i), pHit.ipsi.post(i), 'o' ,'color',Color.Ipsi, 'MarkerSize',MarkSize_CritDpr_small, 'MarkerFaceColor',Color.Ipsi,'LineWidth', 2); hold on;
plot(pFA.difficult.contra.post(i), pHit.contra.post(i), 'o','color',Color.Contra , 'MarkerSize',MarkSize_CritDpr_small,'markerfacecolor',Color.Contra,'LineWidth', 2); hold on;
line([pFA.difficult.ipsi.pre(i),pFA.difficult.ipsi.post(i)], [pHit.ipsi.pre(i),pHit.ipsi.post(i)],'Color',[Color.Ipsi, 0.3] , 'MarkerSize',MarkSize_CritDpr_small,'markerfacecolor',[1 1 1],'LineWidth', 2); hold on;
line([pFA.difficult.contra.pre(i),pFA.difficult.contra.post(i)], [pHit.contra.pre(i),pHit.contra.post(i)],'Color',[Color.Contra 0.3] , 'MarkerSize',MarkSize_CritDpr_small,'markerfacecolor',[1 1 1],'LineWidth', 2); hold on;

%Fill the circle which is post
end
plot([nanmean(pFA.difficult.ipsi.pre),nanmean(pFA.difficult.ipsi.post)], [nanmean(pHit.ipsi.pre),nanmean(pHit.ipsi.post)], 'o-','color',Color.Ipsi , 'MarkerSize',MarkSize_GraphFAR_HR,'markerfacecolor',[1 1 1],'LineWidth', LineWidthSize); hold on;
plot([nanmean(pFA.difficult.contra.pre),nanmean(pFA.difficult.contra.post)], [nanmean(pHit.contra.pre),nanmean(pHit.contra.post)], 'o-','color',Color.Contra , 'MarkerSize',MarkSize_GraphFAR_HR,'markerfacecolor',[1 1 1],'LineWidth', LineWidthSize); hold on;
plot(nanmean(pFA.difficult.ipsi.post), nanmean(pHit.ipsi.post), 'o','color',Color.Ipsi , 'MarkerSize',MarkSize_GraphFAR_HR-1,'markerfacecolor',Color.Ipsi,'LineWidth', 2); hold on;
plot(nanmean(pFA.difficult.contra.post), nanmean(pHit.contra.post), 'o','color',Color.Contra , 'MarkerSize',MarkSize_GraphFAR_HR-1,'markerfacecolor',Color.Contra,'LineWidth', 2); hold on;
legend('ipsi', 'contra', 'Location', 'NorthWest')

%title('difficult')

line([0 1],[1 0],'Color',[0 0 0],'LineStyle',':');
set(gca,'ylim',[0 1],'xlim',[0 1],'fontsize',fs)
xlabel( 'FA rate','fontsize',fs,'fontweight','b', 'Interpreter', 'none' );
ylabel( 'Hitrate','fontsize',fs,'fontweight','b', 'Interpreter', 'none' );
axis square


   if pFA.difficult.ipsi.pvalue < 0.05
        y1 = nanmean(pFA.difficult.ipsi.pre) ;
        y2 = nanmean(pFA.difficult.ipsi.post);
        if strcmp(monkey,'Cornelius' )
        ymax = min(nanmean(pHit.ipsi.pre) ,nanmean(pHit.ipsi.post))-0.3;
        else
        ymax = max(nanmean(pHit.ipsi.pre) ,nanmean(pHit.ipsi.post));

        end
        ext_sigline([y1,y2],pFA.difficult.ipsi.pvalue,[],ymax ,'x', Color.Ipsi); hold on;
    end
    if pHit.ipsi.pvalue < 0.05
        y1 = nanmean(pHit.ipsi.pre) ;
        y2 = nanmean(pHit.ipsi.post);
        ymax = max(nanmean(pFA.difficult.ipsi.pre) ,nanmean(pFA.difficult.ipsi.post));
        ext_sigline([y1,y2],pHit.ipsi.pvalue,[],ymax +0.1,'y',Color.Ipsi); hold on;
    end   
    if pFA.difficult.contra.pvalue < 0.05
        y1 = nanmean(pFA.difficult.contra.pre) ;
        y2 = nanmean(pFA.difficult.contra.post);
        ymax = min(nanmean(pHit.contra.pre) ,nanmean(pHit.contra.post));
        ext_sigline([y1,y2],pFA.difficult.contra.pvalue,[],ymax-0.1,'x', Color.Contra); hold on;
    end
    if pHit.contra.pvalue < 0.05
        y1 = nanmean(pHit.contra.pre) ;
        y2 = nanmean(pHit.contra.post);
        ymax = max(nanmean(pFA.difficult.contra.pre) ,nanmean(pFA.difficult.contra.post));
        ext_sigline([y1,y2],pHit.contra.pvalue,[],ymax +0.1,'y',Color.Contra); hold on;
    end 
    
    
ha(indPos) = subplot(1,2,2);
for i = 1: length(dprime.difficult.contra.pre)
plot([dprime.difficult.contra.pre(i),dprime.difficult.contra.post(i)],[criterion.difficult.contra.pre(i),criterion.difficult.contra.post(i)], 'o','color',[Color.Contra 0.4] , 'MarkerSize',MarkSize_CritDpr_small,'markerfacecolor',[1 1 1 ]); hold on;
plot([dprime.difficult.ipsi.pre(i),dprime.difficult.ipsi.post(i)], [-criterion.difficult.ipsi.pre(i),-criterion.difficult.ipsi.post(i)], 'o','color',[Color.Ipsi 0.4] , 'MarkerSize',MarkSize_CritDpr_small,'markerfacecolor',[1 1 1 ]); hold on;% reverse direction of criterion for ipsi
plot(dprime.difficult.contra.post(i),criterion.difficult.contra.post(i), 'o','color',Color.Contra , 'MarkerSize',MarkSize_CritDpr_small,'markerfacecolor',Color.Contra); hold on;
plot(dprime.difficult.ipsi.post(i),-criterion.difficult.ipsi.post(i), 'o','color',Color.Ipsi , 'MarkerSize',MarkSize_CritDpr_small,'markerfacecolor',Color.Ipsi); hold on;% reverse direction of criterion for ipsi
end
plot([nanmean(dprime.difficult.contra.pre),nanmean(dprime.difficult.contra.post)],[nanmean(criterion.difficult.contra.pre) ,nanmean(criterion.difficult.contra.post)], 'o-','color',Color.Contra , 'MarkerSize',MarkSize_CritDpr,'markerfacecolor',[1 1 1 ], 'LineWidth',LineWidthSize); hold on;
plot([nanmean(dprime.difficult.ipsi.pre),nanmean(dprime.difficult.ipsi.post)],[ -nanmean(criterion.difficult.ipsi.pre) ,-nanmean(criterion.difficult.ipsi.post)], 'o-','color',Color.Ipsi , 'MarkerSize',MarkSize_CritDpr,'markerfacecolor',[1 1 1 ], 'LineWidth',LineWidthSize); hold on;% reverse direction of criterion for ipsi

plot(nanmean(dprime.difficult.contra.post),nanmean(criterion.difficult.contra.post), 'o','color',Color.Contra , 'MarkerSize',MarkSize_CritDpr,'markerfacecolor',Color.Contra); hold on;
plot(nanmean(dprime.difficult.ipsi.post),-nanmean(criterion.difficult.ipsi.post), 'o','color',Color.Ipsi , 'MarkerSize',MarkSize_CritDpr,'markerfacecolor',Color.Ipsi); hold on;% reverse direction of criterion for ipsi

if dprime.difficult.ipsi.pvalue < 0.05
    ymax = max(nanmean(criterion.difficult.ipsi.pre) ,nanmean([criterion.difficult.ipsi.post]));
    ext_sigline([nanmean([dprime.difficult.ipsi.pre]),nanmean([dprime.difficult.ipsi.post])],dprime.difficult.ipsi.pvalue,[ ],(ymax -1 ),'x', Color.Ipsi); hold on;
end
if criterion.difficult.ipsi.pvalue < 0.05
    y1 = -1*nanmean(criterion.difficult.ipsi.pre) ;
    y2 = -1*nanmean([criterion.difficult.ipsi.post]);
    ymax = max(nanmean(dprime.difficult.ipsi.pre) ,nanmean([dprime.difficult.ipsi.post])) ;
    ext_sigline([y1,y2],criterion.difficult.ipsi.pvalue,[],ymax+ 1,'y', Color.Ipsi); hold on;
end


if dprime.difficult.contra.pvalue < 0.05
    ymax = max(nanmean(criterion.difficult.contra.pre) ,nanmean([criterion.difficult.contra.post]));
    ext_sigline([nanmean([dprime.difficult.contra.pre]),nanmean([dprime.difficult.contra.post])],dprime.difficult.contra.pvalue,[ ],(ymax+ 0.8 ),'x', Color.Contra); hold on;
end
if criterion.difficult.contra.pvalue < 0.05
    y1 = nanmean(criterion.difficult.contra.pre) ;
    y2 = nanmean([criterion.difficult.contra.post]);
    ymax = max(nanmean(dprime.difficult.contra.pre) ,nanmean([dprime.difficult.contra.post])) ;
    ext_sigline([y1,y2],criterion.difficult.contra.pvalue,[],ymax+ 0.7,'y', Color.Contra); hold on;
end
axis square
xlabel('sensitivity','fontsize',fs,'fontweight','b', 'Interpreter', 'none')
ylabel('criterion','fontsize',fs,'fontweight','b', 'Interpreter', 'none')
set(gca,'ylim',[-3 3],'xlim',xlim_SDT_diff ,'fontsize',fs)%[-1 4]
        
text(xlim_SDT_diff(1) +0.1,-2.8, 'more Contra (ipsi:NoGo, contra:Go)', 'Color', 'k','fontsize',14)
text(xlim_SDT_diff(1) +0.1,2.8, 'less Contra (ipsi:Go, contra:NoGo)', 'Color', 'k','fontsize',14)

%% save the graphs
if SaveGraph
h= figure(1);
%savefig(h, [path_SaveFig ,filesep,'fig',filesep,  'DoubleStimuli_Crit_Dprime_' monkey,'_' experiment ,'_' folder '.fig'])
print(h,[path_SaveFig ,filesep,  ['png', ControlFolder]  ,filesep,  'DoubleSameStimuli_difficult_HR_FAR_Dprime_Criterion' monkey, '.png'], '-dpng') %'_',experiment ,  folder
set(h,'Renderer','Painters');
set(h,'PaperPositionMode','auto')
compl_filename =  [path_SaveFig , ['ai', ControlFolder]  ,filesep,  'DoubleSameStimuli_difficult_HR_FAR_Dprime_Criterion' monkey,'.ai'] ; %'_' experiment ,'_', folder
print(h,'-depsc',compl_filename);
close all;
end

  