function SDT_singleStimuli(monkey, folder,dataset,folder_baseline,dataset_baseline,experiment,path_SaveFig,path_SaveData,path_getData, ResponseBias_DisplayColor)
%Input
close all;
SaveGraph = 1;
NonParametri = 0;
roundValue = 3;
SaveTable =1; 
%% Colors For Plotting
Color.Ipsi = [0 0 1 ];
Color.Contra = [1 0 1 ];
Color.Black = [0 0 0]; 
% Size of Writing
fs = 28; % font size
MarkSize_GraphFAR_HR = 12;
MarkSize = 20;

if strcmp(experiment, 'Microstimulation')
    xlim_SDT_easy = [1 5];
    xlim_SDT_diff = [0 4];

else
xlim_SDT_easy = [1 5];
xlim_SDT_diff = [0 4];
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
    SET.TimeLine        = {'stim_off', 'go_signal'};
   
elseif strcmp(experiment, 'Inactivation')
    Files = dir([path_getData, monkey ,'\behavior\' folder , filesep, '*.mat' ]); 
    load([path_getData , monkey ,'\behavior\' folder, filesep, Files(end).name ])
    disp( ['analyzed dataset:  ',  filesep,folder, filesep, Files(end).name ])

    SET.TimeLine        = {'stim_off', 'go_signal'};
    if isempty(folder_baseline)== 0
        Bline =load([path_getData , monkey ,'\behavior\' folder_baseline, filesep, dataset_baseline ]);
    end
    
elseif strcmp(experiment, 'Microstimulation')
    load([path_getData ,filesep, monkey , filesep, folder, filesep, dataset ])
    if strcmp(monkey, 'Curius')
        if  strcmp(folder, 'STIM_early_20161020_until_20161209')
            SET.TimeLine        = {'stim_off', 'minus_80ms'};
                        Stimulation = 'Early'; 

        elseif  strcmp(folder, 'STIM_late_20161020_until_20161209')
            SET.TimeLine        = {'stim_off', 'plus_80ms'};
                        Stimulation = 'Late'; 

        end
    elseif strcmp(monkey, 'Cornelius')
        if  strcmp(folder, 'STIM_early_20171126_until_20180123')
            SET.TimeLine        = {'stim_off', 'minus_250ms'};
                        Stimulation = 'Early'; 

        elseif  strcmp(folder, 'STIM_late_20180109_until_20180123')
            SET.TimeLine        = {'stim_off', 'minus_50ms'};
                        Stimulation = 'Late'; 

        end
    end
end

%% extract the data from the datafile
% single trials
%targets : successfull.. as saccaded (Hit) & unsuccessful as staying(Miss)
%contra = ipsi
% go_signal = pre
SET.maxSess =  length(data_struct_per_session.single_targ.ipsi.stim_off);
SET.positionExp         = {'ipsi', 'contra'};
SET.TimeLineExp         = {'pre', 'post'};
SET.Difficulty          = {'easy', 'difficult'};

Data = ExtractDataFromStruct(data_struct_per_session, SET);

%% ONLY FOR INACTIVATION: extract the Data for the Baseline-session
if strcmp(experiment, 'Inactivation') && isempty(folder_baseline) == 0
    maxSess =  length(Bline.data_struct_per_session.single_targ.ipsi.stim_off);
    positionExp     = {'ipsi', 'contra'};
    SET.TimeLineExp         = {'pre', 'post'};
    SET.Difficulty          = {'easy', 'difficult'};
    
    
    for indSess = 1:   maxSess
        for indPos = 1:length(positionExp)
            for indTime = 2
                if   strcmp(positionExp{indPos}, 'contra')
                    Data.Ctr_Session.(positionExp{indPos}).(SET.TimeLineExp{1})(indSess)               = Bline.data_struct_per_session.single_targ.(positionExp{indPos}).(SET.TimeLine{indTime}){1,indSess}.session;

                    Data.NumTotalTargets.(positionExp{indPos}).(SET.TimeLineExp{1})(indSess)               = Bline.data_struct_per_session.single_targ.(positionExp{indPos}).(SET.TimeLine{indTime}){1,indSess}.num_total_trials;
                    Data.Hit.(positionExp{indPos}).(SET.TimeLineExp{1})(indSess)                           = Bline.data_struct_per_session.single_targ.(positionExp{indPos}).(SET.TimeLine{indTime}){1,indSess}.num_choice_contra;
                    Data.Miss.(positionExp{indPos}).(SET.TimeLineExp{1})(indSess)                          = Bline.data_struct_per_session.single_targ.(positionExp{indPos}).(SET.TimeLine{indTime}){1,indSess}.num_choice_fixation;
                    Data.NumTotalDistractor.easy.(positionExp{indPos}).(SET.TimeLineExp{1})(indSess)       = Bline.data_struct_per_session.single_distr.(positionExp{indPos}).(SET.TimeLine{indTime}){2,indSess}.num_total_trials;
                    Data.NumTotalDistractor.difficult.(positionExp{indPos}).(SET.TimeLineExp{1})(indSess)  = Bline.data_struct_per_session.single_distr.(positionExp{indPos}).(SET.TimeLine{indTime}){1,indSess}.num_total_trials;
                    Data.CR.easy.(positionExp{indPos}).(SET.TimeLineExp{1})(indSess)                       = Bline.data_struct_per_session.single_distr.(positionExp{indPos}).(SET.TimeLine{indTime}){2,indSess}.num_choice_fixation;
                    Data.CR.difficult.(positionExp{indPos}).(SET.TimeLineExp{1})(indSess)                  = Bline.data_struct_per_session.single_distr.(positionExp{indPos}).(SET.TimeLine{indTime}){1,indSess}.num_choice_fixation;
                    Data.FA.easy.(positionExp{indPos}).(SET.TimeLineExp{1})(indSess)                       = Bline.data_struct_per_session.single_distr.(positionExp{indPos}).(SET.TimeLine{indTime}){2,indSess}.num_choice_contra;
                    Data.FA.difficult.(positionExp{indPos}).(SET.TimeLineExp{1})(indSess)                  = Bline.data_struct_per_session.single_distr.(positionExp{indPos}).(SET.TimeLine{indTime}){1,indSess}.num_choice_contra;
                    
                elseif     strcmp(positionExp{indPos}, 'ipsi')
                    Data.Ctr_Session.(positionExp{indPos}).(SET.TimeLineExp{1})(indSess)              = Bline.data_struct_per_session.single_targ.(positionExp{indPos}).(SET.TimeLine{indTime}){1,indSess}.session;
                    Data.NumTotalTargets.(positionExp{indPos}).(SET.TimeLineExp{1})(indSess)               =  Bline.data_struct_per_session.single_targ.(positionExp{indPos}).(SET.TimeLine{indTime}){1,indSess}.num_total_trials;
                    Data.Hit.(positionExp{indPos}).(SET.TimeLineExp{1})(indSess)                           =  Bline.data_struct_per_session.single_targ.(positionExp{indPos}).(SET.TimeLine{indTime}){1,indSess}.num_choice_ipsi;
                    Data.Miss.(positionExp{indPos}).(SET.TimeLineExp{1})(indSess)                          =  Bline.data_struct_per_session.single_targ.(positionExp{indPos}).(SET.TimeLine{indTime}){1,indSess}.num_choice_fixation;
                    Data.NumTotalDistractor.easy.(positionExp{indPos}).(SET.TimeLineExp{1})(indSess)       =  Bline.data_struct_per_session.single_distr.(positionExp{indPos}).(SET.TimeLine{indTime}){2,indSess}.num_total_trials;
                    Data.NumTotalDistractor.difficult.(positionExp{indPos}).(SET.TimeLineExp{1})(indSess)  =  Bline.data_struct_per_session.single_distr.(positionExp{indPos}).(SET.TimeLine{indTime}){1,indSess}.num_total_trials;
                    Data.CR.easy.(positionExp{indPos}).(SET.TimeLineExp{1})(indSess)                       =  Bline.data_struct_per_session.single_distr.(positionExp{indPos}).(SET.TimeLine{indTime}){2,indSess}.num_choice_fixation;
                    Data.CR.difficult.(positionExp{indPos}).(SET.TimeLineExp{1})(indSess)                  =  Bline.data_struct_per_session.single_distr.(positionExp{indPos}).(SET.TimeLine{indTime}){1,indSess}.num_choice_fixation;
                    Data.FA.easy.(positionExp{indPos}).(SET.TimeLineExp{1})(indSess)                       =  Bline.data_struct_per_session.single_distr.(positionExp{indPos}).(SET.TimeLine{indTime}){2,indSess}.num_choice_ipsi;
                    Data.FA.difficult.(positionExp{indPos}).(SET.TimeLineExp{1})(indSess)                  =  Bline.data_struct_per_session.single_distr.(positionExp{indPos}).(SET.TimeLine{indTime}){1,indSess}.num_choice_ipsi;
                    
                end
            end
        end
    end
    
    
end
%%  avoid 0 or Inf probabilities
%use the approach regardless wheather or not extreme values are obtained
for indPos = 1:length(SET.positionExp ) %contra vs ipsi selection
    for indTime = 1:length(SET.TimeLineExp) %pre vs post
        Data.Hit.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime})            = Data.Hit.(SET.positionExp{indPos}).(SET.TimeLineExp{indTime}) + 0.5 ;
        Data.Miss.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime})           = Data.Miss.(SET.positionExp{indPos}).(SET.TimeLineExp{indTime}) + 0.5 ;
        Data.FA.easy.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime})        = Data.FA.easy.(SET.positionExp{indPos}).(SET.TimeLineExp{indTime}) + 0.5 ;
        Data.FA.difficult.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime})   = Data.FA.difficult.(SET.positionExp{indPos}).(SET.TimeLineExp{indTime}) + 0.5;
        Data.CR.easy.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime})        = Data.CR.easy.(SET.positionExp{indPos}).(SET.TimeLineExp{indTime})  + 0.5;
        Data.CR.difficult.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime})   = Data.CR.difficult.(SET.positionExp{indPos}).(SET.TimeLineExp{indTime})  + 0.5;
    end
end

%% calculate d'prime & criterion
for indPos = 1:length(SET.positionExp ) %contra vs ipsi selection
    for indTime = 1:length(SET.TimeLine) %pre vs post
        %target is presented
        Data.pHit.(SET.positionExp{indPos}).(SET.TimeLineExp{indTime})=    Data.Hit.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime}) ./ ...
            sum(Data.Hit.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime})  + Data.Miss.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime}),1);
        %distractor is presented
        Data.pFA.easy.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime}) = Data.FA.easy.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime}) ./ ...
            sum(Data.FA.easy.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime})+  Data.CR.easy.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime}),1);
        
        Data.pFA.difficult.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime}) = Data.FA.difficult.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime}) ./ ...
            sum(Data.FA.difficult.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime})+  Data.CR.difficult.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime}),1);
        %distractor is presented
        Data.pMiss.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime}) = Data.Miss.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime}) ./ ...
            sum(Data.Hit.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime})  + Data.Miss.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime}),1);
        Data.pCR.easy.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime}) = Data.CR.easy.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime}) ./ ...
            sum(Data.FA.easy.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime})  + Data.CR.easy.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime}),1);
        Data.pCR.difficult.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime}) = Data.CR.difficult.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime}) ./ ...
            sum(Data.FA.difficult.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime})  + Data.CR.difficult.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime}),1);        
        
         Data.pFix.easy.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime}) = (Data.CR.easy.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime}) +  Data.Miss.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime}) ) ./ ...
            sum(Data.FA.easy.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime})  + Data.CR.easy.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime}) +Data.Hit.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime})  + Data.Miss.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime}),1);
        Data.pFix.difficult.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime}) = (Data.CR.difficult.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime}) + Data.Miss.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime}) )./ ...
            sum(Data.FA.difficult.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime})  + Data.CR.difficult.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime}) +Data.Hit.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime})  + Data.Miss.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime}),1);         
        
        
        
        [Data.dprime.easy.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime}),Data.beta.easy.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime}),Data.criterion.easy.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime})] = ...
            testsim_dprime( Data.pHit.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime}), Data.pFA.easy.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime}));
        [Data.dprime.difficult.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime}),Data.beta.difficult.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime}),Data.criterion.difficult.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime})] = ...
            testsim_dprime( Data.pHit.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime}), Data.pFA.difficult.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime}));
        
        Data.z_dprime.easy.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime}) = zscore(Data.dprime.easy.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime}));
        Data.z_dprime.difficult.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime}) = zscore(Data.dprime.difficult.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime}));
        
        Data.z_criterion.easy.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime}) = zscore(Data.criterion.easy.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime}));
        Data.z_criterion.difficult.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime}) = zscore(Data.criterion.difficult.(SET.positionExp {indPos}).(SET.TimeLineExp{indTime}));
        
        
    end
end



%% %%%%%%%%  SAVE THE DATASET %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T = table([nan;nan], [nan;nan], [nan;nan]) ; 
Table = [];
        all_DV = {'criterion', 'dprime'};
        for indVar = 1: length(all_DV)
            for indPos = 1: length(SET.positionExp)
                for indTime = 1: length(SET.TimeLineExp)
                    for  indDiff = 1: length(SET.Difficulty)
                        Values = []; Hemifield = []; TimeLineExperiment = []; StimulusTyp = []; DisDifficulty = []; Monkey = [];  DependentVariable = []; Date = []; 
                        Experiment = [];
                        
                        Values = Data.(all_DV{indVar}).(SET.Difficulty {indDiff}).(SET.positionExp {indPos}).(SET.TimeLineExp {indTime})';
                        Hemifield   = repmat(SET.positionExp(indPos), length(Values),1);
                        TimeLineExperiment    =  repmat(SET.TimeLineExp(indTime), length(Values),1);
                        StimulusTyp = repmat({'SingleStim'}, length(Values),1);
                        DisDifficulty = repmat(SET.Difficulty (indDiff), length(Values),1);
                        Monkey = repmat(monkey, length(Values),1);

                        if strcmp(experiment, 'Inactivation')&&strcmp((SET.TimeLineExp {indTime}), 'post')
                            Date = Data.Session.(SET.positionExp{indPos}).(SET.TimeLineExp{indTime})';
                            Experiment   = repmat({'Ina'}, length(Date),1);
                            
                        elseif strcmp(experiment, 'Inactivation') && strcmp((SET.TimeLineExp {indTime}), 'pre')
                            TimeLineExperiment = repmat({'post'}, length(Values),1);
                            Date = Data.Ctr_Session.(SET.positionExp {indPos}).(SET.TimeLineExp {indTime})';
                            Experiment   = repmat({'Ctr'}, length(Date),1);
                       
                        elseif strcmp(experiment, 'Microstimulation') 
                            Date        = Data.Session.(SET.positionExp {indPos}).(SET.TimeLineExp {indTime})';
                            if strcmp(Stimulation, 'Late')
                            Experiment    = repmat({'Late'}, length(Date),1);
                            elseif strcmp(Stimulation, 'Early')
                            Experiment    = repmat({'Early'}, length(Date),1);  
                            end
                            if  strcmp((SET.TimeLineExp {indTime}), 'pre')
                            TimeLineExperiment = repmat({'off'}, length(Values),1);
                            else
                            TimeLineExperiment = repmat(SET.TimeLine(2), length(Values),1);  
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
filename =[path_SaveFig,['Table', ControlFolder], filesep, monkey, '_SDTvar_SingleStim' ] ;  
writetable(Table, filename , 'Delimiter', ' ')
disp(['SAVED   ', filename])

end



%% Statistic - 
  TablePwerte = []; 
if strcmp(experiment, 'Inactivation') && isempty(folder_baseline) == 0 && NonParametri == 0
  %  independent t-test for ctr vs. inactivation session
    for indPos = 1:length(SET.positionExp ) %contra vs ipsi selection
            [h,p,ci,stats]  = ttest2(Data.dprime.easy.(SET.positionExp {indPos}).pre, Data.dprime.easy.(SET.positionExp {indPos}).post);
            Data.dprime.easy.(SET.positionExp {indPos}).pvalue = p;
            Data.dprime.easy.(SET.positionExp {indPos}).tstat = stats.tstat;           
            disp(['dprime ', 'easy ' ,(SET.positionExp {indPos}) ,  '  medPre =', num2str(round(median(Data.dprime.easy.(SET.positionExp {indPos}).pre),2)),' ' ,...
                'medPost =', num2str(round(median(Data.dprime.easy.(SET.positionExp {indPos}).post),2)),' p=' ,num2str(round(p,2)),' t=', num2str(round(Data.dprime.easy.(SET.positionExp {indPos}).tstat,2))] )
        TabPwerte = [ ];
        TabPwerte = table({monkey}, {'SingleStim'}, {'dprime'},{ 'easy'} , (SET.positionExp(indPos)) ,round(p,roundValue),{'t'} , round(stats.tstat,3), {'indepTtest'} );
        TablePwerte = [ TablePwerte; TabPwerte ];             
%             T.Var3(indPos) = p; 
%             T.Var2 = table(SET.positionExp {indPos}); 
%             table(SET.positionExp {indPos}, 'easy', Data.dprime.easy.(SET.positionExp {indPos}).pvalue);
%              %STATS=mwwtest(Data.dprime.easy.(SET.positionExp {indPos}).pre,Data.dprime.easy.(SET.positionExp {indPos}).post)
            
            [h,p,ci,stats] = ttest2(Data.criterion.easy.(SET.positionExp {indPos}).pre, Data.criterion.easy.(SET.positionExp {indPos}).post);
            Data.criterion.easy.(SET.positionExp {indPos}).pvalue = p;
            Data.criterion.easy.(SET.positionExp {indPos}).tstat = stats.tstat;
            disp(['criterion ', 'easy ' ,(SET.positionExp {indPos}) ,  '  medPre =', num2str(round(median(Data.criterion.easy.(SET.positionExp {indPos}).pre),2)),' ' ,...
                'medPost =', num2str(round(median(Data.criterion.easy.(SET.positionExp {indPos}).post),2)),' p=' ,num2str(round(p,2)),...
                ' t=', num2str(round(Data.criterion.easy.(SET.positionExp {indPos}).tstat,2))] )
        TabPwerte = [ ];
        TabPwerte = table({monkey}, {'SingleStim'}, {'criterion'},{ 'easy'} , (SET.positionExp(indPos)) ,round(p,roundValue),{'t'} , round(stats.tstat,3), {'indepTtest'} );
        TablePwerte = [ TablePwerte; TabPwerte ]; 
        
             [h,p,ci,stats]= ttest2(Data.dprime.difficult.(SET.positionExp {indPos}).pre, Data.dprime.difficult.(SET.positionExp {indPos}).post);
            Data.dprime.difficult.(SET.positionExp {indPos}).pvalue = p;
            Data.dprime.difficult.(SET.positionExp {indPos}).tstat = stats.tstat;
            disp(['dprime ', 'difficult ' , (SET.positionExp {indPos}) ,  '  medPre =', num2str(round(median(Data.dprime.difficult.(SET.positionExp {indPos}).pre),2)),' ' ,...
                'medPost =', num2str(round(median(Data.dprime.difficult.(SET.positionExp {indPos}).post),2)),' p=' ,num2str(round(p,2)) ' t=', num2str(round(Data.dprime.difficult.(SET.positionExp {indPos}).tstat,2))])
        TabPwerte = [ ];
        TabPwerte = table({monkey}, {'SingleStim'}, {'dprime'},{ 'difficult'} , (SET.positionExp(indPos)) ,round(p,roundValue),{'t'} , round(stats.tstat,3), {'indepTtest'} );
        TablePwerte = [ TablePwerte; TabPwerte ]; 

             [h,p,ci,stats] = ttest2(Data.criterion.difficult.(SET.positionExp {indPos}).pre, Data.criterion.difficult.(SET.positionExp {indPos}).post);
            Data.criterion.difficult.(SET.positionExp {indPos}).pvalue = p;
            Data.criterion.difficult.(SET.positionExp {indPos}).tstat = stats.tstat;
            disp(['criterion ', 'difficult '  , (SET.positionExp {indPos}) ,  '  medPre =', num2str(round(median(Data.criterion.difficult.(SET.positionExp {indPos}).pre),2)),' ' ,...
                'medPost =', num2str(round(median(Data.criterion.difficult.(SET.positionExp {indPos}).post),2)),' p='  ,num2str(round(p,2)) ' t=', num2str(round(Data.criterion.difficult.(SET.positionExp {indPos}).tstat,2))] )
         TabPwerte = [ ];
        TabPwerte = table({monkey}, {'SingleStim'}, {'criterion'},{ 'difficult'} , (SET.positionExp(indPos)) ,round(p,roundValue),{'t'} , round(stats.tstat,3), {'indepTtest'} );
        TablePwerte = [ TablePwerte; TabPwerte ]; 
                   

    end
    for indPos = 1:length(SET.positionExp ) %contra vs ipsi selection
            [h,p,ci,stat] = ttest2(Data.pHit.(SET.positionExp {indPos}).pre, Data.pHit.(SET.positionExp {indPos}).post,'tail','both','alpha',0.05);
            Data.pHit.(SET.positionExp {indPos}).pvalue = p;
            Data.pHit.(SET.positionExp{indPos}).tstat = stat.tstat;
            Data.pHit.(SET.positionExp {indPos}).r = norminv(p)/sqrt(length(Data.pHit.(SET.positionExp {indPos}).post) + length(Data.pHit.(SET.positionExp {indPos}).pre)); 
            median_Pre  = median(Data.pHit.(SET.positionExp {indPos}).pre) ; 
            median_Post = median(Data.pHit.(SET.positionExp {indPos}).post) ; 
            disp(['pHit ' , (SET.positionExp {indPos}) ,'  medPre =',num2str(median_Pre) ,'  medPost =',  num2str(median_Post),  ' p=' ,num2str(round(p,2)),' t=', num2str(round(Data.pHit.(SET.positionExp {indPos}).tstat,2))] )

            
           [h,p,ci,stat] = ttest2(Data.pFA.easy.(SET.positionExp {indPos}).pre, Data.pFA.easy.(SET.positionExp {indPos}).post,'tail','both','alpha',0.05);
            Data.pFA.easy.(SET.positionExp {indPos}).pvalue = p;
            Data.pFA.easy.(SET.positionExp{indPos}).tstat = stat.tstat;
            Data.pFA.easy.(SET.positionExp {indPos}).r = norminv(p)/sqrt(length(Data.pFA.easy.(SET.positionExp {indPos}).post) + length(Data.pFA.easy.(SET.positionExp {indPos}).pre)); 
            median_Pre  = median(Data.pFA.easy.(SET.positionExp {indPos}).pre) ; 
            median_Post = median(Data.pFA.easy.(SET.positionExp {indPos}).post) ;
            disp(['pFA ', 'easy ' , (SET.positionExp {indPos}) ,'  medPre =',num2str(median_Pre) ,'  medPost =',  num2str(median_Post),  ' p=' ,num2str(round(p,2)),'_t=', num2str(round(Data.pFA.easy.(SET.positionExp {indPos}).tstat,2))] )

            
            [h,p,ci,stat] = ttest2(Data.pFA.difficult.(SET.positionExp {indPos}).pre, Data.pFA.difficult.(SET.positionExp {indPos}).post,'tail','both','alpha',0.05);
            Data.pFA.difficult.(SET.positionExp {indPos}).pvalue = p;
            Data.pFA.difficult.(SET.positionExp{indPos}).tstat = stat.tstat;
            Data.pFA.difficult.(SET.positionExp {indPos}).r = norminv(p)/sqrt(length(Data.pFA.difficult.(SET.positionExp {indPos}).post) + length(Data.pFA.difficult.(SET.positionExp {indPos}).pre)); 
            median_Pre  = median(Data.pFA.difficult.(SET.positionExp {indPos}).pre) ; 
            median_Post = median(Data.pFA.difficult.(SET.positionExp {indPos}).post) ;
            disp(['pFA ', 'difficult ' , (SET.positionExp {indPos}) ,'  medPre =',num2str(median_Pre) ,'  medPost =',  num2str(median_Post),  ' p=' ,num2str(round(p,2)),'_t=', num2str(round(Data.pFA.difficult.(SET.positionExp {indPos}).tstat,2))] )

            [h,p,ci,stat] = ttest2(Data.pMiss.(SET.positionExp {indPos}).pre, Data.pMiss.(SET.positionExp {indPos}).post,'tail','both','alpha',0.05);
            Data.pMiss.(SET.positionExp {indPos}).pvalue = p;
            [h,p,ci,stat] = ttest2(Data.pCR.easy.(SET.positionExp {indPos}).pre, Data.pCR.easy.(SET.positionExp {indPos}).post,'tail','both','alpha',0.05);
            
            Data.pCR.easy.(SET.positionExp {indPos}).pvalue = p;
           [h,p,ci,stat] = ttest2(Data.pCR.difficult.(SET.positionExp {indPos}).pre, Data.pCR.difficult.(SET.positionExp {indPos}).post);
            Data.pCR.difficult.(SET.positionExp {indPos}).pvalue = p;      
            
            [h,p,ci,stat] = ttest2(Data.pFix.easy.(SET.positionExp {indPos}).pre, Data.pFix.easy.(SET.positionExp {indPos}).post,'tail','both','alpha',0.05);
            Data.pFix.easy.(SET.positionExp {indPos}).pvalue = p;
            Data.pFix.easy.(SET.positionExp {indPos}).r = norminv(p)/sqrt(length(Data.pFix.easy.(SET.positionExp {indPos}).post) + length(Data.pFix.easy.(SET.positionExp {indPos}).pre)); 
            median_Pre  = median(Data.pFix.easy.(SET.positionExp {indPos}).pre) ; 
            median_Post = median(Data.pFix.easy.(SET.positionExp {indPos}).post) ;
            disp(['pFix ', 'easy ' , (SET.positionExp {indPos}) ,'  medPre =',num2str(median_Pre) ,'  medPost =',  num2str(median_Post),  ' p=' ,num2str(round(p,2)),'_r=', num2str(round(Data.pFix.easy.(SET.positionExp {indPos}).r,2))] )

            [h,p,ci,stat] = ttest2(Data.pFix.difficult.(SET.positionExp {indPos}).pre, Data.pFix.difficult.(SET.positionExp {indPos}).post,'tail','both','alpha',0.05);
            Data.pFix.difficult.(SET.positionExp {indPos}).pvalue = p;  
            Data.pFix.difficult.(SET.positionExp {indPos}).r = norminv(p)/sqrt(length(Data.pFix.difficult.(SET.positionExp {indPos}).post) + length(Data.pFix.difficult.(SET.positionExp {indPos}).pre)); 
            median_Pre  = median(Data.pFix.difficult.(SET.positionExp {indPos}).pre) ; 
            median_Post = median(Data.pFix.difficult.(SET.positionExp {indPos}).post) ;
            disp(['pFix ', 'difficult ' , (SET.positionExp {indPos}) ,'  medPre =',num2str(median_Pre) ,'  medPost =',  num2str(median_Post),  ' p=' ,num2str(round(p,2)),'_r=', num2str(round(Data.pFix.difficult.(SET.positionExp {indPos}).r,2))] )

            
    end
else
    for indPos = 1:length(SET.positionExp ) %contra vs ipsi selection
        [p,h,stat] = ranksum(Data.dprime.easy.(SET.positionExp {indPos}).pre, Data.dprime.easy.(SET.positionExp {indPos}).post);
        Data.dprime.easy.(SET.positionExp {indPos}).pvalue = p;
        disp(['dprime ', 'easy ' ,(SET.positionExp {indPos}) ,  '  medPre =', num2str(round(median(Data.dprime.easy.(SET.positionExp {indPos}).pre),2)),' ' ,...
            'medPost =', num2str(round(median(Data.dprime.easy.(SET.positionExp {indPos}).post),2)),' p=' ,num2str(round(p,2))] )
        
        TabPwerte = [ ];
        TabPwerte = table({monkey }, {'SingleStim'}, {'dprime'},{ 'easy'} , (SET.positionExp(indPos)) ,round(p,roundValue),{'U'} , stat.ranksum, {'ranksum '} );
        TablePwerte = [ TablePwerte; TabPwerte ];
        
        
        [p,h,stat] = ranksum(Data.criterion.easy.(SET.positionExp {indPos}).pre, Data.criterion.easy.(SET.positionExp {indPos}).post);
        Data.criterion.easy.(SET.positionExp {indPos}).pvalue = p;
        TabPwerte = [ ];
        TabPwerte = table({monkey }, {'SingleStim'}, {'criterion'},{ 'easy'} , (SET.positionExp(indPos)) ,round(p,roundValue),{'U'} , stat.ranksum, {'ranksum '} );
        TablePwerte = [ TablePwerte; TabPwerte ];
        
        [p,h,stat] = ranksum(Data.dprime.difficult.(SET.positionExp {indPos}).pre, Data.dprime.difficult.(SET.positionExp {indPos}).post);
        Data.dprime.difficult.(SET.positionExp {indPos}).pvalue = p;
        disp(['dprime ', 'difficult ' , (SET.positionExp {indPos}) ,  '  medPre =', num2str(round(median(Data.dprime.difficult.(SET.positionExp {indPos}).pre),2)),' ' ,...
            'medPost =', num2str(round(median(Data.dprime.difficult.(SET.positionExp {indPos}).post),2)),' p=' ,num2str(round(p,2)) ])
        TabPwerte = [ ];
        TabPwerte = table({monkey }, {'SingleStim'}, {'dprime'},{ 'difficult'} , (SET.positionExp(indPos)) ,round(p,roundValue),{'U'} , stat.ranksum, {'ranksum '} );
        TablePwerte = [ TablePwerte; TabPwerte ];
        
        [p,h,stat] = ranksum(Data.criterion.difficult.(SET.positionExp {indPos}).pre, Data.criterion.difficult.(SET.positionExp {indPos}).post);
        Data.criterion.difficult.(SET.positionExp {indPos}).pvalue = p;
       TabPwerte = [ ];
        TabPwerte = table({monkey }, {'SingleStim'}, {'criterion'},{ 'difficult'} , (SET.positionExp(indPos)) ,round(p,roundValue),{'U'} , stat.ranksum, {'ranksum '} );
        TablePwerte = [ TablePwerte; TabPwerte ];
 
    end
    for indPos = 1:length(SET.positionExp ) %contra vs ipsi selection
        [H,p,CI,STATS] = ttest(Data.pHit.(SET.positionExp {indPos}).pre, Data.pHit.(SET.positionExp {indPos}).post,'tail','both','alpha',0.05);
        Data.pHit.(SET.positionExp {indPos}).pvalue = p;
        TabPwerte = [ ];
        TabPwerte = table({monkey}, {'SingleStim'}, {'pHit'},{ 'nan'} , (SET.positionExp(indPos)) ,round(p,roundValue),{'t'} , round(STATS.tstat,3), {'depTtest'} );
        TablePwerte = [ TablePwerte; TabPwerte ];
        [H,p,CI,STATS] =ttest(Data.pFA.easy.(SET.positionExp {indPos}).pre, Data.pFA.easy.(SET.positionExp {indPos}).post,'tail','both','alpha',0.05);
        Data.pFA.easy.(SET.positionExp {indPos}).pvalue = p;
        TabPwerte = [ ];
        TabPwerte = table({monkey}, {'SingleStim'}, {'pFA'},{ 'easy'} , (SET.positionExp(indPos)) ,round(p,roundValue),{'t'} , round(STATS.tstat,3), {'depTtest'} );
        TablePwerte = [ TablePwerte; TabPwerte ];
        [H,p,CI,STATS] = ttest(Data.pFA.difficult.(SET.positionExp {indPos}).pre, Data.pFA.difficult.(SET.positionExp {indPos}).post,'tail','both','alpha',0.05);
        Data.pFA.difficult.(SET.positionExp {indPos}).pvalue = p;
        
        TabPwerte = [ ];
        TabPwerte = table({monkey}, {'SingleStim'}, {'pFA'},{ 'difficult'} , (SET.positionExp(indPos)) ,round(p,roundValue),{'t'} , round(STATS.tstat,3), {'depTtest'} );
        TablePwerte = [ TablePwerte; TabPwerte ];
        
        
        [h,p,stat] = ttest(Data.pMiss.(SET.positionExp {indPos}).pre, Data.pMiss.(SET.positionExp {indPos}).post);
        Data.pMiss.(SET.positionExp {indPos}).pvalue = p;
        [h,p,stat] = ttest(Data.pCR.easy.(SET.positionExp {indPos}).pre, Data.pCR.easy.(SET.positionExp {indPos}).post);
        Data.pCR.easy.(SET.positionExp {indPos}).pvalue = p;
        [h,p,stat] = ttest(Data.pCR.difficult.(SET.positionExp {indPos}).pre, Data.pCR.difficult.(SET.positionExp {indPos}).post);
        Data.pCR.difficult.(SET.positionExp {indPos}).pvalue = p;
        
        
        [h,p,stat] = ttest(Data.pFix.easy.(SET.positionExp {indPos}).pre, Data.pFix.easy.(SET.positionExp {indPos}).post);
        Data.pFix.easy.(SET.positionExp {indPos}).pvalue = p;
        [h,p,stat] = ttest(Data.pFix.difficult.(SET.positionExp {indPos}).pre, Data.pFix.difficult.(SET.positionExp {indPos}).post);
        Data.pFix.difficult.(SET.positionExp {indPos}).pvalue = p;

 
    end
end
if strcmp(experiment, 'Inactivation') && SaveTable == 1;
    filename =[path_SaveFig,['Table', ControlFolder], filesep, monkey, '_SDTvarPwerte_SingleStim.xlsx' ] ;
    writetable(TablePwerte,filename,'Sheet',1,  'Range' ,'A1' )
    disp(['SAVED   ', filename])
    
elseif strcmp(experiment, 'Microstimulation')
    
   % filename =[addtoDropbox,filesep, monkey,filesep,   'SDTvariables_Microstimulation' , filesep, experiment, filesep, 'Stat' ,filesep, monkey,'_', Stimulation, '_SDTvarPwerte_SingleStim.xlsx' ] ;
    writetable(TablePwerte,filename,'Sheet',1,  'Range' ,'A1' )
    
    disp(['SAVED  ' , path_SaveFig ,filesep,  'Stat' ,filesep, monkey,'_' , Stimulation,  '_SDTvarPwerte_SingleStim' ])
    
end

%% Determine the if the dprime differs against the H0?
for indPos = 1:length(SET.positionExp ) %contra vs ipsi selection
    for indTime = 1:length(SET.TimeLine) %pre vs post
        [p,h,stat] = signrank(Data.dprime.easy.(SET.positionExp {indPos}).pre);
        Data.dprime.easy.(SET.positionExp {indPos}).P_Bias_pre = p;
        [p,h,stat] = signrank(Data.dprime.easy.(SET.positionExp {indPos}).post);
        Data.dprime.easy.(SET.positionExp {indPos}).P_Bias_post = p;
        [p,h,stat] = signrank(Data.dprime.difficult.(SET.positionExp {indPos}).pre);
        Data.dprime.difficult.(SET.positionExp {indPos}).P_Bias_pre = p;
        [p,h,stat] = signrank( Data.dprime.difficult.(SET.positionExp {indPos}).post);
        Data.dprime.difficult.(SET.positionExp {indPos}).P_Bias_post = p;
    end
end
%% Determine the response Bias by testing the criterion again the H0?
for indPos = 1:length(SET.positionExp ) %contra vs ipsi selection
    for indTime = 1:length(SET.TimeLine) %pre vs post
        [p,h,stat] = signrank(Data.criterion.easy.(SET.positionExp {indPos}).pre);
        Data.criterion.easy.(SET.positionExp {indPos}).P_Bias_pre = p;
        [p,h,stat] = signrank(Data.criterion.easy.(SET.positionExp {indPos}).post);
        Data.criterion.easy.(SET.positionExp {indPos}).P_Bias_post = p;
        [p,h,stat] = signrank(Data.criterion.difficult.(SET.positionExp {indPos}).pre);
        Data.criterion.difficult.(SET.positionExp {indPos}).P_Bias_pre = p;
        [p,h,stat] = signrank(Data.criterion.difficult.(SET.positionExp {indPos}).post);
        Data.criterion.difficult.(SET.positionExp {indPos}).P_Bias_post = p;
    end
end

%% CRITERION - WHICH RESPONSE BIAS?

for indPos = 1:length(SET.positionExp ) %contra vs ipsi selection
    if ResponseBias_DisplayColor == 1
        Color_MoreContra    = [0 1 1]; %cyan
        Color_LessContra    = [1 0 1]; %magenta
        Color_Neutral       = [0.5 0.5 0.5];
    elseif ResponseBias_DisplayColor == 0
        Color_MoreContra    = [0 0 0];
        Color_LessContra    = [0 0 0];
        Color_Neutral       = [0 0 0];
    end
    for indDiff = 1: length(SET.Difficulty)
        if  Data.criterion.(SET.Difficulty{indDiff}).(SET.positionExp {indPos}).P_Bias_pre < 0.05
            if strcmp(SET.positionExp {indPos}, 'contra')
                if nanmean(Data.criterion.(SET.Difficulty{indDiff}).(SET.positionExp {indPos}).pre) < 0 %negative Criterion
                    Data.criterion.(SET.Difficulty{indDiff}).(SET.positionExp {indPos}).BiasDirection_pre = 'moreContra' ;
                    Data.criterion.(SET.Difficulty{indDiff}).(SET.positionExp {indPos}).Bias_pre = 'Go' ;
                    Data.criterion.(SET.Difficulty{indDiff}).(SET.positionExp {indPos}).Bias_pre_Color = Color_MoreContra ;
                    
                elseif nanmean(Data.criterion.(SET.Difficulty{indDiff}).(SET.positionExp {indPos}).pre) > 0
                    Data.criterion.(SET.Difficulty{indDiff}).(SET.positionExp {indPos}).BiasDirection_pre = 'lessContra' ;
                    Data.criterion.(SET.Difficulty{indDiff}).(SET.positionExp {indPos}).Bias_pre = 'NoGo' ;
                    Data.criterion.(SET.Difficulty{indDiff}).(SET.positionExp {indPos}).Bias_pre_Color = Color_LessContra ;
                    
                end
                
            elseif  strcmp(SET.positionExp {indPos}, 'ipsi')
                if nanmean(-Data.criterion.(SET.Difficulty{indDiff}).(SET.positionExp {indPos}).pre) < 0 %negative Criterion
                    Data.criterion.(SET.Difficulty{indDiff}).(SET.positionExp {indPos}).BiasDirection_pre = 'moreContra' ;
                    Data.criterion.(SET.Difficulty{indDiff}).(SET.positionExp {indPos}).Bias_pre = 'NoGo' ;
                    Data.criterion.(SET.Difficulty{indDiff}).(SET.positionExp {indPos}).Bias_pre_Color = Color_MoreContra ;
                    
                elseif nanmean(-Data.criterion.(SET.Difficulty{indDiff}).(SET.positionExp {indPos}).pre) > 0
                    Data.criterion.(SET.Difficulty{indDiff}).(SET.positionExp {indPos}).BiasDirection_pre = 'lessContra' ;
                    Data.criterion.(SET.Difficulty{indDiff}).(SET.positionExp {indPos}).Bias_pre = 'Go' ;
                    Data.criterion.(SET.Difficulty{indDiff}).(SET.positionExp {indPos}).Bias_pre_Color = Color_LessContra ;
                    
                end
                
            end
        else
            Data.criterion.(SET.Difficulty{indDiff}).(SET.positionExp {indPos}).BiasDirection_pre = 'neutral' ;
            Data.criterion.(SET.Difficulty{indDiff}).(SET.positionExp {indPos}).Bias_pre = 'neutral' ;
            Data.criterion.(SET.Difficulty{indDiff}).(SET.positionExp {indPos}).Bias_pre_Color = Color_Neutral ;
            
        end
    end
end
%%% POST
for indPos = 1:length(SET.positionExp ) %contra vs ipsi selection
    if ResponseBias_DisplayColor == 1
        Color_MoreContra    = [0 1 1]; %cyan
        Color_LessContra    = [1 0 1]; %magenta
        Color_Neutral       = [0.5 0.5 0.5];
    elseif ResponseBias_DisplayColor == 0
        Color_MoreContra    = [0 0 1];
        Color_LessContra    = [0 0 1];
        Color_Neutral       = [0 0 1];
    end
    for indDiff = 1: length(SET.Difficulty)
        if  Data.criterion.(SET.Difficulty{indDiff}).(SET.positionExp {indPos}).P_Bias_post < 0.05
            if strcmp(SET.positionExp {indPos}, 'contra')
                if nanmean(Data.criterion.(SET.Difficulty{indDiff}).(SET.positionExp {indPos}).post) < 0 %negative Criterion
                    Data.criterion.(SET.Difficulty{indDiff}).(SET.positionExp {indPos}).BiasDirection_post = 'moreContra' ;
                    Data.criterion.(SET.Difficulty{indDiff}).(SET.positionExp {indPos}).Bias_post = 'Go' ;
                    Data.criterion.(SET.Difficulty{indDiff}).(SET.positionExp {indPos}).Bias_post_Color = Color_MoreContra ;
                    
                elseif nanmean(Data.criterion.(SET.Difficulty{indDiff}).(SET.positionExp {indPos}).post) > 0
                    Data.criterion.(SET.Difficulty{indDiff}).(SET.positionExp {indPos}).BiasDirection_post = 'lessContra' ;
                    Data.criterion.(SET.Difficulty{indDiff}).(SET.positionExp {indPos}).Bias_post = 'NoGo' ;
                    Data.criterion.(SET.Difficulty{indDiff}).(SET.positionExp {indPos}).Bias_post_Color = Color_LessContra ;
                    
                end
            elseif  strcmp(SET.positionExp {indPos}, 'ipsi')
                if nanmean(-Data.criterion.(SET.Difficulty{indDiff}).(SET.positionExp {indPos}).post) < 0 %negative Criterion
                    Data.criterion.(SET.Difficulty{indDiff}).(SET.positionExp {indPos}).BiasDirection_post = 'moreContra' ;
                    Data.criterion.(SET.Difficulty{indDiff}).(SET.positionExp {indPos}).Bias_post = 'NoGo' ;
                    Data.criterion.(SET.Difficulty{indDiff}).(SET.positionExp {indPos}).Bias_post_Color = Color_MoreContra ;
                    
                elseif nanmean(-Data.criterion.(SET.Difficulty{indDiff}).(SET.positionExp {indPos}).post) > 0
                    Data.criterion.(SET.Difficulty{indDiff}).(SET.positionExp {indPos}).BiasDirection_post = 'lessContra' ;
                    Data.criterion.(SET.Difficulty{indDiff}).(SET.positionExp {indPos}).Bias_post = 'Go' ;
                    Data.criterion.(SET.Difficulty{indDiff}).(SET.positionExp {indPos}).Bias_post_Color = Color_LessContra ;
                    
                end
                
            end
        else
            Data.criterion.(SET.Difficulty{indDiff}).(SET.positionExp {indPos}).BiasDirection_post = 'neutral' ;
            Data.criterion.(SET.Difficulty{indDiff}).(SET.positionExp {indPos}).Bias_post = 'neutral' ;
            Data.criterion.(SET.Difficulty{indDiff}).(SET.positionExp {indPos}).Bias_post_Color = Color_Neutral ;
        end
    end
end






%% GRAPH



%% Fig1: Criterion vs Dprime
figure('Position',[200 200 1200 900],'PaperPositionMode','auto'); % ,'PaperOrientation','landscape'
set(gcf,'Name','Dprime vs Criterion');
set(gcf,'Color',[1 1 1]);

MarkSize = 15;

ha(indPos) = subplot(1,2,1);
for i = 1: length(Data.dprime.easy.contra.pre)
    plot([Data.dprime.easy.contra.pre(i),Data.dprime.easy.contra.post(i)],[Data.criterion.easy.contra.pre(i),Data.criterion.easy.contra.post(i)], 'o','color',[Color.Contra, 0.4] , 'MarkerSize',MarkSize-8,'markerfacecolor',[1 1 1 ]); hold on;
    plot([Data.dprime.easy.ipsi.pre(i),Data.dprime.easy.ipsi.post(i)], [-Data.criterion.easy.ipsi.pre(i),-Data.criterion.easy.ipsi.post(i)], 'o','color',[Color.Ipsi, 0.4] , 'MarkerSize',MarkSize-8,'markerfacecolor',[1 1 1 ]); hold on;% reverse direction of criterion for ipsi
    plot(Data.dprime.easy.contra.post(i),Data.criterion.easy.contra.post(i), 'o','color',Color.Contra , 'MarkerSize',MarkSize-8,'markerfacecolor',Color.Contra); hold on;
    plot(Data.dprime.easy.ipsi.post(i),-Data.criterion.easy.ipsi.post(i), 'o','color',Color.Ipsi , 'MarkerSize',MarkSize-8,'markerfacecolor',Color.Ipsi); hold on;% reverse direction of criterion for ipsi
    
end


plot([nanmean(Data.dprime.easy.contra.pre),nanmean(Data.dprime.easy.contra.post)],[nanmean(Data.criterion.easy.contra.pre),nanmean(Data.criterion.easy.contra.post)], 'o-','color',Color.Contra , 'MarkerSize',MarkSize,'markerfacecolor',[1 1 1 ]); hold on;
plot([nanmean(Data.dprime.easy.ipsi.pre),nanmean(Data.dprime.easy.ipsi.post)], [-nanmean(Data.criterion.easy.ipsi.pre),-nanmean(Data.criterion.easy.ipsi.post)], 'o-','color',Color.Ipsi , 'MarkerSize',MarkSize,'markerfacecolor',[1 1 1 ]); hold on;% reverse direction of criterion for ipsi

plot(nanmean(Data.dprime.easy.contra.post),nanmean(Data.criterion.easy.contra.post), 'o','color',Color.Contra , 'MarkerSize',MarkSize,'markerfacecolor',Color.Contra ); hold on;
plot(nanmean(Data.dprime.easy.ipsi.post),-nanmean(Data.criterion.easy.ipsi.post), 'o','color',Color.Ipsi , 'MarkerSize',MarkSize,'markerfacecolor',Color.Ipsi); hold on;% reverse direction of criterion for ipsi
if Data.dprime.easy.ipsi.pvalue < 0.05
    ymax = max(nanmean(Data.criterion.easy.ipsi.pre) ,nanmean([Data.criterion.easy.ipsi.post]));
    ext_sigline([nanmean([Data.dprime.easy.ipsi.pre]),nanmean([Data.dprime.easy.ipsi.post])],Data.dprime.easy.ipsi.pvalue,[ ],(ymax +0.5 ),'x', Color.Ipsi); hold on;
end
if Data.criterion.easy.ipsi.pvalue < 0.05
    y1 = -1*nanmean(Data.criterion.easy.ipsi.pre) ;
    y2 = -1*nanmean([Data.criterion.easy.ipsi.post]);
    ymax = max(nanmean(Data.dprime.easy.ipsi.pre) ,nanmean([Data.dprime.easy.ipsi.post])) +0.5 ;
    ext_sigline([y1,y2],Data.criterion.easy.ipsi.pvalue,[],ymax,'y', Color.Ipsi); hold on;
end


if Data.dprime.easy.contra.pvalue < 0.05
    ymax = max(nanmean(Data.criterion.easy.contra.pre) ,nanmean([Data.criterion.easy.contra.post]));
    ext_sigline([nanmean([Data.dprime.easy.contra.pre]),nanmean([Data.dprime.easy.contra.post])],Data.dprime.easy.contra.pvalue,[ ],(ymax +0.5),'x', Color.Contra); hold on;
end
if Data.criterion.easy.contra.pvalue < 0.05
    y1 = nanmean(Data.criterion.easy.contra.pre) ;
    y2 = nanmean([Data.criterion.easy.contra.post]);
    ymax = max(nanmean(Data.dprime.easy.contra.pre) ,nanmean([Data.dprime.easy.contra.post])) +0.5;
    ext_sigline([y1,y2],Data.criterion.easy.contra.pvalue,[],ymax,'y', Color.Contra); hold on;
end
axis square
xlabel('sensitivity','fontsize',fs,'fontweight','b', 'Interpreter', 'none')
ylabel('criterion','fontsize',fs,'fontweight','b', 'Interpreter', 'none')
set(gca,'ylim',[-3 3],'xlim',[0 5],'fontsize',fs)
title(' easy distractor')
text(0.5,-2.8, 'more Contra (ipsi:NoGo, contra:Go)', 'Color', 'k')
text(0.5,2.8, 'less Contra (ipsi:Go, contra:NoGo)', 'Color', 'k')

ha(indPos) = subplot(1,2,2);
for i = 1: length(Data.dprime.difficult.contra.pre)
    plot([Data.dprime.difficult.contra.pre(i),Data.dprime.difficult.contra.post(i)],[Data.criterion.difficult.contra.pre(i),Data.criterion.difficult.contra.post(i)], 'o','color',[Color.Contra, 0.4] , 'MarkerSize',MarkSize-8,'markerfacecolor',[1 1 1 ]); hold on;
    plot([Data.dprime.difficult.ipsi.pre(i),Data.dprime.difficult.ipsi.post(i)], [-Data.criterion.difficult.ipsi.pre(i),-Data.criterion.difficult.ipsi.post(i)], 'o','color',[Color.Ipsi, 0.4] , 'MarkerSize',MarkSize-8,'markerfacecolor',[1 1 1 ]); hold on;% reverse direction of criterion for ipsi
    plot(Data.dprime.difficult.contra.post(i),Data.criterion.difficult.contra.post(i), 'o','color',Color.Contra , 'MarkerSize',MarkSize-8,'markerfacecolor',Color.Contra); hold on;
    plot(Data.dprime.difficult.ipsi.post(i),-Data.criterion.difficult.ipsi.post(i), 'o','color',Color.Ipsi , 'MarkerSize',MarkSize-8,'markerfacecolor',Color.Ipsi); hold on;% reverse direction of criterion for ipsi
    
end
plot([nanmean(Data.dprime.difficult.contra.pre),nanmean(Data.dprime.difficult.contra.post)],[nanmean(Data.criterion.difficult.contra.pre),nanmean(Data.criterion.difficult.contra.post)], 'o-','color',Color.Contra , 'MarkerSize',MarkSize,'markerfacecolor',[1 1 1 ]); hold on;
plot([nanmean(Data.dprime.difficult.ipsi.pre),nanmean(Data.dprime.difficult.ipsi.post)],[-nanmean(Data.criterion.difficult.ipsi.pre),-nanmean(Data.criterion.difficult.ipsi.post)], 'o-','color',Color.Ipsi , 'MarkerSize',MarkSize,'markerfacecolor',[1 1 1 ]); hold on;% reverse direction of criterion for ipsi

plot(nanmean(Data.dprime.difficult.contra.post),nanmean(Data.criterion.difficult.contra.post), 'o','color',Color.Contra , 'MarkerSize',MarkSize,'markerfacecolor',Color.Contra); hold on;
plot(nanmean(Data.dprime.difficult.ipsi.post),-nanmean(Data.criterion.difficult.ipsi.post), 'o','color',Color.Ipsi , 'MarkerSize',MarkSize,'markerfacecolor',Color.Ipsi); hold on;% reverse direction of criterion for ipsi
if Data.dprime.difficult.ipsi.pvalue < 0.05
    ymax = max(nanmean(Data.criterion.difficult.ipsi.pre) ,nanmean([Data.criterion.difficult.ipsi.post]));
    ext_sigline([nanmean([Data.dprime.difficult.ipsi.pre]),nanmean([Data.dprime.difficult.ipsi.post])],Data.dprime.difficult.ipsi.pvalue,[ ],(ymax ),'x', Color.Ipsi); hold on;
end
if Data.criterion.difficult.ipsi.pvalue < 0.05
    y1 = -1*nanmean(Data.criterion.difficult.ipsi.pre) ;
    y2 = -1*nanmean([Data.criterion.difficult.ipsi.post]);
    ymax = max(nanmean(Data.dprime.difficult.ipsi.pre) ,nanmean([Data.dprime.difficult.ipsi.post])) ;
    ext_sigline([y1,y2],Data.criterion.difficult.ipsi.pvalue,[],ymax,'y', Color.Ipsi); hold on;
end


if Data.dprime.difficult.contra.pvalue < 0.05
    ymax = max(nanmean(Data.criterion.difficult.contra.pre) ,nanmean([Data.criterion.difficult.contra.post]));
    ext_sigline([nanmean([Data.dprime.difficult.contra.pre]),nanmean([Data.dprime.difficult.contra.post])],Data.dprime.difficult.contra.pvalue,[ ],(ymax),'x', Color.Contra); hold on;
end
if Data.criterion.difficult.contra.pvalue < 0.05
    y1 = nanmean(Data.criterion.difficult.contra.pre) ;
    y2 = nanmean([Data.criterion.difficult.contra.post]);
    ymax = max(nanmean(Data.dprime.difficult.contra.pre) ,nanmean([Data.dprime.difficult.contra.post])) ;
    ext_sigline([y1,y2],Data.criterion.difficult.contra.pvalue,[],ymax +0.6,'y', Color.Contra); hold on;
end
axis square
xlabel('sensitivity','fontsize',fs,'fontweight','b', 'Interpreter', 'none')
ylabel('criterion','fontsize',fs,'fontweight','b', 'Interpreter', 'none')
set(gca,'ylim',[-3 3],'xlim',[-1 4],'fontsize',fs)
text(-0.5,-2.8, 'more Contra (ipsi:NoGo, contra:Go)', 'Color', 'k')
text(-0.5,2.8, 'less Contra (ipsi:Go, contra:NoGo)', 'Color', 'k')


%% save the graphs
if SaveGraph
    h = figure(1);
    %savefig(h, [path_SaveFig , 'fig',filesep, 'SingleStimuli_Crit_Dprime_' monkey,'_' experiment ,'_' folder '.fig'])
    print(h,[path_SaveFig , ['png', ControlFolder] ,filesep,  'SingleStimuli_Crit_Dprime' monkey,'_' , '.png'], '-dpng')%experiment,'_' folder
    set(h,'Renderer','Painters');
    set(h,'PaperPositionMode','auto')
    compl_filename =  [path_SaveFig , ['ai', ControlFolder] ,filesep,  'SingleStimuli_Crit_Dprime' monkey,'_' experiment ,'_', folder ,'.ai'] ;
    print(h,'-depsc',compl_filename);
    close all;
end


%% graph - From the paper Luo & Maunsell (2018)

figure('Position',[200 200 1200 900],'PaperPositionMode','auto'); % ,'PaperOrientation','landscape'
set(gcf,'Name',['Difficult_Graph1HR_FAR_Graph2Sensitivity_Criterion',ControlFolder ]);
set(gcf,'Color',[1 1 1]);


ha(2) = subplot(1,2,1);
%ipsi Post

for i = 1: length(Data.pFA.difficult.ipsi.pre)
    plot(Data.pFA.difficult.ipsi.post(i), Data.pHit.ipsi.post(i), 'o' ,'color',Color.Ipsi, 'MarkerSize',MarkSize_GraphFAR_HR_PerSession, 'MarkerFaceColor',Color.Ipsi,'LineWidth', 2); hold on;
    plot(Data.pFA.difficult.contra.post(i), Data.pHit.contra.post(i), 'o','color',Color.Contra , 'MarkerSize',MarkSize_GraphFAR_HR_PerSession,'markerfacecolor',Color.Contra ,'LineWidth', 2); hold on;
    
    line([Data.pFA.difficult.ipsi.pre(i),Data.pFA.difficult.ipsi.post(i)], [Data.pHit.ipsi.pre(i),Data.pHit.ipsi.post(i)],'Color',[Color.Ipsi 0.3] , 'MarkerSize',MarkSize_GraphFAR_HR_PerSession,'markerfacecolor',[1 1 1],'LineWidth', 2); hold on;
    line([Data.pFA.difficult.contra.pre(i),Data.pFA.difficult.contra.post(i)], [Data.pHit.contra.pre(i),Data.pHit.contra.post(i)],'Color',[Color.Contra 0.3] , 'MarkerSize',MarkSize_GraphFAR_HR_PerSession,'markerfacecolor',[1 1 1],'LineWidth', 2); hold on;
    
    %Fill the circle which is post
end
plot([nanmean(Data.pFA.difficult.ipsi.pre),nanmean(Data.pFA.difficult.ipsi.post)], [nanmean(Data.pHit.ipsi.pre),nanmean(Data.pHit.ipsi.post)], 'o-','color',Color.Ipsi , 'MarkerSize',MarkSize_GraphFAR_HR,'markerfacecolor',[1 1 1],'LineWidth', LineWidthSize); hold on;
plot([nanmean(Data.pFA.difficult.contra.pre),nanmean(Data.pFA.difficult.contra.post)], [nanmean(Data.pHit.contra.pre),nanmean(Data.pHit.contra.post)], 'o-','color',Color.Contra , 'MarkerSize',MarkSize_GraphFAR_HR,'markerfacecolor',[1 1 1],'LineWidth', LineWidthSize); hold on;
plot(nanmean(Data.pFA.difficult.ipsi.post), nanmean(Data.pHit.ipsi.post), 'o','color',Color.Ipsi , 'MarkerSize',MarkSize_GraphFAR_HR-1,'markerfacecolor',Color.Ipsi,'LineWidth', 2); hold on;
plot(nanmean(Data.pFA.difficult.contra.post), nanmean(Data.pHit.contra.post), 'o','color',Color.Contra , 'MarkerSize',MarkSize_GraphFAR_HR-1,'markerfacecolor',Color.Contra,'LineWidth', 2); hold on;
%legend('ipsi', 'contra', 'Location', 'South')


line([0 1],[1 0],'Color',[0 0 0],'LineStyle',':');
set(gca,'ylim',[0 1],'xlim',[0 1],'fontsize',fs)
xlabel( 'Data.FA rate','fontsize',fs,'fontweight','b', 'Interpreter', 'none' );
ylabel( 'Hitrate','fontsize',fs,'fontweight','b', 'Interpreter', 'none' );
axis square


if Data.pFA.difficult.ipsi.pvalue < 0.05
    y1 = nanmean(Data.pFA.difficult.ipsi.pre) ;
    y2 = nanmean(Data.pFA.difficult.ipsi.post);
    ymax = min(nanmean(Data.pHit.ipsi.pre) ,nanmean(Data.pHit.ipsi.post));
    ext_sigline([y1,y2],Data.pFA.difficult.ipsi.pvalue,[],ymax -0.3,'x', Color.Ipsi); hold on;
end
if Data.pHit.ipsi.pvalue < 0.05
    y1 = nanmean(Data.pHit.ipsi.pre) ;
    y2 = nanmean(Data.pHit.ipsi.post);
    ymax = max(nanmean(Data.pFA.difficult.ipsi.pre) ,nanmean(Data.pFA.difficult.ipsi.post));
    ext_sigline([y1,y2],Data.pHit.ipsi.pvalue,[],ymax +0.2,'y',Color.Ipsi); hold on;
end
if Data.pFA.difficult.contra.pvalue < 0.05
    y1 = nanmean(Data.pFA.difficult.contra.pre) ;
    y2 = nanmean(Data.pFA.difficult.contra.post);
    ymax = min(nanmean(Data.pHit.contra.pre) ,nanmean(Data.pHit.contra.post));
    ext_sigline([y1,y2],Data.pFA.difficult.contra.pvalue,[],ymax -0.2,'x', Color.Contra); hold on;
end
if Data.pHit.contra.pvalue < 0.05
    y1 = nanmean(Data.pHit.contra.pre) ;
    y2 = nanmean(Data.pHit.contra.post);
    ymax = max(nanmean(Data.pFA.difficult.contra.pre) ,nanmean(Data.pFA.difficult.contra.post));
    ext_sigline([y1,y2],Data.pHit.contra.pvalue,[],ymax +0.2,'y',Color.Contra); hold on;
end


ha(indPos) = subplot(1,2,2);
MarkSize = 15;
for i = 1: length(Data.dprime.difficult.contra.pre)
    plot([Data.dprime.difficult.contra.pre(i),Data.dprime.difficult.contra.post(i)],[Data.criterion.difficult.contra.pre(i),Data.criterion.difficult.contra.post(i)], 'o','color',[Color.Contra 0.4] , 'MarkerSize',MarkSize_CritDpr_small,'markerfacecolor',[1 1 1 ]); hold on;
    plot([Data.dprime.difficult.ipsi.pre(i),Data.dprime.difficult.ipsi.post(i)], [-Data.criterion.difficult.ipsi.pre(i),-Data.criterion.difficult.ipsi.post(i)], 'o','color',[Color.Ipsi 0.4] , 'MarkerSize',MarkSize_CritDpr_small,'markerfacecolor',[1 1 1 ]); hold on;% reverse direction of criterion for ipsi
    plot(Data.dprime.difficult.contra.post(i),Data.criterion.difficult.contra.post(i), 'o','color',Color.Contra , 'MarkerSize',MarkSize_CritDpr_small,'markerfacecolor',Color.Contra); hold on;
    plot(Data.dprime.difficult.ipsi.post(i),-Data.criterion.difficult.ipsi.post(i), 'o','color',Color.Ipsi , 'MarkerSize',MarkSize_CritDpr_small,'markerfacecolor',Color.Ipsi); hold on;% reverse direction of criterion for ipsi
end
plot([nanmean(Data.dprime.difficult.contra.pre),nanmean(Data.dprime.difficult.contra.post)],[nanmean(Data.criterion.difficult.contra.pre) ,nanmean(Data.criterion.difficult.contra.post)], 'o-','color',Color.Contra ,'LineWidth',LineWidthSize , 'MarkerSize',MarkSize_CritDpr,'markerfacecolor',[1 1 1 ]); hold on;
plot([nanmean(Data.dprime.difficult.ipsi.pre),nanmean(Data.dprime.difficult.ipsi.post)],[ -nanmean(Data.criterion.difficult.ipsi.pre) ,-nanmean(Data.criterion.difficult.ipsi.post)], 'o-','color',Color.Ipsi ,'LineWidth',LineWidthSize , 'MarkerSize',MarkSize_CritDpr,'markerfacecolor',[1 1 1 ]); hold on;% reverse direction of criterion for ipsi

plot(nanmean(Data.dprime.difficult.contra.post),nanmean(Data.criterion.difficult.contra.post), 'o','color',Color.Contra , 'MarkerSize',MarkSize_CritDpr,'markerfacecolor',Color.Contra); hold on;
plot(nanmean(Data.dprime.difficult.ipsi.post),-nanmean(Data.criterion.difficult.ipsi.post), 'o','color',Color.Ipsi , 'MarkerSize',MarkSize_CritDpr,'markerfacecolor',Color.Ipsi); hold on;% reverse direction of criterion for ipsi

if Data.dprime.difficult.ipsi.pvalue < 0.05
    ymax = max(-1*nanmean(Data.criterion.difficult.ipsi.pre) ,-1*nanmean([Data.criterion.difficult.ipsi.post]));
    ext_sigline([nanmean([Data.dprime.difficult.ipsi.pre]),nanmean([Data.dprime.difficult.ipsi.post])],Data.dprime.difficult.ipsi.pvalue,[ ],(ymax+ 0.5 ),'x', Color.Ipsi); hold on;
end
if Data.criterion.difficult.ipsi.pvalue < 0.05
    y1 = -1*nanmean(Data.criterion.difficult.ipsi.pre) ;
    y2 = -1*nanmean([Data.criterion.difficult.ipsi.post]);
    ymax = max(nanmean(Data.dprime.difficult.ipsi.pre) ,nanmean([Data.dprime.difficult.ipsi.post])) ;
    ext_sigline([y1,y2],Data.criterion.difficult.ipsi.pvalue,[],ymax+ 0.6,'y', Color.Ipsi); hold on;
end


if Data.dprime.difficult.contra.pvalue < 0.05
    ymax = max(nanmean(Data.criterion.difficult.contra.pre) ,nanmean([Data.criterion.difficult.contra.post]));
    ext_sigline([nanmean([Data.dprime.difficult.contra.pre]),nanmean([Data.dprime.difficult.contra.post])],Data.dprime.difficult.contra.pvalue,[ ],(ymax+ 0.8 ),'x', Color.Contra); hold on;
end
if Data.criterion.difficult.contra.pvalue < 0.05
    y1 = nanmean(Data.criterion.difficult.contra.pre) ;
    y2 = nanmean([Data.criterion.difficult.contra.post]);
    ymax = max(nanmean(Data.dprime.difficult.contra.pre) ,nanmean([Data.dprime.difficult.contra.post])) ;
    ext_sigline([y1,y2],Data.criterion.difficult.contra.pvalue,[],ymax+ 0.6,'y', Color.Contra); hold on;
end
axis square
xlabel('sensitivity','fontsize',fs,'fontweight','b', 'Interpreter', 'none')
ylabel('criterion','fontsize',fs,'fontweight','b', 'Interpreter', 'none')
set(gca,'ylim',[-3 3],'xlim',[xlim_SDT_diff],'fontsize',fs)

text(xlim_SDT_diff(1) +0.1,-2.8, 'more Contra (ipsi:NoGo, contra:Go)', 'Color', 'k','fontsize',18)
text(xlim_SDT_diff(1) +0.1,2.8, 'less Contra (ipsi:Go, contra:NoGo)', 'Color', 'k','fontsize',18)

%% save the graphs
if SaveGraph
    h = figure(1);
    %savefig(h, [path_SaveFig ,filesep,'fig',filesep,  'DoubleStimuli_Crit_Dprime_' monkey,'_' experiment ,'_' folder '.fig'])
    print(h,[path_SaveFig ,filesep,  ['png', ControlFolder]  ,filesep,  'SingleStimuli_difficult_HR_FAR_Dprime_Criterion' monkey , '.png'], '-dpng')% '_' experiment ,'_', folder
    set(h,'Renderer','Painters');
    set(h,'PaperPositionMode','auto')
    compl_filename =  [path_SaveFig , ['ai', ControlFolder]  ,filesep,  'SingleStimuli_difficult_HR_FAR_Dprime_Criterion' monkey,'.ai'] ;
    print(h,'-depsc',compl_filename);
    close all;
end


%% EASY - TALK GRAPHS
figure('Position',[200 200 1200 900],'PaperPositionMode','auto'); % ,'PaperOrientation','landscape'
set(gcf,'Name',['Easy_Graph1HR_FAR_Graph2Sensitivity_Criterion',ControlFolder ]);
set(gcf,'Color',[1 1 1]);
 
ha(2) = subplot(1,2,1);
for i = 1: length(Data.pFA.easy.ipsi.pre)
    plot(Data.pFA.easy.ipsi.post(i), Data.pHit.ipsi.post(i), 'o' ,'color',Color.Ipsi, 'MarkerSize',MarkSize_GraphFAR_HR_PerSession, 'MarkerFaceColor',Color.Ipsi ,'LineWidth', 2); hold on;
    plot(Data.pFA.easy.contra.post(i), Data.pHit.contra.post(i), 'o','color',Color.Contra , 'MarkerSize',MarkSize_GraphFAR_HR_PerSession,'markerfacecolor',Color.Contra,'LineWidth', 2); hold on;
    
    line([Data.pFA.easy.ipsi.pre(i),Data.pFA.easy.ipsi.post(i)], [Data.pHit.ipsi.pre(i),Data.pHit.ipsi.post(i)],'Color',[Color.Ipsi 0.3] , 'MarkerSize',MarkSize_GraphFAR_HR_PerSession,'markerfacecolor',[1 1 1],'LineWidth', 1); hold on;
    line([Data.pFA.easy.contra.pre(i),Data.pFA.easy.contra.post(i)], [Data.pHit.contra.pre(i),Data.pHit.contra.post(i)],'Color',[Color.Contra 0.3] , 'MarkerSize',MarkSize_GraphFAR_HR_PerSession,'markerfacecolor',[1 1 1],'LineWidth', 1); hold on;
    
    %Fill the circle which is post
end
plot([nanmean(Data.pFA.easy.ipsi.pre),nanmean(Data.pFA.easy.ipsi.post)], [nanmean(Data.pHit.ipsi.pre),nanmean(Data.pHit.ipsi.post)], 'o-','color',Color.Ipsi , 'MarkerSize',MarkSize_GraphFAR_HR,'markerfacecolor',[1 1 1],'LineWidth',LineWidthSize); hold on;
plot(nanmean(Data.pFA.easy.ipsi.post), nanmean(Data.pHit.ipsi.post), 'o','color',Color.Ipsi , 'MarkerSize',MarkSize_GraphFAR_HR-1,'markerfacecolor',Color.Ipsi,'LineWidth', 2); hold on;

plot([nanmean(Data.pFA.easy.contra.pre),nanmean(Data.pFA.easy.contra.post)], [nanmean(Data.pHit.contra.pre),nanmean(Data.pHit.contra.post)], 'o-','color',Color.Contra, 'MarkerSize',MarkSize_GraphFAR_HR,'markerfacecolor',[1 1 1],'LineWidth', LineWidthSize); hold on;
plot(nanmean(Data.pFA.easy.contra.post), nanmean(Data.pHit.contra.post), 'o','color',Color.Contra , 'MarkerSize',MarkSize_GraphFAR_HR-1,'markerfacecolor',Color.Contra,'LineWidth', 2); hold on;

%legend('ipsi', 'contra', 'Location', 'South')


line([0 1],[1 0],'Color',[0 0 0],'LineStyle',':');
set(gca,'ylim',[0 1],'xlim',[0 1],'fontsize',fs)
xlabel( 'FA rate','fontsize',fs,'fontweight','b', 'Interpreter', 'none' );
ylabel( 'Hitrate','fontsize',fs,'fontweight','b', 'Interpreter', 'none' );
axis square


if Data.pFA.easy.ipsi.pvalue < 0.05
    y1 = nanmean(Data.pFA.easy.ipsi.pre) ;
    y2 = nanmean(Data.pFA.easy.ipsi.post);
    ymax = min(nanmean(Data.pHit.ipsi.pre) ,nanmean(Data.pHit.ipsi.post));
    ext_sigline([y1,y2],Data.pFA.easy.ipsi.pvalue,[],ymax -0.3,'x', Color.Ipsi); hold on;
end
if Data.pHit.ipsi.pvalue < 0.05
    y1 = nanmean(Data.pHit.ipsi.pre) ;
    y2 = nanmean(Data.pHit.ipsi.post);
    ymax = max(nanmean(Data.pFA.easy.ipsi.pre) ,nanmean(Data.pFA.easy.ipsi.post));
    ext_sigline([y1,y2],Data.pHit.ipsi.pvalue,[],ymax +0.2,'y',Color.Ipsi); hold on;
end
if Data.pFA.easy.contra.pvalue < 0.05
    y1 = nanmean(Data.pFA.easy.contra.pre) ;
    y2 = nanmean(Data.pFA.easy.contra.post);
    ymax = min(nanmean(Data.pHit.contra.pre) ,nanmean(Data.pHit.contra.post));
    ext_sigline([y1,y2],Data.pFA.easy.contra.pvalue,[],ymax -0.2,'x', Color.Contra); hold on;
end
if Data.pHit.contra.pvalue < 0.05
    y1 = nanmean(Data.pHit.contra.pre) ;
    y2 = nanmean(Data.pHit.contra.post);
    ymax = max(nanmean(Data.pFA.easy.contra.pre) ,nanmean(Data.pFA.easy.contra.post));
    ext_sigline([y1,y2],Data.pHit.contra.pvalue,[],ymax +0.2,'y',Color.Contra); hold on;
end


ha(indPos) = subplot(1,2,2);

for i = 1: length(Data.dprime.easy.contra.pre)
    plot([Data.dprime.easy.contra.pre(i),Data.dprime.easy.contra.post(i)],[Data.criterion.easy.contra.pre(i),Data.criterion.easy.contra.post(i)], 'o','color',[Color.Contra, 0.4] , 'MarkerSize',MarkSize_CritDpr_small,'markerfacecolor',[1 1 1 ]); hold on;
    plot([Data.dprime.easy.ipsi.pre(i),Data.dprime.easy.ipsi.post(i)], [-Data.criterion.easy.ipsi.pre(i),-Data.criterion.easy.ipsi.post(i)], 'o','color',[Color.Ipsi, 0.4] , 'MarkerSize',MarkSize_CritDpr_small,'markerfacecolor',[1 1 1 ]); hold on;% reverse direction of criterion for ipsi
    plot(Data.dprime.easy.contra.post(i),Data.criterion.easy.contra.post(i), 'o','color',Color.Contra , 'MarkerSize',MarkSize_CritDpr_small,'markerfacecolor',Color.Contra); hold on;
    plot(Data.dprime.easy.ipsi.post(i),-Data.criterion.easy.ipsi.post(i), 'o','color',Color.Ipsi , 'MarkerSize',MarkSize_CritDpr_small,'markerfacecolor',Color.Ipsi); hold on;% reverse direction of criterion for ipsi
end
plot([nanmean(Data.dprime.easy.ipsi.pre),nanmean(Data.dprime.easy.ipsi.post)],[ -nanmean(Data.criterion.easy.ipsi.pre) ,-nanmean(Data.criterion.easy.ipsi.post)], 'o-','color',Color.Ipsi , 'MarkerSize',MarkSize_CritDpr,'LineWidth', LineWidthSize ,'markerfacecolor',[1 1 1 ]); hold on;% reverse direction of criterion for ipsi
plot(nanmean(Data.dprime.easy.ipsi.post),-nanmean(Data.criterion.easy.ipsi.post), 'o','color',Color.Ipsi , 'MarkerSize',MarkSize_CritDpr,'markerfacecolor',Color.Ipsi); hold on;% reverse direction of criterion for ipsi

plot([nanmean(Data.dprime.easy.contra.pre),nanmean(Data.dprime.easy.contra.post)],[nanmean(Data.criterion.easy.contra.pre) ,nanmean(Data.criterion.easy.contra.post)], 'o-','color',Color.Contra , 'MarkerSize',MarkSize_CritDpr,'LineWidth',LineWidthSize , 'markerfacecolor',[1 1 1 ]); hold on;
plot(nanmean(Data.dprime.easy.contra.post),nanmean(Data.criterion.easy.contra.post), 'o','color',Color.Contra , 'MarkerSize',MarkSize_CritDpr,'markerfacecolor',Color.Contra); hold on;


ymax_Cri_ipsi = max(-1*nanmean(Data.criterion.easy.ipsi.pre) ,-1*nanmean([Data.criterion.easy.ipsi.post]));
ymax_Cri_contra = max(nanmean(Data.criterion.easy.contra.pre) ,nanmean([Data.criterion.easy.contra.post]));

if Data.dprime.easy.ipsi.pvalue < 0.05
    if ymax_Cri_ipsi > ymax_Cri_contra
        ymax_Cri_ipsi = ymax_Cri_ipsi + 0.5 ;
    else
        ymax_Cri_ipsi = ymax_Cri_ipsi - 0.5 ;
        
    end
    ext_sigline([nanmean([Data.dprime.easy.ipsi.pre]),nanmean([Data.dprime.easy.ipsi.post])],Data.dprime.easy.ipsi.pvalue,[ ],(ymax_Cri_ipsi +0.5 ),'x', Color.Ipsi); hold on;
end
if Data.criterion.easy.ipsi.pvalue < 0.05
    y1 = -1*nanmean(Data.criterion.easy.ipsi.pre) ;
    y2 = -1*nanmean([Data.criterion.easy.ipsi.post]);
    ymax = max(nanmean(Data.dprime.easy.ipsi.pre) ,nanmean([Data.dprime.easy.ipsi.post])) ;
    ext_sigline([y1,y2],Data.criterion.easy.ipsi.pvalue,[],ymax+ 0.5,'y',Color.Ipsi); hold on;
end


if Data.dprime.easy.contra.pvalue < 0.05
    if ymax_Cri_ipsi > ymax_Cri_contra
        ymax_Cri_contra = ymax_Cri_contra - 0.8 ;
    else
        ymax_Cri_contra = ymax_Cri_contra + 0.8 ;
        
    end
    ext_sigline([nanmean([Data.dprime.easy.contra.pre]),nanmean([Data.dprime.easy.contra.post])],Data.dprime.easy.contra.pvalue,[ ],ymax_Cri_contra,'x', Color.Contra); hold on;
end
if Data.criterion.easy.contra.pvalue < 0.05
    y1 = nanmean(Data.criterion.easy.contra.pre) ;
    y2 = nanmean([Data.criterion.easy.contra.post]);
    ymax = max(nanmean(Data.dprime.easy.contra.pre) ,nanmean([Data.dprime.easy.contra.post])) ;
    ext_sigline([y1,y2],Data.criterion.easy.contra.pvalue,[],ymax+ 0.3,'y', Color.Contra); hold on;
end
axis square
xlabel('Sensitivity','fontsize',fs,'fontweight','b', 'Interpreter', 'none')
ylabel('Criterion','fontsize',fs,'fontweight','b', 'Interpreter', 'none')
set(gca,'ylim',[-3 3],'xlim',xlim_SDT_easy,'fontsize',fs)


text(xlim_SDT_easy(1) +0.1,-2.8, 'more Contra (ipsi:NoGo, contra:Go)', 'Color', 'k','fontsize',18)
text(xlim_SDT_easy(1) + 0.1,2.8, 'less Contra (ipsi:Go, contra:NoGo)', 'Color', 'k','fontsize',18)

%% save the graphs
if SaveGraph
    h= figure(1);
    %savefig(h, [path_SaveFig ,filesep,'fig',filesep,  'DoubleStimuli_Crit_Dprime_' monkey,'_' experiment ,'_' folder '.fig'])
    print(h,[path_SaveFig ,filesep,  ['png', ControlFolder]  ,filesep,  'SingleStimuli_easy_HR_FAR_Dprime_Criterion' monkey, '.png'], '-dpng') %'_' experiment,'_' folder
    set(h,'Renderer','Painters');
    set(h,'PaperPositionMode','auto')
    compl_filename =  [path_SaveFig , ['ai', ControlFolder]  ,filesep,  'SingleStimuli_easy_HR_FAR_Dprime_Criterion' monkey,'.ai'] ; %'_' experiment ,'_', folder ,
    print(h,'-depsc',compl_filename);
    close all;
end


%% FIXATION
figure('Position',[200 200 1200 900],'PaperPositionMode','auto'); % ,'PaperOrientation','landscape'
set(gcf,'Name',['Fixation - Miss',ControlFolder ]);
set(gcf,'Color',[1 1 1]);

ha(1) = subplot(2,2,1); %1 % 5

plot([1;2],[nanmean(Data.pFix.easy.ipsi.pre); nanmean(Data.pFix.easy.ipsi.post)], 'x-','color',Color.Black ,'MarkerSize',10,'LineWidth',2,'markerfacecolor',[1 1 1 ]); hold on;
%plot([1;2],[nanmean(Data.pMiss.ipsi.pre); nanmean(Data.pMiss.ipsi.post)], 'x-','color',Color.Black ,'MarkerSize',10,'LineWidth',2,'markerfacecolor',[1 1 1 ]); hold on;
plot([1;2],[nanmean(Data.pHit.ipsi.pre); nanmean(Data.pHit.ipsi.post)], 'o-','color',Color.Ipsi ,'MarkerSize',10,'LineWidth',2,'markerfacecolor',[1 1 1 ]); hold on;
plot([1;2],[nanmean(Data.pFA.easy.ipsi.pre); nanmean(Data.pFA.easy.ipsi.post)], 's-','color',Color.Ipsi ,'MarkerSize',10,'LineWidth',2,'markerfacecolor',[1 1 1 ]); hold on;

plot(2, nanmean(Data.pHit.ipsi.post), 'o','color',Color.Ipsi ,'MarkerSize',10,'markerfacecolor',Color.Ipsi); hold on;
plot(2, nanmean(Data.pFix.easy.ipsi.post), 'x','color',Color.Black,'MarkerSize',10,'markerfacecolor',Color.Ipsi); hold on;
%plot(2, nanmean(Data.pMiss.ipsi.post), 'x','color',Color.Black,'MarkerSize',10,'markerfacecolor',Color.Ipsi); hold on;
plot(2, nanmean(Data.pFA.easy.ipsi.post), 's','color',Color.Ipsi ,'MarkerSize',10,'markerfacecolor',Color.Ipsi); hold on;

%legend('tar','fix', 'dis', 'Location', 'SouthEast')
for i = 1: length(Data.pFA.easy.ipsi.pre)
    plot([0.9;2.1],[Data.pHit.ipsi.pre(i); Data.pHit.ipsi.post(i)], 'o','color',Color.Ipsi ,'MarkerSize',8,'markerfacecolor',[1 1 1 ]); hold on;
    %plot([0.9;2.1],[Data.pMiss.ipsi.pre(i); Data.pMiss.ipsi.post(i)], 'X','color',Color.Black,'MarkerSize',8,'markerfacecolor',[1 1 1 ]); hold on;
    plot([0.9;2.1],[Data.pFix.easy.ipsi.pre(i); Data.pFix.easy.ipsi.post(i)], 'X','color',Color.Black,'MarkerSize',8,'markerfacecolor',[1 1 1 ]); hold on;
    plot([0.8;2.2],[Data.pFA.easy.ipsi.pre(i); Data.pFA.easy.ipsi.post(i)], 's','color',Color.Ipsi ,'MarkerSize',8,'markerfacecolor',[1 1 1 ]); hold on;
end

% if Data.pMiss.ipsi.pvalue < 0.05
%     ymax = nanmean([nanmean(Data.pMiss.ipsi.pre) ,nanmean(Data.pMiss.ipsi.post)]);
%     ext_sigline([1,2], Data.pMiss.ipsi.pvalue,[],ymax +0.1,'x',Color.Black )
% end

if Data.pFix.easy.ipsi.pvalue < 0.05
    ymax = nanmean([nanmean(Data.pFix.easy.ipsi.pre) ,nanmean(Data.pFix.easy.ipsi.post)]);
    ext_sigline([1,2], Data.pFix.easy.ipsi.pvalue,[],ymax +0.1,'x',Color.Black )
end
if Data.pHit.ipsi.pvalue < 0.05
    ymax = nanmean([nanmean(Data.pHit.ipsi.pre) ,nanmean(Data.pHit.ipsi.post)]);
    ext_sigline([1.4 1.6], Data.pHit.ipsi.pvalue,[],ymax ,'x', Color.Ipsi )
end

if Data.pFA.easy.ipsi.pvalue < 0.05
    ymax = nanmean([nanmean(Data.pFA.easy.ipsi.pre) ,nanmean(Data.pFA.easy.ipsi.post)]);
    ext_sigline([1.4 1.6], Data.pFA.easy.ipsi.pvalue,[],ymax,'x', Color.Ipsi )
end
ylabel( 'Selection (easy)','fontsize',20,'fontweight','b', 'Interpreter', 'none' );
set(gca,'xlim',[0 3],'ylim',[0 1],'Xtick',1:2,'XTickLabel',xLegend,'fontsize',20);

ha(1) = subplot(2,2,2); %1 % 5
%plot([1;2],[nanmean(Data.pMiss.contra.pre); nanmean(Data.pMiss.contra.post)], 'x-','color',Color.Black ,'MarkerSize',10,'LineWidth',2,'markerfacecolor',[1 1 1 ]); hold on;
plot([1;2],[nanmean(Data.pHit.contra.pre); nanmean(Data.pHit.contra.post)], 'o-','color',Color.Contra ,'MarkerSize',10,'LineWidth',2,'markerfacecolor',[1 1 1 ]); hold on;
plot([1;2],[nanmean(Data.pFix.easy.contra.pre); nanmean(Data.pFix.easy.contra.post)], 'x-','color',Color.Black ,'MarkerSize',10,'LineWidth',2,'markerfacecolor',[1 1 1 ]); hold on;
plot([1;2],[nanmean(Data.pFA.easy.contra.pre); nanmean(Data.pFA.easy.contra.post)], 's-','color',Color.Contra ,'MarkerSize',10,'LineWidth',2,'markerfacecolor',[1 1 1 ]); hold on;

plot(2, nanmean(Data.pHit.contra.post), 'o','color',Color.Contra,'MarkerSize',10,'markerfacecolor',Color.Contra); hold on;
plot(2, nanmean(Data.pFix.easy.contra.post), 'X','color',Color.Black ,'MarkerSize',10,'markerfacecolor',Color.Black); hold on;
%plot(2, nanmean(Data.pMiss.contra.post), 'X','color',Color.Black ,'MarkerSize',10,'markerfacecolor',Color.Black); hold on;
plot(2, nanmean(Data.pFA.easy.contra.post), 's','color',Color.Contra ,'MarkerSize',10,'markerfacecolor',Color.Contra); hold on;
legend('tar','fix', 'dis', 'Location', 'NorthEast')


for i = 1: length(Data.pFA.easy.contra.pre)
    plot([0.9;2.1],[Data.pFix.easy.contra.pre(i); Data.pFix.easy.contra.post(i)], 'X','color',Color.Black ,'MarkerSize',8,'markerfacecolor',[1 1 1]); hold on;
  %  plot([0.9;2.1],[Data.pMiss.contra.pre(i); Data.pMiss.contra.post(i)], 'X','color',Color.Black ,'MarkerSize',8,'markerfacecolor',[1 1 1]); hold on;
    plot([0.9;2.1],[Data.pHit.contra.pre(i); Data.pHit.contra.post(i)], 'o','color',Color.Contra ,'MarkerSize',8,'markerfacecolor',[1 1 1]); hold on;
    plot([0.8;2.2],[Data.pFA.easy.contra.pre(i); Data.pFA.easy.contra.post(i)], 's','color',Color.Contra ,'MarkerSize',8,'markerfacecolor',[1 1 1]); hold on;
end

if Data.pFix.easy.contra.pvalue < 0.05
    ymax = nanmean([nanmean(Data.pFix.easy.contra.pre) ,nanmean(Data.pFix.easy.contra.post)]);
    ext_sigline([1,2], Data.pFix.easy.contra.pvalue,[],ymax +0.1,'x', Color.Black )
end

% if Data.pMiss.contra.pvalue < 0.05
%     ymax = nanmean([nanmean(Data.pMiss.contra.pre) ,nanmean(Data.pMiss.contra.post)]);
%     ext_sigline([1,2], Data.pMiss.contra.pvalue,[],ymax +0.1,'x', Color.Black )
% end
if Data.pHit.contra.pvalue < 0.05
    ymax = nanmean([nanmean(Data.pHit.contra.pre) ,nanmean(Data.pHit.contra.post)]);
    ext_sigline([1.4 1.6], Data.pHit.contra.pvalue,[],ymax,'x', Color.Contra )
end
if Data.pFA.easy.contra.pvalue < 0.05
    ymax = nanmean([nanmean(Data.pFA.easy.contra.pre) ,nanmean(Data.pFA.easy.contra.post)]);
    ext_sigline([1.4 1.6], Data.pFA.easy.contra.pvalue,[],ymax,'x', Color.Contra )
end

ylabel( 'Selection (easy)','fontsize',20,'fontweight','b', 'Interpreter', 'none' );
set(gca,'xlim',[0 3],'ylim',[0 1],'Xtick',1:2,'XTickLabel',xLegend,'fontsize',20);



ha(1) = subplot(2,2,3); %1 % 5
plot([1;2],[nanmean(Data.pFix.difficult.ipsi.pre); nanmean(Data.pFix.difficult.ipsi.post)], 'x-','color',Color.Black ,'MarkerSize',10,'LineWidth',2,'markerfacecolor',[1 1 1 ]); hold on;
%plot([1;2],[nanmean(Data.pMiss.ipsi.pre); nanmean(Data.pMiss.ipsi.post)], 'x-','color',Color.Black ,'MarkerSize',10,'LineWidth',2,'markerfacecolor',[1 1 1 ]); hold on;
plot([1;2],[nanmean(Data.pHit.ipsi.pre); nanmean(Data.pHit.ipsi.post)], 'o-','color',Color.Ipsi ,'MarkerSize',10,'LineWidth',2,'markerfacecolor',[1 1 1 ]); hold on;
plot([1;2],[nanmean(Data.pFA.difficult.ipsi.pre); nanmean(Data.pFA.difficult.ipsi.post)], 's-','color',Color.Ipsi ,'MarkerSize',10,'LineWidth',2,'markerfacecolor',[1 1 1 ]); hold on;

plot(2, nanmean(Data.pHit.ipsi.post), 'o','color',Color.Ipsi ,'MarkerSize',10,'markerfacecolor',Color.Ipsi); hold on;
plot(2, nanmean(Data.pFix.difficult.ipsi.post), 'X','color',Color.Black ,'MarkerSize',10,'markerfacecolor',Color.Black); hold on;
%plot(2, nanmean(Data.pMiss.ipsi.post), 'X','color',Color.Black ,'MarkerSize',10,'markerfacecolor',Color.Black); hold on;
plot(2, nanmean(Data.pFA.difficult.ipsi.post), 's','color',Color.Ipsi ,'MarkerSize',10,'markerfacecolor',Color.Ipsi); hold on;

%legend('tar','fix', 'dis', 'Location', 'SouthEast')
for i = 1: length(Data.pFA.difficult.ipsi.pre)
    plot([0.9;2.1],[Data.pHit.ipsi.pre(i); Data.pHit.ipsi.post(i)], 'o','color',Color.Ipsi ,'MarkerSize',8,'markerfacecolor',[1 1 1 ]); hold on;
    plot([0.9;2.1],[Data.pFix.difficult.ipsi.pre(i); Data.pFix.difficult.ipsi.post(i)], 'X','color',Color.Black,'MarkerSize',8,'markerfacecolor',[1 1 1 ]); hold on;
%    plot([0.9;2.1],[Data.pMiss.ipsi.pre(i); Data.pMiss.ipsi.post(i)], 'o','color',Color.Black,'MarkerSize',8,'markerfacecolor',[1 1 1 ]); hold on;
    plot([0.8;2.2],[Data.pFA.difficult.ipsi.pre(i); Data.pFA.difficult.ipsi.post(i)], 's','color',Color.Ipsi ,'MarkerSize',8,'markerfacecolor',[1 1 1 ]); hold on;
end

if Data.pFix.difficult.ipsi.pvalue < 0.05
    ymax = nanmean([nanmean(Data.pFix.difficult.ipsi.pre) ,nanmean(Data.pFix.difficult.ipsi.post)]);
    ext_sigline([1,2], Data.pFix.difficult.ipsi.pvalue,[],ymax +0.1,'x', Color.Black )
end

% if Data.pMiss.ipsi.pvalue < 0.05
%     ymax = nanmean([nanmean(Data.pMiss.ipsi.pre) ,nanmean(Data.pMiss.ipsi.post)]);
%     ext_sigline([1,2], Data.pMiss.ipsi.pvalue,[],ymax +0.1,'x', Color.Black )
% end

if Data.pHit.ipsi.pvalue < 0.05
    ymax = nanmean([nanmean(Data.pHit.ipsi.pre) ,nanmean(Data.pHit.ipsi.post)]);
    ext_sigline([1.4 1.6], Data.pHit.ipsi.pvalue,[],ymax ,'x', Color.Ipsi )
end

if Data.pFA.difficult.ipsi.pvalue < 0.05
    ymax = nanmean([nanmean(Data.pFA.difficult.ipsi.pre) ,nanmean(Data.pFA.difficult.ipsi.post)]);
    ext_sigline([1.4 1.6], Data.pFA.difficult.ipsi.pvalue,[],ymax,'x', Color.Ipsi )
end
ylabel( 'Selection (difficult)','fontsize',20,'fontweight','b', 'Interpreter', 'none' );
set(gca,'xlim',[0 3],'ylim',[0 1],'Xtick',1:2,'XTickLabel',xLegend,'fontsize',20);

ha(1) = subplot(2,2,4); %1 % 5
plot([1;2],[nanmean(Data.pHit.contra.pre); nanmean(Data.pHit.contra.post)], 'o-','color',Color.Contra ,'MarkerSize',10,'LineWidth',2,'markerfacecolor',[1 1 1 ]); hold on;
plot([1;2],[nanmean(Data.pFix.difficult.contra.pre); nanmean(Data.pFix.difficult.contra.post)], 'x-','color',Color.Black ,'MarkerSize',10,'LineWidth',2,'markerfacecolor',[1 1 1 ]); hold on;
%plot([1;2],[nanmean(Data.pMiss.contra.pre); nanmean(Data.pMiss.contra.post)], 'x-','color',Color.Black ,'MarkerSize',10,'LineWidth',2,'markerfacecolor',[1 1 1 ]); hold on;
plot([1;2],[nanmean(Data.pFA.difficult.contra.pre); nanmean(Data.pFA.difficult.contra.post)], 's-','color',Color.Contra ,'MarkerSize',10,'LineWidth',2,'markerfacecolor',[1 1 1 ]); hold on;
plot(2, nanmean(Data.pHit.contra.post), 'o','color',Color.Contra,'MarkerSize',10,'markerfacecolor',Color.Contra); hold on;
plot(2, nanmean(Data.pFix.difficult.contra.post), 'X','color',Color.Black ,'MarkerSize',10,'markerfacecolor',Color.Black); hold on;
%plot(2, nanmean(Data.pMiss.contra.post), 'X','color',Color.Black ,'MarkerSize',10,'markerfacecolor',Color.Black); hold on;
plot(2, nanmean(Data.pFA.difficult.contra.post), 's','color',Color.Contra ,'MarkerSize',10,'markerfacecolor',Color.Contra); hold on;
legend('tar','fix', 'dis', 'Location', 'NorthEast')


for i = 1: length(Data.pFA.difficult.contra.pre)
    plot([0.9;2.1],[Data.pFix.difficult.contra.pre(i); Data.pFix.difficult.contra.post(i)], 'X','color',Color.Black ,'MarkerSize',8,'markerfacecolor',[1 1 1]); hold on;
 %   plot([0.9;2.1],[Data.pMiss.contra.pre(i); Data.pMiss.contra.post(i)], 'X','color',Color.Black ,'MarkerSize',8,'markerfacecolor',[1 1 1]); hold on;
    plot([0.9;2.1],[Data.pHit.contra.pre(i); Data.pHit.contra.post(i)], 'o','color',Color.Contra ,'MarkerSize',8,'markerfacecolor',[1 1 1]); hold on;
    plot([0.8;2.2],[Data.pFA.difficult.contra.pre(i); Data.pFA.difficult.contra.post(i)], 's','color',Color.Contra ,'MarkerSize',8,'markerfacecolor',[1 1 1]); hold on;
end

if Data.pFix.difficult.contra.pvalue < 0.05
    ymax = nanmean([nanmean(Data.pFix.difficult.contra.pre) ,nanmean(Data.pFix.difficult.contra.post)]);
    ext_sigline([1,2], Data.pFix.difficult.contra.pvalue,[],ymax +0.1,'x', Color.Black )
end


% if Data.pMiss.contra.pvalue < 0.05
%     ymax = nanmean([nanmean(Data.pMiss.contra.pre) ,nanmean(Data.pMiss.contra.post)]);
%     ext_sigline([1,2], Data.pMiss.contra.pvalue,[],ymax +0.1,'x', Color.Black )
% end

if Data.pHit.contra.pvalue < 0.05
    ymax = nanmean([nanmean(Data.pHit.contra.pre) ,nanmean(Data.pHit.contra.post)]);
    ext_sigline([1.4 1.6], Data.pHit.contra.pvalue,[],ymax,'x', Color.Contra )
end
if Data.pFA.difficult.contra.pvalue < 0.05
    ymax = nanmean([nanmean(Data.pFA.difficult.contra.pre) ,nanmean(Data.pFA.difficult.contra.post)]);
    ext_sigline([1.4 1.6], Data.pFA.difficult.contra.pvalue,[],ymax,'x', Color.Contra)
end

ylabel( 'Selection (difficult)','fontsize',20,'fontweight','b', 'Interpreter', 'none' );
set(gca,'xlim',[0 3],'ylim',[0 1],'Xtick',1:2,'XTickLabel',xLegend,'fontsize',20);


if SaveGraph
    h = figure(1);
    % savefig(h, [path_SaveFig ,'fig', filesep,  'TarDistStimuli_HR_FAR_' monkey,'_' experiment ,'_' folder ,'.fig'])
    print(h,[path_SaveFig ,filesep,  ['png', ControlFolder] ,filesep,  'SingleStimuli_SELECTION' monkey,'_' experiment ,'_' folder, '.png'], '-dpng')
    set(h,'Renderer','Painters');
    set(h,'PaperPositionMode','auto')
    compl_filename =  [path_SaveFig , ['ai', ControlFolder]  ,filesep,  'SingleStimuli_SELECTION' monkey,'_' experiment ,'_', folder ,'.ai'] ;
    print(h,'-depsc',compl_filename);
    close all;
end
