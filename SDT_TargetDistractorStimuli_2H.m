function SDT_TargetDistractorStimuli_2H(monkey, folder,dataset,folder_baseline,dataset_baseline,experiment,path_SaveFig,path_SaveData,path_getData, ResponseBias_DisplayColor)
close all;

SaveGraph =1; 
SaveTable =1; 

roundValue = 3; 
NonParametri = 0; 
%% Colors For Plotting
Color.Ipsi      = [0 0 1 ];
Color.Contra    = [1 0 1 ];
Color.Black     = [0 0 0 ];
% Size of Writing
fs = 28; % font size
MarkSize_GraphFAR_HR = 12;
xlim_SDT_easy = [1 5];
xlim_SDT_diff = [0 4];

MarkSize_GraphFAR_HR_PerSession = 7;
MarkSize_GraphFAR_HR = 17;
MarkSize_CritDpr_small = 8;
MarkSize_CritDpr = 20;
LineWidthSize = 3;

%% graph
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
    load([path_getData,  monkey ,'\behavior\' folder, filesep, Files(end).name ])
    disp( ['analyzed dataset:  ',  filesep,folder, filesep, Files(end).name ])
    
    SET.TimeLine        = {'stim_off', 'go_signal'};
    if isempty(folder_baseline)== 0
        Bline =load([path_getData, monkey ,'\behavior\' folder_baseline, filesep, dataset_baseline ]);
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
maxSess =  length(data_struct_per_session.single_targ.ipsi.stim_off);
positionExp     = {'ipsi', 'contra'};
TimeLineExp         = {'pre', 'post'};
Difficulty          = { 'difficult','easy',};

%% the datastructure is that position "ipsi" -> Target is on the ipsilateral side and in the same category is the distractor "contra"
% Miss and CR are the same ... Fixation
for indSess = 1:   maxSess
    for indDiff = 1:length(Difficulty)
        for indPos = 1:length(positionExp)
            
            for indTime = 1:length(SET.TimeLine)
                if   strcmp(positionExp{indPos}, 'contra')
                    Session.(Difficulty{indDiff}).(positionExp{indPos}).(TimeLineExp{indTime})(indSess)               = data_struct_per_session.targ_distr_2HF.(positionExp{indPos}).(SET.TimeLine{indTime}){indDiff,indSess}.session;
                    
                    NumTotalTargets.(Difficulty{indDiff}).(positionExp{indPos}).(TimeLineExp{indTime})(indSess)               = data_struct_per_session.targ_distr_2HF.(positionExp{indPos}).(SET.TimeLine{indTime}){indDiff,indSess}.num_total_trials;
                    %Target is on the contralateral side
                    Hit.(Difficulty{indDiff}).(positionExp{indPos}).(TimeLineExp{indTime})(indSess)                           = data_struct_per_session.targ_distr_2HF.(positionExp{indPos}).(SET.TimeLine{indTime}){indDiff,indSess}.num_choice_contra;
                    Miss.(Difficulty{indDiff}).(positionExp{indPos}).(TimeLineExp{indTime})(indSess)                          = data_struct_per_session.targ_distr_2HF.(positionExp{indPos}).(SET.TimeLine{indTime}){indDiff,indSess}.num_choice_fixation;
                    % FA is also on the contralateral side
                    NumTotalDistractor.(Difficulty{indDiff}).(positionExp{indPos}).(TimeLineExp{indTime})(indSess)       = data_struct_per_session.targ_distr_2HF.ipsi.(SET.TimeLine{indTime}){indDiff,indSess}.num_total_trials;
                    CR.(Difficulty{indDiff}).(positionExp{indPos}).(TimeLineExp{indTime})(indSess)                       = data_struct_per_session.targ_distr_2HF.ipsi.(SET.TimeLine{indTime}){indDiff,indSess}.num_choice_fixation;
                    FA.(Difficulty{indDiff}).(positionExp{indPos}).(TimeLineExp{indTime})(indSess)                       = data_struct_per_session.targ_distr_2HF.ipsi.(SET.TimeLine{indTime}){indDiff,indSess}.num_choice_contra;
                elseif     strcmp(positionExp{indPos}, 'ipsi')
                    Session.(Difficulty{indDiff}).(positionExp{indPos}).(TimeLineExp{indTime})(indSess)               = data_struct_per_session.targ_distr_2HF.(positionExp{indPos}).(SET.TimeLine{indTime}){indDiff,indSess}.session;
                    NumTotalTargets.(Difficulty{indDiff}).(positionExp{indPos}).(TimeLineExp{indTime})(indSess)               = data_struct_per_session.targ_distr_2HF.(positionExp{indPos}).(SET.TimeLine{indTime}){indDiff,indSess}.num_total_trials;
                    Hit.(Difficulty{indDiff}).(positionExp{indPos}).(TimeLineExp{indTime})(indSess)                           = data_struct_per_session.targ_distr_2HF.(positionExp{indPos}).(SET.TimeLine{indTime}){indDiff,indSess}.num_choice_ipsi;
                    Miss.(Difficulty{indDiff}).(positionExp{indPos}).(TimeLineExp{indTime})(indSess)                          = data_struct_per_session.targ_distr_2HF.(positionExp{indPos}).(SET.TimeLine{indTime}){indDiff,indSess}.num_choice_fixation;
                    NumTotalDistractor.(Difficulty{indDiff}).(positionExp{indPos}).(TimeLineExp{indTime})(indSess)       = data_struct_per_session.targ_distr_2HF.contra.(SET.TimeLine{indTime}){indDiff,indSess}.num_total_trials;
                    CR.(Difficulty{indDiff}).(positionExp{indPos}).(TimeLineExp{indTime})(indSess)                       = data_struct_per_session.targ_distr_2HF.contra.(SET.TimeLine{indTime}){indDiff,indSess}.num_choice_fixation;
                    FA.(Difficulty{indDiff}).(positionExp{indPos}).(TimeLineExp{indTime})(indSess)                       = data_struct_per_session.targ_distr_2HF.contra.(SET.TimeLine{indTime}){indDiff,indSess}.num_choice_ipsi;
                end
            end
        end
    end
end

%% ONLY FOR INACTIVATION: extract the Data for the Baseline-session
if strcmp(experiment, 'Inactivation') && isempty(folder_baseline) == 0
    maxSess =  length(Bline.data_struct_per_session.single_targ.ipsi.stim_off);
    
    for indSess = 1:   maxSess
        for indDiff = 1:length(Difficulty)
            for indPos = 1:length(positionExp)
                for indTime = 2
                    if   strcmp(positionExp{indPos}, 'contra')
                        Ctr_Session.(Difficulty{indDiff}).(positionExp{indPos}).(TimeLineExp{1})(indSess)               = Bline.data_struct_per_session.targ_distr_2HF.combined.(SET.TimeLine{indTime}){indDiff,indSess}.session;
                        NumTotalTargets.(Difficulty{indDiff}).(positionExp{indPos}).(TimeLineExp{1})(indSess)               = Bline.data_struct_per_session.targ_distr_2HF.(positionExp{indPos}).(SET.TimeLine{indTime}){indDiff,indSess}.num_total_trials;
                        Hit.(Difficulty{indDiff}).(positionExp{indPos}).(TimeLineExp{1})(indSess)                           = Bline.data_struct_per_session.targ_distr_2HF.(positionExp{indPos}).(SET.TimeLine{indTime}){indDiff,indSess}.num_choice_contra;
                        Miss.(Difficulty{indDiff}).(positionExp{indPos}).(TimeLineExp{1})(indSess)                          = Bline.data_struct_per_session.targ_distr_2HF.(positionExp{indPos}).(SET.TimeLine{indTime}){indDiff,indSess}.num_choice_fixation;
                        % logical about ipsi
                        NumTotalDistractor.(Difficulty{indDiff}).(positionExp{indPos}).(TimeLineExp{1})(indSess)       = Bline.data_struct_per_session.targ_distr_2HF.ipsi.(SET.TimeLine{indTime}){indDiff,indSess}.num_total_trials;
                        CR.(Difficulty{indDiff}).(positionExp{indPos}).(TimeLineExp{1})(indSess)                       = Bline.data_struct_per_session.targ_distr_2HF.ipsi.(SET.TimeLine{indTime}){indDiff,indSess}.num_choice_fixation;
                        FA.(Difficulty{indDiff}).(positionExp{indPos}).(TimeLineExp{1})(indSess)                       = Bline.data_struct_per_session.targ_distr_2HF.ipsi.(SET.TimeLine{indTime}){indDiff,indSess}.num_choice_contra;
                        
                    elseif     strcmp(positionExp{indPos}, 'ipsi')
                        Ctr_Session.(Difficulty{indDiff}).(positionExp{indPos}).(TimeLineExp{1})(indSess)               = Bline.data_struct_per_session.targ_distr_2HF.(positionExp{indPos}).(SET.TimeLine{indTime}){indDiff,indSess}.session;
                        NumTotalTargets.(Difficulty{indDiff}).(positionExp{indPos}).(TimeLineExp{1})(indSess)               = Bline.data_struct_per_session.targ_distr_2HF.(positionExp{indPos}).(SET.TimeLine{indTime}){indDiff,indSess}.num_total_trials;
                        Hit.(Difficulty{indDiff}).(positionExp{indPos}).(TimeLineExp{1})(indSess)                           = Bline.data_struct_per_session.targ_distr_2HF.(positionExp{indPos}).(SET.TimeLine{indTime}){indDiff,indSess}.num_choice_ipsi;
                        Miss.(Difficulty{indDiff}).(positionExp{indPos}).(TimeLineExp{1})(indSess)                          = Bline.data_struct_per_session.targ_distr_2HF.(positionExp{indPos}).(SET.TimeLine{indTime}){indDiff,indSess}.num_choice_fixation;
                        NumTotalDistractor.(Difficulty{indDiff}).(positionExp{indPos}).(TimeLineExp{1})(indSess)       = Bline.data_struct_per_session.targ_distr_2HF.contra.(SET.TimeLine{indTime}){indDiff,indSess}.num_total_trials;
                        CR.(Difficulty{indDiff}).(positionExp{indPos}).(TimeLineExp{1})(indSess)                       = Bline.data_struct_per_session.targ_distr_2HF.contra.(SET.TimeLine{indTime}){indDiff,indSess}.num_choice_fixation;
                        FA.(Difficulty{indDiff}).(positionExp{indPos}).(TimeLineExp{1})(indSess)                       = Bline.data_struct_per_session.targ_distr_2HF.contra.(SET.TimeLine{indTime}){indDiff,indSess}.num_choice_ipsi;
                        
                    end
                end
            end
        end
    end
    
end
%%  avoid 0 or Inf probabilities
%use the approach regardless wheather or not extreme values are obtained
for indDiff = 1:length(Difficulty)
    for indPos = 1:length(positionExp) %contra vs ipsi selection
        for indTime = 1:length(SET.TimeLine) %pre vs post
            Hit.(Difficulty{indDiff}).(positionExp{indPos}).(TimeLineExp{indTime})            = Hit.(Difficulty{indDiff}).(positionExp{indPos}).(TimeLineExp{indTime}) + 0.5 ;
            Miss.(Difficulty{indDiff}).(positionExp{indPos}).(TimeLineExp{indTime})           = Miss.(Difficulty{indDiff}).(positionExp{indPos}).(TimeLineExp{indTime}) + 0.5 ;
            FA.(Difficulty{indDiff}).(positionExp{indPos}).(TimeLineExp{indTime})        = FA.(Difficulty{indDiff}).(positionExp{indPos}).(TimeLineExp{indTime}) + 0.5 ;
            CR.(Difficulty{indDiff}).(positionExp{indPos}).(TimeLineExp{indTime})        = CR.(Difficulty{indDiff}).(positionExp{indPos}).(TimeLineExp{indTime})  + 0.5;
        end
    end
end
%% calculate d'prime & criterion
% for ipsi vs contra
% control vs inactivation
% difficult vs easy
for indPos = 1:length(positionExp) %contra vs ipsi selection
    for indTime = 1:length(SET.TimeLine) %pre vs post
        %target is presented  for example the contralateral side % the
        %FA is one the ipsiversive side
        if indPos == 1 ; indPos2 = 2 ; else indPos2 = 1 ; end
        
        pHit.easy.(positionExp{indPos}).(TimeLineExp{indTime})=    Hit.easy.(positionExp{indPos}).(TimeLineExp{indTime}) ./ ...
            sum(Hit.easy.(positionExp{indPos}).(TimeLineExp{indTime})  + Miss.easy.(positionExp{indPos}).(TimeLineExp{indTime})+ FA.easy.(positionExp{indPos2}).(TimeLineExp{indTime}),1);
        
        pHit.difficult.(positionExp{indPos}).(TimeLineExp{indTime})=    Hit.difficult.(positionExp{indPos}).(TimeLineExp{indTime}) ./ ...
            sum(Hit.difficult.(positionExp{indPos}).(TimeLineExp{indTime})  + Miss.difficult.(positionExp{indPos}).(TimeLineExp{indTime})+ FA.difficult.(positionExp{indPos2}).(TimeLineExp{indTime}),1);
        
        %distractor is presented
        pFA.easy.(positionExp{indPos}).(TimeLineExp{indTime}) = FA.easy.(positionExp{indPos}).(TimeLineExp{indTime}) ./ ...
            sum(FA.easy.(positionExp{indPos}).(TimeLineExp{indTime})+  CR.easy.(positionExp{indPos}).(TimeLineExp{indTime})+Hit.easy.(positionExp{indPos2}).(TimeLineExp{indTime}),1);
        
        pFA.difficult.(positionExp{indPos}).(TimeLineExp{indTime}) = FA.difficult.(positionExp{indPos}).(TimeLineExp{indTime}) ./ ...
            sum(FA.difficult.(positionExp{indPos}).(TimeLineExp{indTime})+  CR.difficult.(positionExp{indPos}).(TimeLineExp{indTime})+ Hit.difficult.(positionExp{indPos2}).(TimeLineExp{indTime}),1);
        
        
        
        %distractor is presented
        pMiss.easy.(positionExp{indPos}).(TimeLineExp{indTime}) = Miss.easy.(positionExp{indPos}).(TimeLineExp{indTime}) ./ ...
            sum(Hit.easy.(positionExp{indPos}).(TimeLineExp{indTime})  + Miss.easy.(positionExp{indPos}).(TimeLineExp{indTime})+ FA.easy.(positionExp{indPos2}).(TimeLineExp{indTime}),1);
        
        pMiss.difficult.(positionExp{indPos}).(TimeLineExp{indTime}) = Miss.difficult.(positionExp{indPos}).(TimeLineExp{indTime}) ./ ...
            sum(Hit.difficult.(positionExp{indPos}).(TimeLineExp{indTime})  + Miss.difficult.(positionExp{indPos}).(TimeLineExp{indTime})+ FA.difficult.(positionExp{indPos2}).(TimeLineExp{indTime}),1);
        
        %SDT-variables calculated spatially Hitrate (contra) & False
        %Alarm rate (contra)
        [d_prime.easy.(positionExp{indPos}).(TimeLineExp{indTime}),beta.easy.(positionExp{indPos}).(TimeLineExp{indTime}),criterion.easy.(positionExp{indPos}).(TimeLineExp{indTime})] = ...
            testsim_dprime( pHit.easy.(positionExp{indPos}).(TimeLineExp{indTime}), pFA.easy.(positionExp{indPos}).(TimeLineExp{indTime}));
        [d_prime.difficult.(positionExp{indPos}).(TimeLineExp{indTime}),beta.difficult.(positionExp{indPos}).(TimeLineExp{indTime}),criterion.difficult.(positionExp{indPos}).(TimeLineExp{indTime})] = ...
            testsim_dprime( pHit.difficult.(positionExp{indPos}).(TimeLineExp{indTime}), pFA.difficult.(positionExp{indPos}).(TimeLineExp{indTime}));
        %SDT-variables calculated as perception Hitrate (contra) & False
        %Alarm rate (ipsi)
        
        
    end
end

% For Not spatial comparison
% D-T presented which are now used to calculate dprime
TimeLineExp_directCmp         = {'pre_direct', 'post_direct'};
for indPos = 1:length(positionExp) %contra vs ipsi selection
    for indTime = 1:length(SET.TimeLine) %pre vs post
        %target is presented  for example the contralateral side % the
        %FA is one the ipsiversive side
        if indPos == 1 ; indPos2 = 2 ; else indPos2 = 1 ; end
        [d_prime.easy.(positionExp{indPos}).(TimeLineExp_directCmp{indTime}),beta.easy.(positionExp{indPos}).(TimeLineExp_directCmp{indTime}),criterion.easy.(positionExp{indPos}).(TimeLineExp_directCmp{indTime})] = ...
            testsim_dprime( pHit.easy.(positionExp{indPos}).(TimeLineExp{indTime}), pFA.easy.(positionExp{indPos2}).(TimeLineExp{indTime}));
        [d_prime.difficult.(positionExp{indPos}).(TimeLineExp_directCmp{indTime}),beta.difficult.(positionExp{indPos}).(TimeLineExp_directCmp{indTime}),criterion.difficult.(positionExp{indPos}).(TimeLineExp_directCmp{indTime})] = ...
            testsim_dprime( pHit.difficult.(positionExp{indPos}).(TimeLineExp{indTime}), pFA.difficult.(positionExp{indPos2}).(TimeLineExp{indTime}));
    end
end
%% Statistic - Contrast: Pre vs Post 
TablePwerte = [];

if strcmp(experiment, 'Inactivation') && isempty(folder_baseline) == 0 && NonParametri == 0
 disp(['independent t-test (ttest2) is displayed for the',  experiment , monkey, 'Target_Distractor trials'] )
    for indPos = 1:length(positionExp) %contra vs ipsi selection
        [h,p,ci,stats]  = ttest2( d_prime.easy.(positionExp{indPos}).post, d_prime.easy.(positionExp{indPos}).pre);
        d_prime.easy.(positionExp{indPos}).pvalue = p;
        d_prime.easy.(positionExp{indPos}).tstat = stats.tstat;    
        disp(['dprime ', 'easy ' , (positionExp{indPos}) '_',num2str(round(p,2)),'_t=', num2str(round(d_prime.easy.(positionExp{indPos}).tstat,2))] )
 
        
            TabPwerte = [ ];
            TabPwerte = table({monkey}, {'DoubleDiffStim'}, {'dprime'},{ 'easy'} , (positionExp(indPos)) ,round(p,roundValue),{'t'} , round(d_prime.easy.(positionExp{indPos}).tstat,2), {'indepTtest '} );
            TablePwerte = [ TablePwerte; TabPwerte ];
            
            [h,p,ci,stats]  = ttest2( criterion.easy.(positionExp{indPos}).post, criterion.easy.(positionExp{indPos}).pre);
        criterion.easy.(positionExp{indPos}).pvalue = p;
        criterion.easy.(positionExp{indPos}).tstat = stats.tstat;    
        disp(['criterion ', 'easy ' , (positionExp{indPos}) '_',num2str(round(p,2)),'_t=', num2str(round(criterion.easy.(positionExp{indPos}).tstat,2))] )
             TabPwerte = [ ];
            TabPwerte = table({monkey}, {'DoubleDiffStim'}, {'criterion'},{ 'easy'} , (positionExp(indPos)) ,round(p,roundValue),{'t'} , round(criterion.easy.(positionExp{indPos}).tstat,2), {'indepTtest '} );
            TablePwerte = [ TablePwerte; TabPwerte ];
            
            [h,p,ci,stats]  = ttest2( d_prime.difficult.(positionExp{indPos}).post, d_prime.difficult.(positionExp{indPos}).pre);
        d_prime.difficult.(positionExp{indPos}).pvalue = p;
        d_prime.difficult.(positionExp{indPos}).tstat = stats.tstat;    
        disp(['d_prime ', 'difficult ' , (positionExp{indPos}) '_',num2str(round(p,2)),'_t=', num2str(round(d_prime.difficult.(positionExp{indPos}).tstat,2))] )
  
                    TabPwerte = [ ];
            TabPwerte = table({monkey}, {'DoubleDiffStim'}, {'dprime'},{ 'difficult'} , (positionExp(indPos)) ,round(p,roundValue),{'t'} , round(d_prime.difficult.(positionExp{indPos}).tstat,2), {'indepTtest '} );
            TablePwerte = [ TablePwerte; TabPwerte ];
            
            [h,p,ci,stats]  = ttest2( criterion.difficult.(positionExp{indPos}).post, criterion.difficult.(positionExp{indPos}).pre);
        criterion.difficult.(positionExp{indPos}).pvalue = p;
        criterion.difficult.(positionExp{indPos}).tstat = stats.tstat;    
        disp(['criterion ', 'difficult ' , (positionExp{indPos}) '_',num2str(round(p,2)),'_t=', num2str(round(criterion.difficult.(positionExp{indPos}).tstat,2))] )

                    TabPwerte = [ ];
            TabPwerte = table({monkey}, {'DoubleDiffStim'}, {'criterion'},{ 'difficult'} , (positionExp(indPos)) ,round(p,roundValue),{'t'} , round(criterion.difficult.(positionExp{indPos}).tstat,2), {'indepTtest '} );
            TablePwerte = [ TablePwerte; TabPwerte ];
            
        [h,p,ci,stat] = ttest2( pFA.easy.(positionExp{indPos}).post, pFA.easy.(positionExp{indPos}).pre);
        pFA.easy.(positionExp{indPos}).pvalue = p;
        pFA.easy.(positionExp{indPos}).tstat = stat.tstat;
        pFA.easy.(positionExp{indPos}).r = norminv(p)/sqrt(length(pFA.easy.(positionExp{indPos}).post) + length(pFA.easy.(positionExp{indPos}).pre));
        disp(['pFA', 'easy ' , (positionExp{indPos}) '_',num2str(round(p,2)),'_t=', num2str(round(pFA.easy.(positionExp{indPos}).tstat,2))] )
        
        
        [h,p,ci,stat] = ttest2( pHit.easy.(positionExp{indPos}).post, pHit.easy.(positionExp{indPos}).pre);
        pHit.easy.(positionExp{indPos}).pvalue = p;
        pHit.easy.(positionExp{indPos}).tstat = stat.tstat;
        pHit.easy.(positionExp{indPos}).r = norminv(p)/sqrt(length(criterion.difficult.(positionExp{indPos}).post) + length( pHit.easy.(positionExp{indPos}).pre));
        disp(['pHit ', 'easy ' , (positionExp{indPos}) '_',num2str(round(p,2)),'_t=', num2str(round( pHit.easy.(positionExp{indPos}).tstat,2))] )
        
        [h,p,stat] = ttest2( pMiss.easy.(positionExp{indPos}).post, pMiss.easy.(positionExp{indPos}).pre);
        pMiss.easy.(positionExp{indPos}).pvalue = p;
        
        [h,p,ci,stat] = ttest2( pFA.difficult.(positionExp{indPos}).post, pFA.difficult.(positionExp{indPos}).pre);
        pFA.difficult.(positionExp{indPos}).pvalue = p;
        pFA.difficult.(positionExp{indPos}).tstat = stat.tstat;
        pFA.difficult.(positionExp{indPos}).r = norminv(p)/sqrt(length(criterion.difficult.(positionExp{indPos}).post) + length(pFA.difficult.(positionExp{indPos}).pre));
        disp(['pFA ', 'difficult ' , (positionExp{indPos}) '_',num2str(round(p,2)),'_t=', num2str(round(pFA.difficult.(positionExp{indPos}).tstat,2))] )
        
        [h,p,ci,stat] = ttest2( pHit.difficult.(positionExp{indPos}).post, pHit.difficult.(positionExp{indPos}).pre);
        pHit.difficult.(positionExp{indPos}).pvalue = p;
        pHit.difficult.(positionExp{indPos}).tstat = stat.tstat;
        pHit.difficult.(positionExp{indPos}).r = norminv(p)/sqrt(length(pHit.difficult.(positionExp{indPos}).post) + length(criterion.difficult.(positionExp{indPos}).pre));
        disp(['pHit ', 'difficult ' , (positionExp{indPos}) '_',num2str(round(p,2)),'_t=', num2str(round(pHit.difficult.(positionExp{indPos}).tstat,2))] )
        
        [h,p,ci,stat] = ttest2( pMiss.difficult.(positionExp{indPos}).post, pMiss.difficult.(positionExp{indPos}).pre);
        pMiss.difficult.(positionExp{indPos}).pvalue = p;
        
        %Direct Comparison
        [p,h,stat] = ranksum(d_prime.easy.(positionExp{indPos}).pre_direct, d_prime.easy.(positionExp{indPos}).post_direct);
        d_prime.easy.(positionExp{indPos}).pvalue_direct = p;
        [p,h,stat] = ranksum(criterion.easy.(positionExp{indPos}).pre_direct, criterion.easy.(positionExp{indPos}).post_direct);
        criterion.easy.(positionExp{indPos}).pvalue_direct = p;
        [p,h,stat] = ranksum(d_prime.difficult.(positionExp{indPos}).pre_direct, d_prime.difficult.(positionExp{indPos}).post_direct);
        d_prime.difficult.(positionExp{indPos}).pvalue_direct = p;
        [p,h,stat] = ranksum(criterion.difficult.(positionExp{indPos}).pre_direct, criterion.difficult.(positionExp{indPos}).post_direct);
        criterion.difficult.(positionExp{indPos}).pvalue_direct = p;
        
    end
else
    for indPos = 1:length(positionExp) %contra vs ipsi selection
         disp(['non-parametric test (ranksum, independent) is displayed for the',  'experiment'] )

        [p,h,stat] = ranksum(d_prime.easy.(positionExp{indPos}).pre, d_prime.easy.(positionExp{indPos}).post);
        d_prime.easy.(positionExp{indPos}).pvalue = p;
        
        TabPwerte = [ ];
        TabPwerte = table({monkey }, {'DoubleDiffStim'}, {'dprime'},{ 'easy'} , (positionExp(indPos)) ,round(p,roundValue),{'U'} , round(stat.ranksum,2), {'ranksum '} );
        TablePwerte = [ TablePwerte; TabPwerte ];
        
        
        [p,h,stat] = ranksum(criterion.easy.(positionExp{indPos}).pre, criterion.easy.(positionExp{indPos}).post);
        criterion.easy.(positionExp{indPos}).pvalue = p;
        TabPwerte = [ ];
        TabPwerte = table({monkey }, {'DoubleDiffStim'}, {'criterion'},{ 'easy'} , (positionExp(indPos)) ,round(p,roundValue),{'U'} , round(stat.ranksum,2), {'ranksum '} );
        TablePwerte = [ TablePwerte; TabPwerte ];
        
        
        [p,h,stat] = ranksum(d_prime.difficult.(positionExp{indPos}).pre, d_prime.difficult.(positionExp{indPos}).post);
        d_prime.difficult.(positionExp{indPos}).pvalue = p;
        TabPwerte = [ ];
        TabPwerte = table({monkey }, {'DoubleDiffStim'}, {'dprime'},{ 'difficult'} , (positionExp(indPos)) ,round(p,roundValue),{'U'} , round(stat.ranksum,2), {'ranksum '} );
        TablePwerte = [ TablePwerte; TabPwerte ];
        
        
        [p,h,stat] = ranksum(criterion.difficult.(positionExp{indPos}).pre, criterion.difficult.(positionExp{indPos}).post);
        criterion.difficult.(positionExp{indPos}).pvalue = p;
        TabPwerte = [ ];
        TabPwerte = table({monkey }, {'DoubleDiffStim'}, {'criterion'},{ 'difficult'} , (positionExp(indPos)) ,round(p,roundValue),{'U'} , round(stat.ranksum,2), {'ranksum '} );
        TablePwerte = [ TablePwerte; TabPwerte ];
        
        
        [p,h,stat] = ranksum(d_prime.easy.(positionExp{indPos}).pre_direct, d_prime.easy.(positionExp{indPos}).post_direct);
        d_prime.easy.(positionExp{indPos}).pvalue_direct = p;
        TabPwerte = [ ];
        
        
        [p,h,stat] = ranksum(criterion.easy.(positionExp{indPos}).pre_direct, criterion.easy.(positionExp{indPos}).post_direct);
        criterion.easy.(positionExp{indPos}).pvalue_direct = p;
        
        
        [p,h,stat] = ranksum(d_prime.difficult.(positionExp{indPos}).pre_direct, d_prime.difficult.(positionExp{indPos}).post_direct);
        d_prime.difficult.(positionExp{indPos}).pvalue_direct = p;
        
        [p,h,stat] = ranksum(criterion.difficult.(positionExp{indPos}).pre_direct, criterion.difficult.(positionExp{indPos}).post_direct);
        criterion.difficult.(positionExp{indPos}).pvalue_direct = p;
        
        
        
        [H, p,CI,STATS] = ttest(pFA.easy.(positionExp{indPos}).pre, pFA.easy.(positionExp{indPos}).post);
        pFA.easy.(positionExp{indPos}).pvalue = p;
        TabPwerte = [ ];
        TabPwerte = table({monkey}, {'DoubleDiffStim'}, {'pFA'},{ 'easy'} , (positionExp(indPos)) ,round(p,roundValue),{'t'} , round(STATS.tstat,3), {'depTtest'} );
        TablePwerte = [ TablePwerte; TabPwerte ];
        
        [H, p,CI,STATS] = ttest(pHit.easy.(positionExp{indPos}).pre, pHit.easy.(positionExp{indPos}).post);
        pHit.easy.(positionExp{indPos}).pvalue = p;
        TabPwerte = [ ];
        TabPwerte = table({monkey}, {'DoubleDiffStim'}, {'pHit'},{ 'easy'} , (positionExp(indPos)) ,round(p,roundValue),{'t'} , round(STATS.tstat,3), {'depTtest'} );
        TablePwerte = [ TablePwerte; TabPwerte ];
        
        [H, p,CI,STATS]= ttest2(pMiss.easy.(positionExp{indPos}).pre, pMiss.easy.(positionExp{indPos}).post);
        pMiss.easy.(positionExp{indPos}).pvalue = p;
        
        [H, p,CI,STATS]= ttest(pFA.difficult.(positionExp{indPos}).pre, pFA.difficult.(positionExp{indPos}).post);
        pFA.difficult.(positionExp{indPos}).pvalue = p;
        TabPwerte = [ ];
        TabPwerte = table({monkey}, {'DoubleDiffStim'}, {'pFA'},{ 'difficult'} , (positionExp(indPos)) ,round(p,roundValue),{'t'} , round(STATS.tstat,3), {'depTtest'} );
        TablePwerte = [ TablePwerte; TabPwerte ];
        
        [H, p,CI,STATS] = ttest(pHit.difficult.(positionExp{indPos}).pre, pHit.difficult.(positionExp{indPos}).post);
        pHit.difficult.(positionExp{indPos}).pvalue = p;
        TabPwerte = [ ];
        TabPwerte = table({monkey}, {'DoubleDiffStim'}, {'pHit'},{ 'difficult'} , (positionExp(indPos)) ,round(p,roundValue),{'t'} , round(STATS.tstat,3), {'depTtest'} );
        TablePwerte = [ TablePwerte; TabPwerte ];
        
        [H, p,CI,STATS]= ttest(pMiss.difficult.(positionExp{indPos}).pre, pMiss.difficult.(positionExp{indPos}).post);
        pMiss.difficult.(positionExp{indPos}).pvalue = p;
    end
end

%addtoDropbox = 'C:\Users\kkaduk\Dropbox\PhD\Projects\Monkey_Ina_ECG_Respiration\AGit_ECG_Respiration_Ina\PreProcessedData\';
if strcmp(experiment, 'Inactivation')
% filename =[addtoDropbox,filesep, monkey,filesep,   'SDTvariables' , filesep, experiment, filesep, monkey, '_SDTvarPwerte_DoubleDiffStim.xlsx' ] ;  

filename =[path_SaveFig,['Table', ControlFolder], filesep, monkey, '_SDTvarPwerte_DoubleDiffStim.xlsx' ] ;  
writetable(TablePwerte,filename,'Sheet',1,  'Range' ,'A1' )

disp(['SAVED   ', filename])

   
elseif strcmp(experiment, 'Microstimulation')
%filename = [addtoDropbox,filesep, monkey,filesep,  'SDTvariables',  filesep, experiment, filesep, 'Stat' ,filesep, monkey,'_', Stimulation, '_SDTvarPwerte_DoubleDiffStim.xlsx' ] ;
writetable(TablePwerte,filename,'Sheet',1,  'Range' ,'A1' )

disp(['SAVED   ', filename])
end

%% Determine the if the dprime differs against the H0?
for indPos = 1:length(positionExp) %contra vs ipsi selection
    for indTime = 1:length(SET.TimeLine) %pre vs post
        [p,h,stat] = signrank(d_prime.easy.(positionExp{indPos}).pre);
        d_prime.easy.(positionExp{indPos}).P_Bias_pre = p;
        [p,h,stat] = signrank(d_prime.easy.(positionExp{indPos}).post);
        d_prime.easy.(positionExp{indPos}).P_Bias_post = p;
        [p,h,stat] = signrank(d_prime.difficult.(positionExp{indPos}).pre);
        d_prime.difficult.(positionExp{indPos}).P_Bias_pre = p;
        [p,h,stat] = signrank( d_prime.difficult.(positionExp{indPos}).post);
        d_prime.difficult.(positionExp{indPos}).P_Bias_post = p;
    end
end
%% Determine the response Bias by testing the criterion again the H0?
for indPos = 1:length(positionExp) %contra vs ipsi selection
    for indTime = 1:length(SET.TimeLine) %pre vs post
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



%% %%%%%%%%  SAVE THE DATASET %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T = table([nan;nan], [nan;nan], [nan;nan]) ;
Table = [];
Data.dprime = d_prime ;
Data.criterion = criterion ;
Data.Session = Session;
if strcmp(experiment, 'Inactivation');  Data.Ctr_Session = Ctr_Session;end

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
                StimulusTyp     = repmat({'DoubleDiffStim'}, length(Values),1);
                DisDifficulty   = repmat(Difficulty (indDiff), length(Values),1);
                Monkey          = repmat(monkey, length(Values),1);
                
                if strcmp(experiment, 'Inactivation')&&  strcmp((TimeLineExp {indTime}), 'post')
                    Date        = Data.Session.(Difficulty {indDiff}).(positionExp{indPos}).(TimeLineExp{indTime})';
                    Experiment    = repmat({'Ina'}, length(Date),1);
                    
                elseif strcmp(experiment, 'Inactivation')&&  strcmp((TimeLineExp {indTime}), 'pre')
                    TimeLineExperiment = repmat({'post'}, length(Values),1);
                    Date        = Data.Ctr_Session.(Difficulty {indDiff}).(positionExp {indPos}).(TimeLineExp {indTime})';
                    Experiment    = repmat({'Ctr'}, length(Date),1);
                 elseif strcmp(experiment, 'Microstimulation') 
                            Date        = Data.Session.(Difficulty {indDiff}).(positionExp {indPos}).(TimeLineExp {indTime})';
                            if strcmp(Stimulation, 'Late')
                            Experiment    = repmat({'Late'}, length(Date),1);
                            elseif strcmp(Stimulation, 'Early')
                            Experiment    = repmat({'Early'}, length(Date),1);  
                            end
                            if  strcmp((TimeLineExp {indTime}), 'pre')
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
 filename =[path_SaveFig,['Table', ControlFolder], filesep, monkey, '_SDTvar_DoubleDiffStim' ] ;  

writetable(Table, filename , 'Delimiter', ' ')
disp(['SAVED   ', filename])

%writetable(Table, [addtoDropbox, monkey, filesep,'SDTvariables', filesep, experiment, filesep, monkey, '_SDTvar_DoubleDiffStim' ], 'Delimiter', ' ')
end

%% CRITERION - WHICH RESPONSE BIAS?
Color_MoreContra = [0 1 1]; %cyan
Color_LessContra = [1 0 1]; %magenta

for indPos = 1:length(positionExp) %contra vs ipsi selection
    if ResponseBias_DisplayColor == 1
        Color_MoreContra    = [0 1 1]; %cyan
        Color_LessContra    = [1 0 1]; %magenta
        Color_Neutral       = [0.5 0.5 0.5];
    elseif ResponseBias_DisplayColor == 0
        Color_MoreContra    = [0 0 0];
        Color_LessContra    = [0 0 0];
        Color_Neutral       = [0 0 0];
    end
    for indDiff = 1: length(Difficulty)
        if  criterion.(Difficulty{indDiff}).(positionExp{indPos}).P_Bias_pre < 0.05
            if strcmp(positionExp{indPos}, 'contra')
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
        Color_MoreContra    = [0 1 1]; %cyan
        Color_LessContra    = [1 0 1]; %magenta
        Color_Neutral       = [0.5 0.5 0.5];
    elseif ResponseBias_DisplayColor == 0
        Color_MoreContra    = [0 0 1];
        Color_LessContra    = [0 0 1];
        Color_Neutral       = [0 0 1];
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
            criterion.(Difficulty{indDiff}).(positionExp{indPos}).Bias_post_Color = Color_Neutral ;
        end
    end
end





%% Criterion vs Dprime
figure('Position',[200 200 1200 900],'PaperPositionMode','auto'); % ,'PaperOrientation','landscape'
set(gcf,'Name','Dprime vs Criterion');
set(gcf,'Color',[1 1 1]);


ha(indPos) = subplot(1,2,1);
for i = 1: length(d_prime.easy.contra.pre)
    plot([d_prime.easy.contra.pre(i),d_prime.easy.contra.post(i)],[criterion.easy.contra.pre(i),criterion.easy.contra.post(i)], 'o','color',[Color.Contra, 0.4] , 'MarkerSize',MarkSize_GraphFAR_HR_PerSession,'markerfacecolor',[1 1 1 ]); hold on;
    plot([d_prime.easy.ipsi.pre(i),d_prime.easy.ipsi.post(i)], [-criterion.easy.ipsi.pre(i),-criterion.easy.ipsi.post(i)], 'o','color',[Color.Ipsi, 0.4] , 'MarkerSize',MarkSize_GraphFAR_HR_PerSession,'markerfacecolor',[1 1 1 ]); hold on;% reverse direction of criterion for ipsi
    plot(d_prime.easy.contra.post(i),criterion.easy.contra.post(i), 'o','color',Color.Contra , 'MarkerSize',MarkSize_GraphFAR_HR_PerSession,'markerfacecolor',Color.Contra); hold on;
    plot(d_prime.easy.ipsi.post(i),-criterion.easy.ipsi.post(i), 'o','color',Color.Ipsi , 'MarkerSize',MarkSize_GraphFAR_HR_PerSession,'markerfacecolor',Color.Ipsi); hold on;% reverse direction of criterion for ipsi
end
plot([nanmean(d_prime.easy.contra.pre),nanmean(d_prime.easy.contra.post) ],[nanmean(criterion.easy.contra.pre),nanmean(criterion.easy.contra.post)], 'o-','color',Color.Contra , 'MarkerSize',MarkSize_GraphFAR_HR,'markerfacecolor',[1 1 1 ], 'LineWidth',LineWidthSize); hold on;
plot([nanmean(d_prime.easy.ipsi.pre),nanmean(d_prime.easy.ipsi.post)],[-nanmean(criterion.easy.ipsi.pre) ,-nanmean(criterion.easy.ipsi.post)], 'o-','color',Color.Ipsi , 'MarkerSize',MarkSize_GraphFAR_HR,'markerfacecolor',[1 1 1 ], 'LineWidth',LineWidthSize); hold on;% reverse direction of criterion for ipsi

plot(nanmean(d_prime.easy.contra.post),nanmean(criterion.easy.contra.post), 'o','color',Color.Contra , 'MarkerSize',MarkSize_GraphFAR_HR,'markerfacecolor',Color.Contra); hold on;
plot(nanmean(d_prime.easy.ipsi.post),-nanmean(criterion.easy.ipsi.post), 'o','color',Color.Ipsi , 'MarkerSize',MarkSize_GraphFAR_HR,'markerfacecolor',Color.Ipsi); hold on;% reverse direction of criterion for ipsi

if d_prime.easy.ipsi.pvalue < 0.05
    ymax = min(-nanmean(criterion.easy.ipsi.pre) ,-nanmean([criterion.easy.ipsi.post])) -0.3;
    ext_sigline([nanmean([d_prime.easy.ipsi.pre]),nanmean([d_prime.easy.ipsi.post])],d_prime.easy.ipsi.pvalue,[ ],(ymax ),'x', Color.Ipsi); hold on;
end
if criterion.easy.ipsi.pvalue < 0.05
    y1 = -1*nanmean(criterion.easy.ipsi.pre) ;
    y2 = -1*nanmean([criterion.easy.ipsi.post]);
    ymax = max(nanmean(d_prime.easy.ipsi.pre) ,nanmean([d_prime.easy.ipsi.post])) ;
    ext_sigline([y1,y2],criterion.easy.ipsi.pvalue,[],ymax+ 0.2,'y', Color.Ipsi); hold on;
end


if d_prime.easy.contra.pvalue < 0.05
    ymax = max(nanmean(criterion.easy.contra.pre) ,nanmean([criterion.easy.contra.post]));
    ext_sigline([nanmean([d_prime.easy.contra.pre]),nanmean([d_prime.easy.contra.post])],d_prime.easy.contra.pvalue,[ ],(ymax +0.2),'x', Color.Contra); hold on;
end
if criterion.easy.contra.pvalue < 0.05
    y1 = nanmean(criterion.easy.contra.pre) ;
    y2 = nanmean([criterion.easy.contra.post]);
    ymax = max(nanmean(d_prime.easy.contra.pre) ,nanmean([d_prime.easy.contra.post])) ;
    ext_sigline([y1,y2],criterion.easy.contra.pvalue,[],ymax +0.2,'y', Color.Contra); hold on;
end
axis square
xlabel('sensitivity','fontsize',fs,'fontweight','b', 'Interpreter', 'none')
ylabel('criterion','fontsize',fs,'fontweight','b', 'Interpreter', 'none')
set(gca,'ylim',[-3 3],'xlim',[1 5],'fontsize',fs)
title(' easy distractor')
text(xlim_SDT_diff(1) +0.1,-2.8, 'more Contra (ipsi:NoGo, contra:Go)' , 'Color', 'k','fontsize',18)
text(xlim_SDT_diff(1) +0.1,2.8, 'less Contra (ipsi:Go, contra:NoGo)', 'Color', 'k','fontsize',18)

ha(indPos) = subplot(1,2,2);
for i = 1: length(d_prime.difficult.contra.pre)
    plot([d_prime.difficult.contra.pre(i),d_prime.difficult.contra.post(i)],[criterion.difficult.contra.pre(i),criterion.difficult.contra.post(i)], 'o','color',[Color.Contra, 0.4] , 'MarkerSize',MarkSize_CritDpr_small,'markerfacecolor',[1 1 1 ]); hold on;
    plot([d_prime.difficult.ipsi.pre(i),d_prime.difficult.ipsi.post(i)], [-criterion.difficult.ipsi.pre(i),-criterion.difficult.ipsi.post(i)], 'o','color',[Color.Ipsi 0.4] , 'MarkerSize',MarkSize_CritDpr_small,'markerfacecolor',[1 1 1 ]); hold on;% reverse direction of criterion for ipsi
    plot(d_prime.difficult.contra.post(i),criterion.difficult.contra.post(i), 'o','color',Color.Contra , 'MarkerSize',MarkSize_CritDpr_small,'markerfacecolor',Color.Contra); hold on;
    plot(d_prime.difficult.ipsi.post(i),-criterion.difficult.ipsi.post(i), 'o','color',Color.Ipsi , 'MarkerSize',MarkSize_CritDpr_small,'markerfacecolor',Color.Ipsi); hold on;% reverse direction of criterion for ipsi
end
plot([nanmean(d_prime.difficult.contra.pre),nanmean(d_prime.difficult.contra.post)],[nanmean(criterion.difficult.contra.pre) ,nanmean(criterion.difficult.contra.post)], 'o-','color',Color.Contra , 'MarkerSize',MarkSize_CritDpr,'markerfacecolor',[1 1 1 ], 'LineWidth', LineWidthSize); hold on;
plot([nanmean(d_prime.difficult.ipsi.pre),nanmean(d_prime.difficult.ipsi.post)],[ -nanmean(criterion.difficult.ipsi.pre) ,-nanmean(criterion.difficult.ipsi.post)], 'o-','color',Color.Ipsi, 'MarkerSize',MarkSize_CritDpr,'markerfacecolor',[1 1 1 ], 'LineWidth', LineWidthSize); hold on;% reverse direction of criterion for ipsi

plot(nanmean(d_prime.difficult.contra.post),nanmean(criterion.difficult.contra.post), 'o','color',Color.Contra , 'MarkerSize',MarkSize_CritDpr,'markerfacecolor',Color.Contra); hold on;
plot(nanmean(d_prime.difficult.ipsi.post),-nanmean(criterion.difficult.ipsi.post), 'o','color',Color.Ipsi, 'MarkerSize',MarkSize_CritDpr,'markerfacecolor',Color.Ipsi); hold on;% reverse direction of criterion for ipsi

if d_prime.difficult.ipsi.pvalue < 0.05
    ymax = max(nanmean(criterion.difficult.ipsi.pre) ,nanmean([criterion.difficult.ipsi.post]));
    ext_sigline([nanmean([d_prime.difficult.ipsi.pre]),nanmean([d_prime.difficult.ipsi.post])],d_prime.difficult.ipsi.pvalue,[ ],(ymax -1 ),'x', Color.Ipsi); hold on;
end
if criterion.difficult.ipsi.pvalue < 0.05
    y1 = -1*nanmean(criterion.difficult.ipsi.pre) ;
    y2 = -1*nanmean([criterion.difficult.ipsi.post]);
    ymax = max(nanmean(d_prime.difficult.ipsi.pre) ,nanmean([d_prime.difficult.ipsi.post])) ;
    ext_sigline([y1,y2],criterion.difficult.ipsi.pvalue,[],ymax+ 0.5,'y',Color.Ipsi); hold on;
end


if d_prime.difficult.contra.pvalue < 0.05
    ymax = max(nanmean(criterion.difficult.contra.pre) ,nanmean([criterion.difficult.contra.post]));
    ext_sigline([nanmean([d_prime.difficult.contra.pre]),nanmean([d_prime.difficult.contra.post])],d_prime.difficult.contra.pvalue,[ ],(ymax+ 0.8 ),'x', Color.Contra); hold on;
end
if criterion.difficult.contra.pvalue < 0.05
    y1 = nanmean(criterion.difficult.contra.pre) ;
    y2 = nanmean([criterion.difficult.contra.post]);
    ymax = max(nanmean(d_prime.difficult.contra.pre) ,nanmean([d_prime.difficult.contra.post])) ;
    ext_sigline([y1,y2],criterion.difficult.contra.pvalue,[],ymax+ 0.3,'y', Color.Contra); hold on;
end
axis square
xlabel('sensitivity','fontsize',fs,'fontweight','b', 'Interpreter', 'none')
ylabel('criterion','fontsize',fs,'fontweight','b', 'Interpreter', 'none')
set(gca,'ylim',[-3 3],'xlim',[0 4],'fontsize',fs)

text(0.2,-2.8, 'more Contra (ipsi:NoGo, contra:Go)', 'Color', 'k','fontsize',18)
text(0.2,2.8, 'less Contra (ipsi:Go, contra:NoGo)', 'Color', 'k','fontsize',18)

%% save the graphs
if SaveGraph
    h= figure(1);
    %savefig(h, [path_SaveFig ,filesep,'fig',filesep,  'DoubleStimuli_Crit_Dprime_' monkey,'_' experiment ,'_' folder '.fig'])
    print(h,[path_SaveFig ,filesep,  ['png', ControlFolder]  ,filesep,  'TarDistStimuli_Crit_Dprime' monkey,'_' experiment,'_' folder '.png'], '-dpng')
    set(h,'Renderer','Painters');
    set(h,'PaperPositionMode','auto')
    compl_filename =  [path_SaveFig , ['ai', ControlFolder]  ,filesep,  'TarDistStimuli_Crit_Dprime' monkey,'_' experiment ,'_', folder ,'.ai'] ;
    print(h,'-depsc',compl_filename);
    close all;
end


%% Criterion vs Dprime: direct comparison
figure('Position',[200 200 1200 900],'PaperPositionMode','auto'); % ,'PaperOrientation','landscape'
set(gcf,'Name','Dprime vs Criterion - DirectComparison');
set(gcf,'Color',[1 1 1]);




ha(indPos) = subplot(1,2,1);
for i = 1: length(d_prime.easy.contra.pre_direct)
    plot([d_prime.easy.contra.pre_direct(i),d_prime.easy.contra.post_direct(i)],[criterion.easy.contra.pre_direct(i),criterion.easy.contra.post_direct(i)], 'o','color',[Color.Contra, 0.4] , 'MarkerSize',MarkSize_CritDpr_small,'markerfacecolor',[1 1 1 ]); hold on;
    plot([d_prime.easy.ipsi.pre_direct(i),d_prime.easy.ipsi.post_direct(i)], [-criterion.easy.ipsi.pre_direct(i),-criterion.easy.ipsi.post_direct(i)], 'o','color',[Color.Ipsi, 0.4] , 'MarkerSize',MarkSize_CritDpr_small,'markerfacecolor',[1 1 1 ]); hold on;% reverse direction of criterion for ipsi
    plot(d_prime.easy.contra.post_direct(i),criterion.easy.contra.post_direct(i), 'o','color',Color.Contra , 'MarkerSize',MarkSize_CritDpr_small,'markerfacecolor',Color.Contra); hold on;
    plot(d_prime.easy.ipsi.post_direct(i),-criterion.easy.ipsi.post_direct(i), 'o','color',Color.Ipsi , 'MarkerSize',MarkSize_CritDpr_small,'markerfacecolor',Color.Ipsi); hold on;% reverse direction of criterion for ipsi
end
plot([nanmean(d_prime.easy.contra.pre_direct),nanmean(d_prime.easy.contra.post_direct) ],[nanmean(criterion.easy.contra.pre_direct),nanmean(criterion.easy.contra.post_direct)], 'o-','color',Color.Contra, 'MarkerSize',MarkSize_CritDpr,'markerfacecolor',[1 1 1 ]); hold on;
plot([nanmean(d_prime.easy.ipsi.pre_direct),nanmean(d_prime.easy.ipsi.post_direct)],[-nanmean(criterion.easy.ipsi.pre_direct) ,-nanmean(criterion.easy.ipsi.post_direct)], 'o-','color',Color.Ipsi , 'MarkerSize',MarkSize_CritDpr,'markerfacecolor',[1 1 1 ]); hold on;% reverse direction of criterion for ipsi

plot(nanmean(d_prime.easy.contra.post_direct),nanmean(criterion.easy.contra.post_direct), 'o','color',Color.Contra , 'MarkerSize',MarkSize_CritDpr,'markerfacecolor',Color.Contra); hold on;
plot(nanmean(d_prime.easy.ipsi.post_direct),-nanmean(criterion.easy.ipsi.post_direct), 'o','color',Color.Ipsi , 'MarkerSize',MarkSize_CritDpr,'markerfacecolor',Color.Ipsi); hold on;% reverse direction of criterion for ipsi

if d_prime.easy.ipsi.pvalue_direct < 0.05
    ymax = min(-nanmean(criterion.easy.ipsi.pre_direct) ,-nanmean([criterion.easy.ipsi.post_direct])) -0.3;
    ext_sigline([nanmean([d_prime.easy.ipsi.pre_direct]),nanmean([d_prime.easy.ipsi.post_direct])],d_prime.easy.ipsi.pvalue_direct,[ ],(ymax ),'x', Color.Ipsi); hold on;
end
if criterion.easy.ipsi.pvalue_direct < 0.05
    y1 = -1*nanmean(criterion.easy.ipsi.pre_direct) ;
    y2 = -1*nanmean([criterion.easy.ipsi.post_direct]);
    ymax = max(nanmean(d_prime.easy.ipsi.pre_direct) ,nanmean([d_prime.easy.ipsi.post_direct])) ;
    ext_sigline([y1,y2],criterion.easy.ipsi.pvalue_direct,[],ymax+ 0.2,'y', Color.Ipsi); hold on;
end


if d_prime.easy.contra.pvalue_direct < 0.05
    ymax = max(nanmean(criterion.easy.contra.pre_direct) ,nanmean([criterion.easy.contra.post_direct]));
    ext_sigline([nanmean([d_prime.easy.contra.pre_direct]),nanmean([d_prime.easy.contra.post_direct])],d_prime.easy.contra.pvalue_direct,[ ],(ymax +0.2),'x', Color.Contra); hold on;
end
if criterion.easy.contra.pvalue_direct < 0.05
    y1 = nanmean(criterion.easy.contra.pre_direct) ;
    y2 = nanmean([criterion.easy.contra.post_direct]);
    ymax = max(nanmean(d_prime.easy.contra.pre_direct) ,nanmean([d_prime.easy.contra.post_direct])) ;
    ext_sigline([y1,y2],criterion.easy.contra.pvalue_direct,[],ymax +0.2,'y', Color.Contra); hold on;
end
axis square
xlabel('sensitivity','fontsize',fs,'fontweight','b', 'Interpreter', 'none')
ylabel('criterion','fontsize',fs,'fontweight','b', 'Interpreter', 'none')
set(gca,'ylim',[-3 3],'xlim',[1 5],'fontsize',fs)
title(' easy distractor')
text(1.5,-2.8, 'more Contra & ipsi:NoGo, contra:Go' , 'Color', 'k','fontsize',14)
text(1.5,2.8, 'less Contra & ipsi:Go, contra:NoGo', 'Color', 'k','fontsize',14)

ha(indPos) = subplot(1,2,2);
for i = 1: length(d_prime.difficult.contra.pre_direct)
    plot([d_prime.difficult.contra.pre_direct(i),d_prime.difficult.contra.post_direct(i)],[criterion.difficult.contra.pre_direct(i),criterion.difficult.contra.post_direct(i)], 'o','color',[Color.Contra, 0.4] , 'MarkerSize',MarkSize_CritDpr_small,'markerfacecolor',[1 1 1 ]); hold on;
    plot([d_prime.difficult.ipsi.pre_direct(i),d_prime.difficult.ipsi.post_direct(i)], [-criterion.difficult.ipsi.pre_direct(i),-criterion.difficult.ipsi.post_direct(i)], 'o','color',[Color.Ipsi, 0.4] , 'MarkerSize',MarkSize_CritDpr_small,'markerfacecolor',[1 1 1 ]); hold on;% reverse direction of criterion for ipsi
    plot(d_prime.difficult.contra.post_direct(i),criterion.difficult.contra.post_direct(i), 'o','color',Color.Contra , 'MarkerSize',MarkSize_CritDpr_small,'markerfacecolor',Color.Contra); hold on;
    plot(d_prime.difficult.ipsi.post_direct(i),-criterion.difficult.ipsi.post_direct(i), 'o','color',Color.Ipsi , 'MarkerSize',MarkSize_CritDpr_small,'markerfacecolor',Color.Ipsi); hold on;% reverse direction of criterion for ipsi
end
plot([nanmean(d_prime.difficult.contra.pre_direct),nanmean(d_prime.difficult.contra.post_direct)],[nanmean(criterion.difficult.contra.pre_direct) ,nanmean(criterion.difficult.contra.post_direct)], 'o-','color',Color.Contra , 'MarkerSize',MarkSize_CritDpr,'markerfacecolor',[1 1 1 ]); hold on;
plot([nanmean(d_prime.difficult.ipsi.pre_direct),nanmean(d_prime.difficult.ipsi.post_direct)],[ -nanmean(criterion.difficult.ipsi.pre_direct) ,-nanmean(criterion.difficult.ipsi.post_direct)], 'o-','color',Color.Ipsi , 'MarkerSize',MarkSize_CritDpr,'markerfacecolor',[1 1 1 ]); hold on;% reverse direction of criterion for ipsi

plot(nanmean(d_prime.difficult.contra.post_direct),nanmean(criterion.difficult.contra.post_direct), 'o','color',Color.Contra , 'MarkerSize',MarkSize_CritDpr,'markerfacecolor',Color.Contra); hold on;
plot(nanmean(d_prime.difficult.ipsi.post_direct),-nanmean(criterion.difficult.ipsi.post_direct), 'o','color',Color.Ipsi , 'MarkerSize',MarkSize_CritDpr,'markerfacecolor',Color.Ipsi); hold on;% reverse direction of criterion for ipsi

if d_prime.difficult.ipsi.pvalue_direct < 0.05
    ymax = max(nanmean(criterion.difficult.ipsi.pre_direct) ,nanmean([criterion.difficult.ipsi.post_direct]));
    ext_sigline([nanmean([d_prime.difficult.ipsi.pre_direct]),nanmean([d_prime.difficult.ipsi.post_direct])],d_prime.difficult.ipsi.pvalue_direct,[ ],(ymax -1 ),'x', Color.Ipsi); hold on;
end
if criterion.difficult.ipsi.pvalue_direct < 0.05
    y1 = -1*nanmean(criterion.difficult.ipsi.pre_direct) ;
    y2 = -1*nanmean([criterion.difficult.ipsi.post_direct]);
    ymax = max(nanmean(d_prime.difficult.ipsi.pre_direct) ,nanmean([d_prime.difficult.ipsi.post_direct])) ;
    ext_sigline([y1,y2],criterion.difficult.ipsi.pvalue_direct,[],ymax+ 0.2,'y', Color.Ipsi); hold on;
end


if d_prime.difficult.contra.pvalue_direct < 0.05
    ymax = max(nanmean(criterion.difficult.contra.pre_direct) ,nanmean([criterion.difficult.contra.post_direct]));
    ext_sigline([nanmean([d_prime.difficult.contra.pre_direct]),nanmean([d_prime.difficult.contra.post_direct])],d_prime.difficult.contra.pvalue_direct,[ ],(ymax+ 0.8 ),'x', Color.Contra); hold on;
end
if criterion.difficult.contra.pvalue_direct < 0.05
    y1 = nanmean(criterion.difficult.contra.pre_direct) ;
    y2 = nanmean([criterion.difficult.contra.post_direct]);
    ymax = max(nanmean(d_prime.difficult.contra.pre) ,nanmean([d_prime.difficult.contra.post_direct])) ;
    ext_sigline([y1,y2],criterion.difficult.contra.pvalue_direct,[],ymax+ 0.2,'y', Color.Contra); hold on;
end
axis square
xlabel('sensitivity','fontsize',fs,'fontweight','b', 'Interpreter', 'none')
ylabel('criterion','fontsize',fs,'fontweight','b', 'Interpreter', 'none')
set(gca,'ylim',[-3 3],'xlim',[-1 5],'fontsize',fs)

text(0.2,-2.8, 'more Contra & ipsi:NoGo, contra:Go', 'Color', 'k','fontsize',14)
text(0.2,2.8, 'less Contra & ipsi:Go, contra:NoGo', 'Color', 'k','fontsize',14)

if SaveGraph
    %% save the graphs
    h(1) = figure(1);
    % savefig(h, [path_SaveFig ,'fig', filesep  'TarDistStimuli_Crit_Dprime_' monkey,'_' experiment ,'_' folder '.fig'])
    print(h,[path_SaveFig ,filesep,  ['png', ControlFolder]  ,filesep,  'TarDistStimuli_Crit_Dprime_DirectComparison' monkey,'_' experiment,'_' folder '.png'], '-dpng')
    set(h,'Renderer','Painters');
    set(h,'PaperPositionMode','auto')
    compl_filename =  [path_SaveFig , ['ai', ControlFolder]  ,filesep,  'TarDistStimuli_Crit_Dprime_DirectionComparison' monkey,'_' experiment ,'_', folder ,'.ai'] ;
    print(h,'-depsc',compl_filename);
    close all;
end



%% FIXATION
figure('Position',[200 200 1200 900],'PaperPositionMode','auto'); % ,'PaperOrientation','landscape'
set(gcf,'Name',['Fixation - Miss',ControlFolder ]);
set(gcf,'Color',[1 1 1]);

ha(1) = subplot(2,2,1); %1 % 5

plot([1;2],[nanmean(pMiss.easy.ipsi.pre); nanmean(pMiss.easy.ipsi.post)], 'x-','color',Color.Black ,'MarkerSize',10,'LineWidth',2,'markerfacecolor',[1 1 1 ]); hold on;
plot([1;2],[nanmean(pHit.easy.ipsi.pre); nanmean(pHit.easy.ipsi.post)], 'o-','color',Color.Ipsi ,'MarkerSize',10,'LineWidth',2,'markerfacecolor',[1 1 1 ]); hold on;
plot([1;2],[nanmean(pFA.easy.ipsi.pre); nanmean(pFA.easy.ipsi.post)], 's-','color',Color.Ipsi ,'MarkerSize',10,'LineWidth',2,'markerfacecolor',[1 1 1 ]); hold on;

plot(2, nanmean(pHit.easy.ipsi.post), 'o','color',Color.Ipsi ,'MarkerSize',10,'markerfacecolor',Color.Ipsi); hold on;
plot(2, nanmean(pMiss.easy.ipsi.post), 'x','color',Color.Black,'MarkerSize',10,'markerfacecolor',Color.Ipsi); hold on;
plot(2, nanmean(pFA.easy.ipsi.post), 's','color',Color.Ipsi ,'MarkerSize',10,'markerfacecolor',Color.Ipsi); hold on;

%legend('tar','fix', 'dis', 'Location', 'SouthEast')
for i = 1: length(pFA.easy.ipsi.pre)
    plot([0.9;2.1],[pHit.easy.ipsi.pre(i); pHit.easy.ipsi.post(i)], 'o','color',Color.Ipsi ,'MarkerSize',8,'markerfacecolor',[1 1 1 ]); hold on;
    plot([0.9;2.1],[pMiss.easy.ipsi.pre(i); pMiss.easy.ipsi.post(i)], 'X','color',Color.Black,'MarkerSize',8,'markerfacecolor',[1 1 1 ]); hold on;
    plot([0.8;2.2],[pFA.easy.ipsi.pre(i); pFA.easy.ipsi.post(i)], 's','color',Color.Ipsi ,'MarkerSize',8,'markerfacecolor',[1 1 1 ]); hold on;
end

if pMiss.easy.ipsi.pvalue < 0.05
    ymax = nanmean([nanmean(pMiss.easy.ipsi.pre) ,nanmean(pMiss.easy.ipsi.post)]);
    ext_sigline([1,2], pMiss.easy.ipsi.pvalue,[],ymax +0.1,'x',Color.Black )
end

if pHit.easy.ipsi.pvalue < 0.05
    ymax = nanmean([nanmean(pHit.easy.ipsi.pre) ,nanmean(pHit.easy.ipsi.post)]);
    ext_sigline([1.4 1.6], pHit.easy.ipsi.pvalue,[],ymax ,'x', Color.Ipsi )
end

if pFA.easy.ipsi.pvalue < 0.05
    ymax = nanmean([nanmean(pFA.easy.ipsi.pre) ,nanmean(pFA.easy.ipsi.post)]);
    ext_sigline([1.4 1.6], pFA.easy.ipsi.pvalue,[],ymax,'x', Color.Ipsi )
end
ylabel( 'Selection (easy)','fontsize',20,'fontweight','b', 'Interpreter', 'none' );
set(gca,'xlim',[0 3],'ylim',[0 1],'Xtick',1:2,'XTickLabel',xLegend,'fontsize',20);

ha(1) = subplot(2,2,2); %1 % 5
plot([1;2],[nanmean(pHit.easy.contra.pre); nanmean(pHit.easy.contra.post)], 'o-','color',Color.Contra ,'MarkerSize',10,'LineWidth',2,'markerfacecolor',[1 1 1 ]); hold on;
plot([1;2],[nanmean(pMiss.easy.contra.pre); nanmean(pMiss.easy.contra.post)], 'x-','color',Color.Black ,'MarkerSize',10,'LineWidth',2,'markerfacecolor',[1 1 1 ]); hold on;
plot([1;2],[nanmean(pFA.easy.contra.pre); nanmean(pFA.easy.contra.post)], 's-','color',Color.Contra ,'MarkerSize',10,'LineWidth',2,'markerfacecolor',[1 1 1 ]); hold on;
plot(2, nanmean(pHit.easy.contra.post), 'o','color',Color.Contra,'MarkerSize',10,'markerfacecolor',Color.Contra); hold on;
plot(2, nanmean(pMiss.easy.contra.post), 'X','color',Color.Black ,'MarkerSize',10,'markerfacecolor',Color.Black); hold on;
plot(2, nanmean(pFA.easy.contra.post), 's','color',Color.Contra ,'MarkerSize',10,'markerfacecolor',Color.Contra); hold on;
legend('tar','fix', 'dis', 'Location', 'NorthEast')


for i = 1: length(pFA.easy.contra.pre)
    plot([0.9;2.1],[pMiss.easy.contra.pre(i); pMiss.easy.contra.post(i)], 'X','color',Color.Black ,'MarkerSize',8,'markerfacecolor',[1 1 1]); hold on;
    plot([0.9;2.1],[pHit.easy.contra.pre(i); pHit.easy.contra.post(i)], 'o','color',Color.Contra ,'MarkerSize',8,'markerfacecolor',[1 1 1]); hold on;
    plot([0.8;2.2],[pFA.easy.contra.pre(i); pFA.easy.contra.post(i)], 's','color',Color.Contra ,'MarkerSize',8,'markerfacecolor',[1 1 1]); hold on;
end

if pMiss.easy.contra.pvalue < 0.05
    ymax = nanmean([nanmean(pMiss.easy.contra.pre) ,nanmean(pMiss.easy.contra.post)]);
    ext_sigline([1,2], pMiss.easy.contra.pvalue,[],ymax +0.1,'x', Color.Black )
end
if pHit.easy.contra.pvalue < 0.05
    ymax = nanmean([nanmean(pHit.easy.contra.pre) ,nanmean(pHit.easy.contra.post)]);
    ext_sigline([1.4 1.6], pHit.easy.contra.pvalue,[],ymax,'x', Color.Contra )
end
if pFA.easy.contra.pvalue < 0.05
    ymax = nanmean([nanmean(pFA.easy.contra.pre) ,nanmean(pFA.easy.contra.post)]);
    ext_sigline([1.4 1.6], pFA.easy.contra.pvalue,[],ymax,'x', Color.Contra )
end

ylabel( 'Selection (easy)','fontsize',20,'fontweight','b', 'Interpreter', 'none' );
set(gca,'xlim',[0 3],'ylim',[0 1],'Xtick',1:2,'XTickLabel',xLegend,'fontsize',20);



ha(1) = subplot(2,2,3); %1 % 5

plot([1;2],[nanmean(pMiss.difficult.ipsi.pre); nanmean(pMiss.difficult.ipsi.post)], 'x-','color',Color.Black ,'MarkerSize',10,'LineWidth',2,'markerfacecolor',[1 1 1 ]); hold on;
plot([1;2],[nanmean(pHit.difficult.ipsi.pre); nanmean(pHit.difficult.ipsi.post)], 'o-','color',Color.Ipsi ,'MarkerSize',10,'LineWidth',2,'markerfacecolor',[1 1 1 ]); hold on;
plot([1;2],[nanmean(pFA.difficult.ipsi.pre); nanmean(pFA.difficult.ipsi.post)], 's-','color',Color.Ipsi ,'MarkerSize',10,'LineWidth',2,'markerfacecolor',[1 1 1 ]); hold on;

plot(2, nanmean(pHit.difficult.ipsi.post), 'o','color',Color.Ipsi ,'MarkerSize',10,'markerfacecolor',Color.Ipsi); hold on;
plot(2, nanmean(pMiss.difficult.ipsi.post), 'X','color',Color.Black ,'MarkerSize',10,'markerfacecolor',Color.Black); hold on;
plot(2, nanmean(pFA.difficult.ipsi.post), 's','color',Color.Ipsi ,'MarkerSize',10,'markerfacecolor',Color.Ipsi); hold on;

%legend('tar','fix', 'dis', 'Location', 'SouthEast')
for i = 1: length(pFA.difficult.ipsi.pre)
    plot([0.9;2.1],[pHit.difficult.ipsi.pre(i); pHit.difficult.ipsi.post(i)], 'o','color',Color.Ipsi ,'MarkerSize',8,'markerfacecolor',[1 1 1 ]); hold on;
    plot([0.9;2.1],[pMiss.difficult.ipsi.pre(i); pMiss.difficult.ipsi.post(i)], 'X','color',Color.Black,'MarkerSize',8,'markerfacecolor',[1 1 1 ]); hold on;
    plot([0.8;2.2],[pFA.difficult.ipsi.pre(i); pFA.difficult.ipsi.post(i)], 's','color',Color.Ipsi ,'MarkerSize',8,'markerfacecolor',[1 1 1 ]); hold on;
end

if pMiss.difficult.ipsi.pvalue < 0.05
    ymax = nanmean([nanmean(pMiss.difficult.ipsi.pre) ,nanmean(pMiss.difficult.ipsi.post)]);
    ext_sigline([1,2], pMiss.difficult.ipsi.pvalue,[],ymax +0.1,'x', Color.Black )
end

if pHit.difficult.ipsi.pvalue < 0.05
    ymax = nanmean([nanmean(pHit.difficult.ipsi.pre) ,nanmean(pHit.difficult.ipsi.post)]);
    ext_sigline([1.4 1.6], pHit.difficult.ipsi.pvalue,[],ymax ,'x', Color.Ipsi )
end

if pFA.difficult.ipsi.pvalue < 0.05
    ymax = nanmean([nanmean(pFA.difficult.ipsi.pre) ,nanmean(pFA.difficult.ipsi.post)]);
    ext_sigline([1.4 1.6], pFA.difficult.ipsi.pvalue,[],ymax,'x', Color.Ipsi )
end
ylabel( 'Selection (difficult)','fontsize',20,'fontweight','b', 'Interpreter', 'none' );
set(gca,'xlim',[0 3],'ylim',[0 1],'Xtick',1:2,'XTickLabel',xLegend,'fontsize',20);

ha(1) = subplot(2,2,4); %1 % 5

plot([1;2],[nanmean(pHit.difficult.contra.pre); nanmean(pHit.difficult.contra.post)], 'o-','color',Color.Contra ,'MarkerSize',10,'LineWidth',2,'markerfacecolor',[1 1 1 ]); hold on;
plot([1;2],[nanmean(pMiss.difficult.contra.pre); nanmean(pMiss.difficult.contra.post)], 'x-','color',Color.Black ,'MarkerSize',10,'LineWidth',2,'markerfacecolor',[1 1 1 ]); hold on;
plot([1;2],[nanmean(pFA.difficult.contra.pre); nanmean(pFA.difficult.contra.post)], 's-','color',Color.Contra ,'MarkerSize',10,'LineWidth',2,'markerfacecolor',[1 1 1 ]); hold on;
plot(2, nanmean(pHit.difficult.contra.post), 'o','color',Color.Contra,'MarkerSize',10,'markerfacecolor',Color.Contra); hold on;
plot(2, nanmean(pMiss.difficult.contra.post), 'X','color',Color.Black ,'MarkerSize',10,'markerfacecolor',Color.Black); hold on;
plot(2, nanmean(pFA.difficult.contra.post), 's','color',Color.Contra ,'MarkerSize',10,'markerfacecolor',Color.Contra); hold on;
legend('tar','fix', 'dis', 'Location', 'NorthEast')


for i = 1: length(pFA.difficult.contra.pre)
    plot([0.9;2.1],[pHit.difficult.contra.pre(i); pHit.difficult.contra.post(i)], 'o','color',Color.Contra ,'MarkerSize',8,'markerfacecolor',[1 1 1]); hold on;
    plot([0.9;2.1],[pMiss.difficult.contra.pre(i); pMiss.difficult.contra.post(i)], 'X','color',Color.Black ,'MarkerSize',8,'markerfacecolor',[1 1 1]); hold on;
    plot([0.8;2.2],[pFA.difficult.contra.pre(i); pFA.difficult.contra.post(i)], 's','color',Color.Contra ,'MarkerSize',8,'markerfacecolor',[1 1 1]); hold on;
end

if pMiss.difficult.contra.pvalue < 0.05
    ymax = nanmean([nanmean(pMiss.difficult.contra.pre) ,nanmean(pMiss.difficult.contra.post)]);
    ext_sigline([1,2], pMiss.difficult.contra.pvalue,[],ymax +0.1,'x', Color.Black )
end
if pHit.difficult.contra.pvalue < 0.05
    ymax = nanmean([nanmean(pHit.difficult.contra.pre) ,nanmean(pHit.difficult.contra.post)]);
    ext_sigline([1.4 1.6], pHit.difficult.contra.pvalue,[],ymax,'x', Color.Contra )
end
if pFA.difficult.contra.pvalue < 0.05
    ymax = nanmean([nanmean(pFA.difficult.contra.pre) ,nanmean(pFA.difficult.contra.post)]);
    ext_sigline([1.4 1.6], pFA.difficult.contra.pvalue,[],ymax,'x', Color.Contra)
end

ylabel( 'Selection (difficult)','fontsize',20,'fontweight','b', 'Interpreter', 'none' );
set(gca,'xlim',[0 3],'ylim',[0 1],'Xtick',1:2,'XTickLabel',xLegend,'fontsize',20);


if SaveGraph
    h = figure(1);
    % savefig(h, [path_SaveFig ,'fig', filesep,  'TarDistStimuli_HR_FAR_' monkey,'_' experiment ,'_' folder ,'.fig'])
    print(h,[path_SaveFig ,filesep,  ['png', ControlFolder] ,filesep,  'TarDistStimuli_SELECTION' monkey,'.png'], '-dpng') %'_' experiment ,'_' folder,
    set(h,'Renderer','Painters');
    set(h,'PaperPositionMode','auto')
    compl_filename =  [path_SaveFig , ['ai', ControlFolder]  ,filesep,  'TarDistStimuli_SELECTION' monkey,'_' experiment ,'_', folder ,'.ai'] ;
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
    plot(pFA.easy.ipsi.post(i), pHit.easy.ipsi.post(i), 'o' ,'color',Color.Ipsi, 'MarkerSize',MarkSize_GraphFAR_HR_PerSession, 'MarkerFaceColor',Color.Ipsi,'LineWidth', 2); hold on;
    plot(pFA.easy.contra.post(i), pHit.easy.contra.post(i), 'o','color', Color.Contra , 'MarkerSize',MarkSize_GraphFAR_HR_PerSession,'markerfacecolor', Color.Contra,'LineWidth', 2); hold on;
    
    line([pFA.easy.ipsi.pre(i),pFA.easy.ipsi.post(i)], [pHit.easy.ipsi.pre(i),pHit.easy.ipsi.post(i)],'Color',[Color.Ipsi, 0.3] , 'MarkerSize',MarkSize_GraphFAR_HR_PerSession,'markerfacecolor',[1 1 1],'LineWidth', 2); hold on;
    line([pFA.easy.contra.pre(i),pFA.easy.contra.post(i)], [pHit.easy.contra.pre(i),pHit.easy.contra.post(i)],'Color',[Color.Contra, 0.3] , 'MarkerSize',MarkSize_GraphFAR_HR_PerSession,'markerfacecolor',[1 1 1],'LineWidth', 2); hold on;
    
end
%ipsi Post
plot([nanmean(pFA.easy.ipsi.pre),nanmean(pFA.easy.ipsi.post)], [nanmean(pHit.easy.ipsi.pre),nanmean(pHit.easy.ipsi.post)], 'o-','color',Color.Ipsi , 'MarkerSize',MarkSize_GraphFAR_HR,'markerfacecolor',[1 1 1],'LineWidth', 3); hold on;
plot([nanmean(pFA.easy.contra.pre),nanmean(pFA.easy.contra.post)], [nanmean(pHit.easy.contra.pre),nanmean(pHit.easy.contra.post)], 'o-','color', Color.Contra , 'MarkerSize',MarkSize_GraphFAR_HR,'markerfacecolor',[1 1 1],'LineWidth', 3); hold on;
%Fill the circle which is post
plot(nanmean(pFA.easy.ipsi.post), nanmean(pHit.easy.ipsi.post), 'o','color',Color.Ipsi , 'MarkerSize',MarkSize_GraphFAR_HR-1,'markerfacecolor',Color.Ipsi,'LineWidth', 2); hold on;
plot(nanmean(pFA.easy.contra.post), nanmean(pHit.easy.contra.post), 'o','color', Color.Contra , 'MarkerSize',MarkSize_GraphFAR_HR-1,'markerfacecolor', Color.Contra,'LineWidth', 2); hold on;
legend('ipsi', 'contra')

title('easy')

line([0 1],[1 0],'Color',[0 0 0],'LineStyle',':');
set(gca,'ylim',[0 1],'xlim',[0 1],'fontsize',fs)
xlabel( 'FA rate','fontsize',fs,'fontweight','b', 'Interpreter', 'none' );
ylabel( 'Hitrate','fontsize',fs,'fontweight','b', 'Interpreter', 'none' );
axis square


if pFA.easy.ipsi.pvalue < 0.05
    y1 = nanmean(pFA.easy.ipsi.pre) ;
    y2 = nanmean(pFA.easy.ipsi.post);
    ymax = min(nanmean(pHit.easy.ipsi.pre) ,nanmean(pHit.easy.ipsi.post));
    ext_sigline([y1,y2],pFA.easy.ipsi.pvalue,[],ymax -0.2,'x', Color.Ipsi); hold on;
end
if pHit.easy.ipsi.pvalue < 0.05
    y1 = nanmean(pHit.easy.ipsi.pre) ;
    y2 = nanmean(pHit.easy.ipsi.post);
    ymax = max(nanmean(pFA.easy.ipsi.pre) ,nanmean(pFA.easy.ipsi.post));
    ext_sigline([y1,y2],pHit.easy.ipsi.pvalue,[],ymax +0.2,'y',Color.Ipsi); hold on;
end
if pFA.easy.contra.pvalue < 0.05
    y1 = nanmean(pFA.easy.contra.pre) ;
    y2 = nanmean(pFA.easy.contra.post);
    ymax = min(nanmean(pHit.easy.contra.pre) ,nanmean(pHit.easy.contra.post));
    ext_sigline([y1,y2],pFA.easy.contra.pvalue,[],ymax -0.2,'x', Color.Contra); hold on;
end
if pHit.easy.contra.pvalue < 0.05
    y1 = nanmean(pHit.easy.contra.pre) ;
    y2 = nanmean(pHit.easy.contra.post);
    ymax = max(nanmean(pFA.easy.contra.pre) ,nanmean(pFA.easy.contra.post));
    ext_sigline([y1,y2],pHit.easy.contra.pvalue,[],ymax +0.2,'y', Color.Contra); hold on;
end



%%

ha(2) = subplot(1,2,2);
%ipsi Post

for i = 1: length(pFA.easy.ipsi.pre)
    plot(pFA.difficult.ipsi.post(i), pHit.difficult.ipsi.post(i), 'o' ,'color',Color.Ipsi, 'MarkerSize',MarkSize_CritDpr_small, 'MarkerFaceColor',Color.Ipsi,'LineWidth', 2); hold on;
    plot(pFA.difficult.contra.post(i), pHit.difficult.contra.post(i), 'o','color', Color.Contra , 'MarkerSize',MarkSize_CritDpr_small,'markerfacecolor', Color.Contra,'LineWidth', 2); hold on;
    
    line([pFA.difficult.ipsi.pre(i),pFA.difficult.ipsi.post(i)], [pHit.difficult.ipsi.pre(i),pHit.difficult.ipsi.post(i)],'Color',[Color.Ipsi, 0.3] , 'MarkerSize',MarkSize_CritDpr_small,'markerfacecolor',[1 1 1],'LineWidth', 2); hold on;
    line([pFA.difficult.contra.pre(i),pFA.difficult.contra.post(i)], [pHit.difficult.contra.pre(i),pHit.difficult.contra.post(i)],'Color',[ Color.Contra, 0.3] , 'MarkerSize',MarkSize_CritDpr_small,'markerfacecolor',[1 1 1],'LineWidth', 2); hold on;
    
    %Fill the circle which is post
end
plot([nanmean(pFA.difficult.ipsi.pre),nanmean(pFA.difficult.ipsi.post)], [nanmean(pHit.difficult.ipsi.pre),nanmean(pHit.difficult.ipsi.post)], 'o-','color',Color.Ipsi , 'MarkerSize',MarkSize_CritDpr,'markerfacecolor',[1 1 1],'LineWidth', 2); hold on;
plot([nanmean(pFA.difficult.contra.pre),nanmean(pFA.difficult.contra.post)], [nanmean(pHit.difficult.contra.pre),nanmean(pHit.difficult.contra.post)], 'o-','color', Color.Contra , 'MarkerSize',MarkSize_CritDpr,'markerfacecolor',[1 1 1],'LineWidth', 2); hold on;
plot(nanmean(pFA.difficult.ipsi.post), nanmean(pHit.difficult.ipsi.post), 'o','color',Color.Ipsi , 'MarkerSize',MarkSize_CritDpr-1,'markerfacecolor',Color.Ipsi,'LineWidth', 2); hold on;
plot(nanmean(pFA.difficult.contra.post), nanmean(pHit.difficult.contra.post), 'o','color',Color.Contra , 'MarkerSize',MarkSize_CritDpr-1,'markerfacecolor', Color.Contra,'LineWidth', 2); hold on;
legend('ipsi', 'contra', 'Location', 'South')

%title('difficult')

line([0 1],[1 0],'Color',[0 0 0],'LineStyle',':');
set(gca,'ylim',[0 1],'xlim',[0 1],'fontsize',fs)
xlabel( 'FA rate','fontsize',fs,'fontweight','b', 'Interpreter', 'none' );
ylabel( 'Hitrate','fontsize',fs,'fontweight','b', 'Interpreter', 'none' );
axis square


if pFA.difficult.ipsi.pvalue < 0.05
    y1 = nanmean(pFA.difficult.ipsi.pre) ;
    y2 = nanmean(pFA.difficult.ipsi.post);
    ymax = min(nanmean(pHit.difficult.ipsi.pre) ,nanmean(pHit.difficult.ipsi.post));
    ext_sigline([y1,y2],pFA.difficult.ipsi.pvalue,[],ymax -0.3,'x',  Color.Ipsi); hold on;
end
if pHit.difficult.ipsi.pvalue < 0.05
    y1 = nanmean(pHit.difficult.ipsi.pre) ;
    y2 = nanmean(pHit.difficult.ipsi.post);
    ymax = max(nanmean(pFA.difficult.ipsi.pre) ,nanmean(pFA.difficult.ipsi.post));
    ext_sigline([y1,y2],pHit.difficult.ipsi.pvalue,[],ymax +0.2,'y',Color.Ipsi); hold on;
end
if pFA.difficult.contra.pvalue < 0.05
    y1 = nanmean(pFA.difficult.contra.pre) ;
    y2 = nanmean(pFA.difficult.contra.post);
    ymax = min(nanmean(pHit.difficult.contra.pre) ,nanmean(pHit.difficult.contra.post));
    ext_sigline([y1,y2],pFA.difficult.contra.pvalue,[],ymax -0.2,'x', Color.Contra); hold on;
end
if pHit.difficult.contra.pvalue < 0.05
    y1 = nanmean(pHit.difficult.contra.pre) ;
    y2 = nanmean(pHit.difficult.contra.post);
    ymax = max(nanmean(pFA.difficult.contra.pre) ,nanmean(pFA.difficult.contra.post));
    ext_sigline([y1,y2],pHit.difficult.contra.pvalue,[],ymax +0.2,'y', Color.Contra); hold on;
end
if SaveGraph
    h = figure(1);
    % savefig(h, [path_SaveFig ,'fig', filesep,  'TarDistStimuli_HR_FAR_' monkey,'_' experiment ,'_' folder ,'.fig'])
    print(h,[path_SaveFig ,filesep,  ['png', ControlFolder] ,filesep,  'TarDistStimuli_HR_FAR_AfterLuoetal' monkey,'_' experiment ,'_' folder, '.png'], '-dpng')
    set(h,'Renderer','Painters');
    set(h,'PaperPositionMode','auto')
    compl_filename =  [path_SaveFig , ['ai', ControlFolder]  ,filesep,  'TarDistStimuli_HR_FAR_AfterLuoetal.' monkey,'_' experiment ,'_', folder ,'.ai'] ;
    print(h,'-depsc',compl_filename);
    close all;
end



%% FOr the TALK
%% graph - From the paper Luo & Maunsell (2018)
figure('Position',[200 200 1200 900],'PaperPositionMode','auto'); % ,'PaperOrientation','landscape'
set(gcf,'Name',['Difficult_Graph1HR_FAR_Graph2Sensitivity_Criterion',ControlFolder ]);
set(gcf,'Color',[1 1 1]);


ha(2) = subplot(1,2,1);
%ipsi Post

for i = 1: length(pFA.easy.ipsi.pre)
    plot(pFA.difficult.ipsi.post(i), pHit.difficult.ipsi.post(i), 'o' ,'color',Color.Ipsi, 'MarkerSize',MarkSize_GraphFAR_HR_PerSession, 'MarkerFaceColor',Color.Ipsi,'LineWidth', 2); hold on;
    plot(pFA.difficult.contra.post(i), pHit.difficult.contra.post(i), 'o','color',Color.Contra , 'MarkerSize',MarkSize_GraphFAR_HR_PerSession,'markerfacecolor',Color.Contra,'LineWidth', 2); hold on;
    
    line([pFA.difficult.ipsi.pre(i),pFA.difficult.ipsi.post(i)], [pHit.difficult.ipsi.pre(i),pHit.difficult.ipsi.post(i)],'Color',[Color.Ipsi, 0.3] , 'MarkerSize',MarkSize_GraphFAR_HR_PerSession,'markerfacecolor',[1 1 1],'LineWidth', 2); hold on;
    line([pFA.difficult.contra.pre(i),pFA.difficult.contra.post(i)], [pHit.difficult.contra.pre(i),pHit.difficult.contra.post(i)],'Color',[Color.Contra, 0.3] , 'MarkerSize',MarkSize_GraphFAR_HR_PerSession,'markerfacecolor',[1 1 1],'LineWidth', 2); hold on;
    
    %Fill the circle which is post
end
plot([nanmean(pFA.difficult.ipsi.pre),nanmean(pFA.difficult.ipsi.post)], [nanmean(pHit.difficult.ipsi.pre),nanmean(pHit.difficult.ipsi.post)], 'o-','color',Color.Ipsi , 'MarkerSize',MarkSize_GraphFAR_HR,'markerfacecolor',[1 1 1],'LineWidth', LineWidthSize); hold on;
plot([nanmean(pFA.difficult.contra.pre),nanmean(pFA.difficult.contra.post)], [nanmean(pHit.difficult.contra.pre),nanmean(pHit.difficult.contra.post)], 'o-','color',Color.Contra , 'MarkerSize',MarkSize_GraphFAR_HR,'markerfacecolor',[1 1 1],'LineWidth', LineWidthSize); hold on;
plot(nanmean(pFA.difficult.ipsi.post), nanmean(pHit.difficult.ipsi.post), 'o','color',Color.Ipsi , 'MarkerSize',MarkSize_GraphFAR_HR-1,'markerfacecolor',Color.Ipsi,'LineWidth', LineWidthSize); hold on;
plot(nanmean(pFA.difficult.contra.post), nanmean(pHit.difficult.contra.post), 'o','color',Color.Contra, 'MarkerSize',MarkSize_GraphFAR_HR-1,'markerfacecolor',Color.Contra,'LineWidth', LineWidthSize); hold on;
%legend('ipsi', 'contra', 'Location', 'South')

%S('difficult')

line([0 1],[1 0],'Color',[0 0 0],'LineStyle',':');
set(gca,'ylim',[0 1],'xlim',[0 1],'fontsize',fs)
xlabel( 'FA rate','fontsize',fs,'fontweight','b', 'Interpreter', 'none' );
ylabel( 'Hitrate','fontsize',fs,'fontweight','b', 'Interpreter', 'none' );
axis square


if pFA.difficult.ipsi.pvalue < 0.05
    y1 = nanmean(pFA.difficult.ipsi.pre) ;
    y2 = nanmean(pFA.difficult.ipsi.post);
    ymax = min(nanmean(pHit.difficult.ipsi.pre) ,nanmean(pHit.difficult.ipsi.post));
    ext_sigline([y1,y2],pFA.difficult.ipsi.pvalue,[],ymax -0.3,'x',Color.Ipsi); hold on;
end
if pHit.difficult.ipsi.pvalue < 0.05
    y1 = nanmean(pHit.difficult.ipsi.pre) ;
    y2 = nanmean(pHit.difficult.ipsi.post);
    ymax = max(nanmean(pFA.difficult.ipsi.pre) ,nanmean(pFA.difficult.ipsi.post));
    ext_sigline([y1,y2],pHit.difficult.ipsi.pvalue,[],ymax +0.2,'y',Color.Ipsi); hold on;
end
if pFA.difficult.contra.pvalue < 0.05
    y1 = nanmean(pFA.difficult.contra.pre) ;
    y2 = nanmean(pFA.difficult.contra.post);
    ymax = min(nanmean(pHit.difficult.contra.pre) ,nanmean(pHit.difficult.contra.post));
    ext_sigline([y1,y2],pFA.difficult.contra.pvalue,[],ymax -0.2,'x', Color.Contra); hold on;
end
if pHit.difficult.contra.pvalue < 0.05
    y1 = nanmean(pHit.difficult.contra.pre) ;
    y2 = nanmean(pHit.difficult.contra.post);
    ymax = max(nanmean(pFA.difficult.contra.pre) ,nanmean(pFA.difficult.contra.post));
    ext_sigline([y1,y2],pHit.difficult.contra.pvalue,[],ymax +0.2,'y',Color.Contra); hold on;
end


ha(indPos) = subplot(1,2,2);

for i = 1: length(d_prime.difficult.contra.pre)
    plot([d_prime.difficult.contra.pre(i),d_prime.difficult.contra.post(i)],[criterion.difficult.contra.pre(i),criterion.difficult.contra.post(i)], 'o','color',[Color.Contra 0.4] , 'MarkerSize',MarkSize_CritDpr_small,'markerfacecolor',[1 1 1 ]); hold on;
    plot([d_prime.difficult.ipsi.pre(i),d_prime.difficult.ipsi.post(i)], [-criterion.difficult.ipsi.pre(i),-criterion.difficult.ipsi.post(i)], 'o','color',[Color.Ipsi, 0.4] , 'MarkerSize',MarkSize_CritDpr_small,'markerfacecolor',[1 1 1 ]); hold on;% reverse direction of criterion for ipsi
    plot(d_prime.difficult.contra.post(i),criterion.difficult.contra.post(i), 'o','color',Color.Contra, 'MarkerSize',MarkSize_CritDpr_small,'markerfacecolor',Color.Contra); hold on;
    plot(d_prime.difficult.ipsi.post(i),-criterion.difficult.ipsi.post(i), 'o','color',Color.Ipsi , 'MarkerSize',MarkSize_CritDpr_small,'markerfacecolor',Color.Ipsi); hold on;% reverse direction of criterion for ipsi
end
plot([nanmean(d_prime.difficult.contra.pre),nanmean(d_prime.difficult.contra.post)],[nanmean(criterion.difficult.contra.pre) ,nanmean(criterion.difficult.contra.post)], 'o-','color',Color.Contra , 'MarkerSize',MarkSize_CritDpr,'markerfacecolor',[1 1 1 ],'LineWidth',LineWidthSize ); hold on;
plot([nanmean(d_prime.difficult.ipsi.pre),nanmean(d_prime.difficult.ipsi.post)],[ -nanmean(criterion.difficult.ipsi.pre) ,-nanmean(criterion.difficult.ipsi.post)], 'o-','color',Color.Ipsi , 'MarkerSize',MarkSize_CritDpr,'markerfacecolor',[1 1 1 ],'LineWidth',LineWidthSize ); hold on;% reverse direction of criterion for ipsi

plot(nanmean(d_prime.difficult.contra.post),nanmean(criterion.difficult.contra.post), 'o','color',Color.Contra, 'MarkerSize',MarkSize_CritDpr,'markerfacecolor',Color.Contra); hold on;
plot(nanmean(d_prime.difficult.ipsi.post),-nanmean(criterion.difficult.ipsi.post), 'o','color',Color.Ipsi , 'MarkerSize',MarkSize_CritDpr,'markerfacecolor',Color.Ipsi); hold on;% reverse direction of criterion for ipsi

if d_prime.difficult.ipsi.pvalue < 0.05
    ymax = max(-1*nanmean(criterion.difficult.ipsi.pre) ,nanmean(-1*[criterion.difficult.ipsi.post]));
    ext_sigline([nanmean([d_prime.difficult.ipsi.pre]),nanmean([d_prime.difficult.ipsi.post])],d_prime.difficult.ipsi.pvalue,[ ],(ymax -1 ),'x', Color.Ipsi); hold on;
end
if criterion.difficult.ipsi.pvalue < 0.05
    y1 = -1*nanmean(criterion.difficult.ipsi.pre) ;
    y2 = -1*nanmean([criterion.difficult.ipsi.post]);
    ymax = max(nanmean(d_prime.difficult.ipsi.pre) ,nanmean([d_prime.difficult.ipsi.post])) ;
    ext_sigline([y1,y2],criterion.difficult.ipsi.pvalue,[],ymax+ 0.5,'y', Color.Ipsi); hold on;
end


if d_prime.difficult.contra.pvalue < 0.05
    ymax = max(nanmean(criterion.difficult.contra.pre) ,nanmean([criterion.difficult.contra.post]));
    ext_sigline([nanmean([d_prime.difficult.contra.pre]),nanmean([d_prime.difficult.contra.post])],d_prime.difficult.contra.pvalue,[ ],(ymax+ 0.8 ),'x', Color.Contra); hold on;
end
if criterion.difficult.contra.pvalue < 0.05
    y1 = nanmean(criterion.difficult.contra.pre) ;
    y2 = nanmean([criterion.difficult.contra.post]);
    ymax = max(nanmean(d_prime.difficult.contra.pre) ,nanmean([d_prime.difficult.contra.post])) ;
    ext_sigline([y1,y2],criterion.difficult.contra.pvalue,[],ymax+ 0.3,'y', Color.Contra); hold on;
end
axis square
xlabel('sensitivity','fontsize',fs,'fontweight','b', 'Interpreter', 'none')
ylabel('criterion','fontsize',fs,'fontweight','b', 'Interpreter', 'none')
set(gca,'ylim',[-3 3],'xlim',xlim_SDT_diff,'fontsize',fs)

text(xlim_SDT_diff(1)+0.1,-2.8, 'more Contra (ipsi:NoGo, contra:Go)', 'Color', 'k','fontsize',18)
text(xlim_SDT_diff(1)+0.1,2.8, 'less Contra (ipsi:Go, contra:NoGo)', 'Color', 'k','fontsize',18)

%% save the graphs
if SaveGraph
    h= figure(1);
    %savefig(h, [path_SaveFig ,filesep,'fig',filesep,  'DoubleStimuli_Crit_Dprime_' monkey,'_' experiment ,'_' folder '.fig'])
    print(h,[path_SaveFig ,filesep,  ['png', ControlFolder]  ,filesep,  'TarDistStimuli_difficult_HR_FAR_Dprime_Criterion' '.png'], '-dpng') %monkey,'_' experiment,'_' folder
    set(h,'Renderer','Painters');
    set(h,'PaperPositionMode','auto')
    compl_filename =  [path_SaveFig , ['ai', ControlFolder]  ,filesep,  'TarDistStimuli_difficult_HR_FAR_Dprime_Criterion' monkey,'_'  ,'.ai'] ; %experiment ,'_', folder
    print(h,'-depsc',compl_filename);
    close all;
end


%% EASY - TALK GRAPHS
figure('Position',[200 200 1200 900],'PaperPositionMode','auto'); % ,'PaperOrientation','landscape'
set(gcf,'Name',['Easy_Graph1HR_FAR_Graph2Sensitivity_Criterion',ControlFolder ]);
set(gcf,'Color',[1 1 1]);



ha(2) = subplot(1,2,1);
%ipsi Post

for i = 1: length(pFA.easy.ipsi.pre)
    plot(pFA.easy.ipsi.post(i), pHit.easy.ipsi.post(i), 'o' ,'color',Color.Ipsi, 'MarkerSize',MarkSize_GraphFAR_HR_PerSession, 'MarkerFaceColor',Color.Ipsi,'LineWidth', 2); hold on;
    plot(pFA.easy.contra.post(i), pHit.easy.contra.post(i), 'o','color',Color.Contra , 'MarkerSize',MarkSize_GraphFAR_HR_PerSession,'markerfacecolor',Color.Contra,'LineWidth', 2); hold on;
    
    line([pFA.easy.ipsi.pre(i),pFA.easy.ipsi.post(i)], [pHit.easy.ipsi.pre(i),pHit.easy.ipsi.post(i)],'Color',[Color.Ipsi, 0.3] , 'MarkerSize',MarkSize_GraphFAR_HR_PerSession,'markerfacecolor',[1 1 1],'LineWidth', 2); hold on;
    line([pFA.easy.contra.pre(i),pFA.easy.contra.post(i)], [pHit.easy.contra.pre(i),pHit.easy.contra.post(i)],'Color',[Color.Contra, 0.3] , 'MarkerSize',MarkSize_GraphFAR_HR_PerSession,'markerfacecolor',[1 1 1],'LineWidth', 2); hold on;
    
    %Fill the circle which is post
end
plot([nanmean(pFA.easy.ipsi.pre),nanmean(pFA.easy.ipsi.post)], [nanmean(pHit.easy.ipsi.pre),nanmean(pHit.easy.ipsi.post)], 'o-','color',Color.Ipsi , 'MarkerSize',MarkSize_GraphFAR_HR,'markerfacecolor',[1 1 1],'LineWidth', LineWidthSize); hold on;
plot([nanmean(pFA.easy.contra.pre),nanmean(pFA.easy.contra.post)], [nanmean(pHit.easy.contra.pre),nanmean(pHit.easy.contra.post)], 'o-','color',Color.Contra , 'MarkerSize',MarkSize_GraphFAR_HR,'markerfacecolor',[1 1 1],'LineWidth', LineWidthSize); hold on;
plot(nanmean(pFA.easy.ipsi.post), nanmean(pHit.easy.ipsi.post), 'o','color',Color.Ipsi , 'MarkerSize',MarkSize_GraphFAR_HR-1,'markerfacecolor',Color.Ipsi,'LineWidth', LineWidthSize); hold on;
plot(nanmean(pFA.easy.contra.post), nanmean(pHit.easy.contra.post), 'o','color',Color.Contra , 'MarkerSize',MarkSize_GraphFAR_HR-1,'markerfacecolor',Color.Contra,'LineWidth', LineWidthSize); hold on;
%legend('ipsi', 'contra', 'Location', 'South')

%title('easy')

line([0 1],[1 0],'Color',[0 0 0],'LineStyle',':');
set(gca,'ylim',[0 1],'xlim',[0 1],'fontsize',fs)
xlabel( 'FA rate','fontsize',fs,'fontweight','b', 'Interpreter', 'none' );
ylabel( 'Hitrate','fontsize',fs,'fontweight','b', 'Interpreter', 'none' );
axis square


if pFA.easy.ipsi.pvalue < 0.05
    y1 = nanmean(pFA.easy.ipsi.pre) ;
    y2 = nanmean(pFA.easy.ipsi.post);
    ymax = min(nanmean(pHit.easy.ipsi.pre) ,nanmean(pHit.easy.ipsi.post));
    ext_sigline([y1,y2],pFA.easy.ipsi.pvalue,[],ymax -0.3,'x', Color.Ipsi); hold on;
end
if pHit.easy.ipsi.pvalue < 0.05
    y1 = nanmean(pHit.easy.ipsi.pre) ;
    y2 = nanmean(pHit.easy.ipsi.post);
    ymax = max(nanmean(pFA.easy.ipsi.pre) ,nanmean(pFA.easy.ipsi.post));
    ext_sigline([y1,y2],pHit.easy.ipsi.pvalue,[],ymax +0.2,'y',Color.Ipsi); hold on;
end
if pFA.easy.contra.pvalue < 0.05
    y1 = nanmean(pFA.easy.contra.pre) ;
    y2 = nanmean(pFA.easy.contra.post);
    ymax = min(nanmean(pHit.easy.contra.pre) ,nanmean(pHit.easy.contra.post));
    ext_sigline([y1,y2],pFA.easy.contra.pvalue,[],ymax -0.2,'x', Color.Contra); hold on;
end
if pHit.easy.contra.pvalue < 0.05
    y1 = nanmean(pHit.easy.contra.pre) ;
    y2 = nanmean(pHit.easy.contra.post);
    ymax = max(nanmean(pFA.easy.contra.pre) ,nanmean(pFA.easy.contra.post));
    ext_sigline([y1,y2],pHit.easy.contra.pvalue,[],ymax +0.2,'y',Color.Contra); hold on;
end


ha(indPos) = subplot(1,2,2);
for i = 1: length(d_prime.easy.contra.pre)
    plot([d_prime.easy.contra.pre(i),d_prime.easy.contra.post(i)],[criterion.easy.contra.pre(i),criterion.easy.contra.post(i)], 'o','color',[Color.Contra, 0.4] , 'MarkerSize',MarkSize_CritDpr_small,'markerfacecolor',[1 1 1 ]); hold on;
    plot([d_prime.easy.ipsi.pre(i),d_prime.easy.ipsi.post(i)], [-criterion.easy.ipsi.pre(i),-criterion.easy.ipsi.post(i)], 'o','color',[Color.Ipsi, 0.4] , 'MarkerSize',MarkSize_CritDpr_small,'markerfacecolor',[1 1 1 ]); hold on;% reverse direction of criterion for ipsi
    plot(d_prime.easy.contra.post(i),criterion.easy.contra.post(i), 'o','color',Color.Contra , 'MarkerSize',MarkSize_CritDpr_small,'markerfacecolor',Color.Contra); hold on;
    plot(d_prime.easy.ipsi.post(i),-criterion.easy.ipsi.post(i), 'o','color',Color.Ipsi, 'MarkerSize',MarkSize_CritDpr_small,'markerfacecolor',Color.Ipsi); hold on;% reverse direction of criterion for ipsi
end
plot([nanmean(d_prime.easy.ipsi.pre),nanmean(d_prime.easy.ipsi.post)],[ -nanmean(criterion.easy.ipsi.pre) ,-nanmean(criterion.easy.ipsi.post)], 'o-','color',Color.Ipsi , 'MarkerSize',MarkSize_CritDpr,'markerfacecolor',[1 1 1 ],'LineWidth',LineWidthSize ); hold on;% reverse direction of criterion for ipsi
plot(nanmean(d_prime.easy.ipsi.post),-nanmean(criterion.easy.ipsi.post), 'o','color',Color.Ipsi , 'MarkerSize',MarkSize_CritDpr,'markerfacecolor',Color.Ipsi); hold on;% reverse direction of criterion for ipsi
plot([nanmean(d_prime.easy.contra.pre),nanmean(d_prime.easy.contra.post)],[nanmean(criterion.easy.contra.pre) ,nanmean(criterion.easy.contra.post)], 'o-','color',Color.Contra, 'MarkerSize',MarkSize_CritDpr,'markerfacecolor',[1 1 1 ],'LineWidth',LineWidthSize ); hold on;
plot(nanmean(d_prime.easy.contra.post),nanmean(criterion.easy.contra.post), 'o','color',Color.Contra , 'MarkerSize',MarkSize_CritDpr,'markerfacecolor',Color.Contra); hold on;

if d_prime.easy.ipsi.pvalue < 0.05
    ymax = max(-1*nanmean(criterion.easy.ipsi.pre) ,-1*nanmean([criterion.easy.ipsi.post]));
    ext_sigline([nanmean([d_prime.easy.ipsi.pre]),nanmean([d_prime.easy.ipsi.post])],d_prime.easy.ipsi.pvalue,[ ],(ymax + 0.5 ),'x', Color.Ipsi); hold on;
end
if criterion.easy.ipsi.pvalue < 0.05
    y1 = -1*nanmean(criterion.easy.ipsi.pre) ;
    y2 = -1*nanmean([criterion.easy.ipsi.post]);
    ymax = max(nanmean(d_prime.easy.ipsi.pre) ,nanmean([d_prime.easy.ipsi.post])) ;
    ext_sigline([y1,y2],criterion.easy.ipsi.pvalue,[],ymax+ 0.2,'y', Color.Ipsi); hold on;
end


if d_prime.easy.contra.pvalue < 0.05
    ymax = max(nanmean(criterion.easy.contra.pre) ,nanmean([criterion.easy.contra.post]));
    ext_sigline([nanmean([d_prime.easy.contra.pre]),nanmean([d_prime.easy.contra.post])],d_prime.easy.contra.pvalue,[ ],(ymax- 0.8 ),'x', Color.Contra); hold on;
end
if criterion.easy.contra.pvalue < 0.05
    y1 = nanmean(criterion.easy.contra.pre) ;
    y2 = nanmean([criterion.easy.contra.post]);
    ymax = max(nanmean(d_prime.easy.contra.pre) ,nanmean([d_prime.easy.contra.post])) ;
    ext_sigline([y1,y2],criterion.easy.contra.pvalue,[],ymax+ 0.1,'y', Color.Contra); hold on;
end
axis square
xlabel('sensitivity','fontsize',fs,'fontweight','b', 'Interpreter', 'none')
ylabel('criterion','fontsize',fs,'fontweight','b', 'Interpreter', 'none')
set(gca,'ylim',[-3 3],'xlim',xlim_SDT_easy,'fontsize',fs)

text(xlim_SDT_easy(1) + 0.1,-2.8, 'more Contra (ipsi:NoGo, contra:Go)', 'Color', 'k','fontsize',18)
text(xlim_SDT_easy(1) + 0.1,2.8, 'less Contra (ipsi:Go, contra:NoGo)', 'Color', 'k','fontsize',18)

%% save the graphs
if SaveGraph
    h= figure(1);
    %savefig(h, [path_SaveFig ,filesep,'fig',filesep,  'DoubleStimuli_Crit_Dprime_' monkey,'_' experiment ,'_' folder '.fig'])
    print(h,[path_SaveFig ,filesep,  ['png', ControlFolder]  ,filesep,  'TarDistStimuli_easy_HR_FAR_Dprime_Criterion' monkey,'.png'], '-dpng') %'_' experiment,'_' folder
    set(h,'Renderer','Painters');
    set(h,'PaperPositionMode','auto')
    compl_filename =  [path_SaveFig , ['ai', ControlFolder]  ,filesep,  'TarDistStimuli_easy_HR_FAR_Dprime_Criterion' monkey,'_' ,'.ai'] ; %experiment ,'_', folder
    print(h,'-depsc',compl_filename);
    close all;
end

end


