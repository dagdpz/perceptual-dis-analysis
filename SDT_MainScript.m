

%% Inactivation %%%%%%%%% MANUSCRIPT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
monkey = 'Curius'; 
BrainArea = 'dPul';
folder = 'Inactivation_20190729_20190801_20190809_20190814_20190820_20190905_20190913';
dataset = 'Behavior_Select_RT_20201010-2246'; 
folder_baseline = 'Baseline_20190912_20190910_20190903_20190815_20190808_20190806_20190802';
dataset_baseline = 'Behavior_Select_RT_20200405-1703';
experiment = 'Inactivation'; 
task = 'distractorColorTask'; 
path_SaveFig = ['Y:\Projects\Pulv_distractor_spatial_choice\Inactivation\Results',filesep, BrainArea,  filesep, experiment, filesep, monkey,filesep, task,filesep]; 
path_SaveData = ['Y:\Projects\Pulv_distractor_spatial_choice\Inactivation\Data'];
path_getData = ['Y:\Projects\Pulv_Inac_ECG_respiration\Figures\'];

monkey = 'Cornelius'; 
BrainArea = 'dPul';
folder = 'Inactivation_20190124_20190129_20190201_20190207_20190214_20190228_20190314';
dataset = 'Behavior_Select_RT_20201005-2112'; 
folder_baseline = 'Baseline_20190121_20190131_20190213_20190216_20190227_20190304_20190313'; %'Baseline_20190131_20190213_20190216_20190227_20190304_20190313';
dataset_baseline = 'Behavior_Select_RT_20201005-2122'; %'Behavior_20191217-1505.mat';
experiment = 'Inactivation'; 
path_SaveFig = ['Y:\Projects\Pulv_distractor_spatial_choice\Inactivation\Results',filesep, BrainArea, filesep, experiment, filesep, monkey,filesep]; 
path_SaveData = ['Y:\Projects\Pulv_distractor_spatial_choice\Inactivation\Data'];
path_getData = ['Y:\Projects\Pulv_Inac_ECG_respiration\Figures\'];



%% Cmb data sets to compute Criterion, Dprime (save table), Create Figures, Compute ttest (table with t-value,p-values)
ResponseBias_DisplayColor = 0;
SDT_singleStimuli(monkey, folder,dataset,folder_baseline,dataset_baseline, experiment,path_SaveFig,path_SaveData,path_getData,ResponseBias_DisplayColor)
SDT_doubleSameStimuli_2H(monkey, folder,dataset,folder_baseline,dataset_baseline,experiment,path_SaveFig,path_SaveData,path_getData, ResponseBias_DisplayColor)
SDT_TargetDistractorStimuli_2H(monkey, folder,dataset,folder_baseline,dataset_baseline,experiment,path_SaveFig,path_SaveData,path_getData, ResponseBias_DisplayColor)

