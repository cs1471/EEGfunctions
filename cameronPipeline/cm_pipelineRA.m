%% Cameron's main pipeline for EEG-RA anaylsis
% Fall 2015 rotation in MAXLAB
% GU NetID: ccm98
% Date created: 10/5/2015

%% General structure:
% 1) SET-UP: change directory, load in/create important variables
% 2) PREPROCESSING (EITHER WITH OR WITHOUT ICA - SWITCH STATEMENTS BASED ON RUN.analysisType)
% 3) CREATE TRIALSTRUCT TO SORT TRIALS AND CONDITION (A LOGICAL VECTOR BASED ON TRIAL PARAMETERS)
% 4) CONVERT FROM EEGLAB TO FIELDTRIP
% 5) PLOTTING!!
% 6) CLUSTERING (AND BEYOND!)

%% 1. SET-UP

% clear workspace/path; create RUN
clear all
restoredefaultpath
global RUN

% set data path
dataPath = 'C:/Users/Sikoya/Desktop/vwfaEEG/';
cd(dataPath)
RUN.dataPath = dataPath;
addpath([dataPath 'data/'])
addpath([dataPath 'functions/'])

% load any other datasets with info about subjects)
% do this manually for now
subs = {'subj776', 'subj794', 'subj801', 'subj859',...
        'subj880', 'subj882', 'subj893', 'subj895',...
        'subj897', 'subj898', 'subj900'};
RUN.subIDs = subs;
RUN.subNums = { '776', '794', '801',...
                '859', '880', '882',...
                '893', '895', '897',...
                '898', '900'};
RUN.filenames.preprocClara = strcat(subs,...
    'FT_bothStim_correctOnly_500hz0.1hz-30hz_fromContinuous.mat');
RUN.filenames.eegFilename = { 'S776_20140220_vwfa001.raw', 'S794_20140403_vwfa001.raw', 'S801_20140602_vwfa001.raw',...
    'S859_20140428_vwfa001.raw', 'S880_20140624_vwfa001.raw', 'S882_20140717_vwfa001.raw',...
    'S893_20140528_vwfa001.raw', 'S895_20140610_vwfa001.raw', 'S897_20140626_vwfa001.raw',...
    'S898_20140806_vwfa001.raw', 'S900_20140717_vwfa001.raw' };
RUN.filenames.exptFilename = { 's776.2408.mat', 's794.2408.mat', 's801.2408.mat',...
    's859.2408.mat', 's880_vwfa.2408.mat', 's882_vwfa.2408.mat',...
    's893.2408.mat', 's895.2408.mat', 's897_vwfa.2136.mat',...
    's898.2408.mat', 's900_vwfa.2408.mat' };

% make template
load([dataPath, 'Hydrocel_GSN_128_1.0_TRIM_mod.sfp']);
RUN.template = Hydrocel_GSN_128_1_0_TRIM_mod;
RUN.elecLoc.elecpos = RUN.template(:,2:4);
RUN.elecLoc.chanpos = RUN.template(:,2:4);
RUN.elecLoc.label = strread(num2str(RUN.template(:,1)'),'%s')'; % convert to cell array of char
RUN.elecLoc.unit = 'dm';

%save RUN, add some other stuff to the path
% save([RUN.dataPath, 'data/', 'RUN.mat'], 'RUN');
% addpath([dataPath 'fieldtrip-20141113/']);
% addpath([dataPath 'tightfig/']);
% addpath([dataPath 'altmany-export_fig-d966721/']);

%% 1A. Do some of Clara's pipeline
% concatenate all the .raw files into one EEGLAB .set per subject
for iSub = 1:length(RUN.subIDs)
    cm_autocrank(iSub)
end


%% 2. PREPROCESSING (EITHER WITH OR WITHOUT ICA - SWITCH STATEMENTS BASED ON RUN.analysisType)

addpath([RUN.dataPath 'eeglab13_4_4b/']);

% use RUN to set variables for preprocessing
RUN.forICA.filterSettings = {0.5, 70};
RUN.forICA.filterType = 'causal'; % can be either 'causal' or 'non-causal'
RUN.forICA.epochLength = [-0.9 0.7];
RUN.forICA.interp = RUN.maxlabPreproc.interp
RUN.forICA.reref = 'no';
RUN.forICA.noiseThresh = [-125 1500];
RUN.forICA.baseline = [-0.9 -0.7];
% trying to speed up ICA
RUN.forICA.resample = 'yes'; % can be either 'yes' or 'no'
RUN.forICA.newSrate = 250; % set to 550 if not resampling


% new set of RUN parameters for Maxlab
RUN.analysisType = 'maxlab'; %or set to 'cameron'
RUN.maxlabPreproc.filterSettings = {0.1, 30};
RUN.maxlabPreproc.filterType = 'causal'; % can be either 'causal' or 'non-causal'
RUN.maxlabPreproc.epochLength = [-0.9 0.7];
RUN.maxlabPreproc.interp = RUN.forICA.interp;
RUN.maxlabPreproc.reref = 'no';
RUN.maxlabPreproc.noiseThresh = [-125 1500];
RUN.maxlabPreproc.baseline = [-0.9 -0.7];
RUN.maxlabPreproc.artDetection = 'thresh';
RUN.maxlabPreproc.resample = 'no';

eeglab

% loop through each subject
for iSub = 1:length(RUN.subNums)
    
    currentSub = RUN.subNums{iSub}
    
    % load in each subject's data
    EEG = pop_loadset('filename', [currentSub, '_continuous.set'], 'filepath', [RUN.dataPath, 'data/', currentSub, '/']);
    
    % call my function here
    switch RUN.analysisType
        case 'maxlab'
            cm_maxlabPreproc(EEG,currentSub,iSub);
        case 'cameron'
            cm_vwfaICA(EEG,currentSub,iSub);
    end
    
    
end

save([RUN.dataPath, 'data/RUN.mat'], 'RUN');


%% 2A. REGULAR PREPROCESSING (IF ONLY ICA WAS DONE IN PREVIOUS SECTION)
% this will have similar steps being applied to raw data
% also need to remove components

addpath([RUN.dataPath 'eeglab13_4_4b/']);
load([dataPath 'data/RUN.mat'])

% use RUN to set variables for preprocessing
RUN.preproc.filterSettings = {1, 40};
RUN.preproc.filterType = 'causal'; % can be either 'causal' or 'non-causal'
RUN.preproc.stimuli = {'DIN2'};
RUN.preproc.epochLength = [-0.9 0.7];
RUN.preproc.interp = RUN.forICA.interp;
RUN.preproc.reref = 'no'; % can do 'average' or 'average_mastoid' too

RUN.preproc.resample = 'yes'; % can be either 'yes' or 'no'
RUN.preproc.newSrate = 250; % set to 550 if not resampling
RUN.preproc.noiseThresh = [-125 1500];
RUN.preproc.baseline = [-0.9 -0.7];

RUN.preproc.rejComps = 'yes_first'; % 'yes_first' if ICA has been changed
RUN.preproc.artDetection = 'thresh'; % can also be 'thresh'
RUN.preproc.artCriteria = [-100 -100]; % % can be trend and R, or thresh limits
RUN.arfEpochs = [];
RUN.numTrials = [];

% only do preprocessing here if using 'cameron' as RUN.analysisType
switch RUN.analysisType
    case 'cameron'
        
        % loop through each subject
        for iSub = 1:length(RUN.subNums)
            
            currentSub = RUN.subNums{iSub}
            
            % load in each subject's raw data
            EEG = pop_loadset([currentSub, '_continuous.set'], [RUN.dataPath, 'data/', currentSub, '/']);
            
            % load in each subject's ICAstruct
            load([RUN.dataPath, 'data/', currentSub, '/ICAstruct.mat'])
            
            % to pass in ICA stuff into preprocessing function
            RUN.icaweights{iSub} = ICAstruct.icaweights;
            RUN.icasphere{iSub} = ICAstruct.icasphere;
            RUN.icawinv{iSub} = ICAstruct.icawinv;
            RUN.icachansind{iSub} = ICAstruct.icachansind;
            
            % call my function here
            cm_vwfaPreproc(EEG,currentSub,iSub);
            
            
        end
        
end

save([RUN.dataPath, 'data/RUN.mat'], 'RUN');




%% 3. CREATE TRIALSTRUCT TO SORT TRIALS AND CONDITION (A LOGICAL VECTOR BASED ON TRIAL PARAMETERS)

% MAKE TRIALSTRUCT

load([dataPath 'data/RUN.mat'])

for iSub = 1:length(RUN.subNums)
    
    % set up subject number and the filename of their data
    currentSub = RUN.subNums{iSub}

    % add eeglab to path, remove fieldtrip from path
    rmpath([RUN.dataPath, 'fieldtrip-20151012/']);
    addpath([RUN.dataPath 'eeglab13_4_4b/']);
    
    % actually loading in preprocessed data
    switch RUN.analysisType
        case 'maxlab'
%             EEG = pop_loadset([currentSub, '_maxlabPreproc.set'], [RUN.dataPath, 'data/', currentSub, '/']);
            EEG = pop_loadset([currentSub, '_maxlabPreprocCausal.set'], [RUN.dataPath, 'data/', currentSub, '/']);
        case 'cameron'
            EEG = pop_loadset([currentSub, '_preproc.set'], [RUN.dataPath, 'data/', currentSub, '/']);
    end
    
    % create empty trialStruct
    trialStruct = dataset();

    % add columns to trialStruct
    trialStruct = cm_vwfaTrialStruct( RUN, EEG, iSub );
    
    trialStruct.isArtifact = EEG.reject.rejthresh';

    % add number of trials used and rejected to allSubSummary
    sizeStruct = size(trialStruct);
    numTrials = sizeStruct(1);
    numRejected = sum(trialStruct.isArtifact);
    
%     switch RUN.analysisType
%         case 'maxlab'
%             RUN.maxlabNumTrials{iSub} = numTrials - numRejected;
%         case 'cameron'
%             RUN.numTrials{iSub} = numTrials - numRejected;
%     end
    
    % rejected trial indexes (incorrect resp and/or artifact)
    badTrials = zeros(1,length(numTrials));
    for iTrial = 1:numTrials
        if trialStruct.isArtifact(iTrial) == 1 || trialStruct.isCorrect(iTrial) == 0
            badTrials = [badTrials iTrial];
        end
    end
    
    EEG.badTrials = badTrials;
    
    % save trialStruct to RUN and save trialStruct separately too
    switch RUN.analysisType
        case 'maxlab'
            RUN.maxlabTrialStruct{1,iSub} = trialStruct;
%             save([RUN.dataPath, 'data/', currentSub, '/', 'maxlabTrialStruct.mat'], 'trialStruct');
            save([RUN.dataPath, 'data/', currentSub, '/', 'maxlabTrialStructCausal.mat'], 'trialStruct');
        case 'cameron'
            RUN.trialStruct{1,iSub} = trialStruct;
            save([RUN.dataPath, 'data/', currentSub, '/', 'trialStruct.mat'], 'trialStruct');
    end

    % create conditions and convert to fieldtrip
    condition = dataset();
    condition = cm_vwfaCreateConditions(trialStruct, EEG);
    
    % save data in different forms for fieldtrip
    data = eeglab2fieldtrip(EEG,'preprocessing');
    data.condition = condition;
    data.trialStruct = trialStruct;
    
    switch RUN.analysisType
        case 'maxlab'
%             save([RUN.dataPath, 'data/', currentSub, '/', 'maxlabPreprocessedData.mat'], 'data'); 
            save([RUN.dataPath, 'data/', currentSub, '/', 'maxlabPreprocessedDataCausal.mat'], 'data'); 
        case 'cameron'
            save([RUN.dataPath, 'data/', currentSub, '/', 'preprocessedData.mat'], 'data'); 
    end
    
    
    % clear EEG and data, just in case
    clear EEG
    clear data

end

save([RUN.dataPath, 'data/', 'RUN.mat'], 'RUN');




%% 4. CONVERT FROM EEGLAB TO FIELDTRIP

addpath([RUN.dataPath, 'fieldtrip-20151012/']);
rmpath([RUN.dataPath 'eeglab13_4_4b/']);

load([dataPath, 'data\', 'RUN.mat'])

ga = [];
gaCausal = [];
timelock = [];
conditions = [];

% averaging by condition
for iSub = 1:length(RUN.subNums)
    
    % set up subject number and the filename of their data
    currentSub = RUN.subNums{iSub}
    load([RUN.dataPath, 'data/', currentSub, '/', 'maxlabPreprocessedDataCausal.mat']);
         
    fn = fieldnames(data.condition);
    for iConditions= 1:length(fn);
        cfg = [];
        cfg.trials = condition.(fn{iConditions});
        if sum(cfg.trials) ~= 0 % so you don't end up averaging together 0 trials
            dataOnly = rmfield(data,{'trialStruct', 'condition'});
            timelock.(fn{iConditions}){iSub} = ft_timelockanalysis(cfg,dataOnly);
        end
    end
    
%     for iConditions= 1:length(fn);
%         data = eval(fn{iConditions});
%         cfg = [];
%         [dummy, numTrials] = size(data.time);
%         cfg.trials = 'all'; % or: = ones(1,numTrials);
%         timelock.(fn{iConditions}){iSub} = ft_timelockanalysis(cfg,data);
% %         timelock.(fn{iConditions}){iSub} = data;
%     end
  
end

% baseline timelock before making ga
cfgbase = [];
cfgbase.baseline = [-.2 0];
cfgbase.fields = {'avg', 'individual'};
timelock = cm_vwfaBaseline( cfgbase, timelock );
        
% grand average by condition
fn = fieldnames(timelock);
for iConditions = 1:length(fn)
    cfg = [];
    cfg.keepindividual = 'yes';

    gaCausal.(fn{iConditions}) = ft_timelockgrandaverage(cfg, timelock.(fn{iConditions}){:});
    gaCausal.(fn{iConditions}).elec = timelock.(fn{iConditions}){1,1}.elec ;

    % since this is the first time running this (presumably) one needs to
    % fix the channel locations and save it appropriately
    gaCausal = cm_fixChanLocs(gaCausal, fn, iConditions);
    
end


% set up contrasts
contrasts = [];
contrasts.ortho = ga.DWSC;
contrasts.ortho.individual = ga.DWSCcorrect.individual - ga.SWSCcorrect.individual;
contrasts.catsem = ga.DWDC;
contrasts.catsem.individual = ga.DWDCcorrect.individual - ga.DWSCcorrect.individual;
contrasts.subsem = ga.DWDSC;
contrasts.subsem.individual = ga.DWDSCcorrect.individual - ga.DWSCcorrect.individual;
contrasts.DWSCminusSWSC = ga.DWSC;
contrasts.DWSCminusSWSC.individual = ga.DWSCcorrect.individual - ga.SWSCcorrect.individual;
contrasts.DWDCminusDWSC = ga.DWDC;
contrasts.DWDCminusDWSC.individual = ga.DWDCcorrect.individual - ga.DWSCcorrect.individual;
contrasts.DWDSCminusDWSC = ga.DWDC;
contrasts.DWDSCminusDWSC.individual = ga.DWDSCcorrect.individual - ga.DWSCcorrect.individual;

% add field to RUN to define types/locations of channels
RUN.channels = cm_defineChannels(RUN.elecLoc.elecpos, gaCausal.DWDC.label);

% create an 'avg' struct in each ga and contrasts field
gaCausal = cm_createAvg(gaCausal);
contrasts = cm_createAvg(contrasts);

% before moving on, SAVE EVERYTHING
save([RUN.dataPath, 'data/timelock.mat'], 'timelock');
% save([RUN.dataPath, 'data/ga.mat'], 'ga');
save([RUN.dataPath, 'data/gaCausal.mat'], 'gaCausal');
save([RUN.dataPath, 'data/contrasts.mat'], 'contrasts');
save([RUN.dataPath, 'data/RUN.mat'], 'RUN');


%% 5. PLOTTING!!

addpath([dataPath 'data/'])
addpath([RUN.dataPath, 'fieldtrip-20151012/']);
rmpath([RUN.dataPath 'eeglab13_4_4b/']);

load([dataPath 'data/RUN.mat']) % load RUN
load([dataPath 'data/ga.mat']) % load ga
load([dataPath 'data/gaCausal.mat']) % load ga
load([dataPath 'data/contrasts.mat']) % load contrasts
load([dataPath 'data/timelock.mat'])

cfg = [];
cfg.parameter = 'avg';
cfg.hlim = [-0.2 0.8];
if ~isfield(cfg,'hlim')
    cfg.baseline = [-0.9 -0.7];
else
    cfg.baseline = [-0.2 0];
end
cfg.channel = 'all';

% compare SWSC for the two filter settings
figure;
ft_multiplotER(cfg, ga.SWSCcorrect, gaCausal.SWSCcorrect)

% now look at the four conditions for the gaCausal filter
figure;
ft_multiplotER(cfg, gaCausal.SWSCcorrect, gaCausal.DWSCcorrect,...
    gaCausal.DWDCcorrect, gaCausal.DWDSCcorrect)


figure;
ft_multiplotER(cfg, ga.SWSC, ga.DWSC, ga.DWDC, ga.DWDSC)
figure;
ft_multiplotER(cfg, ga.DWSCcorrect, ga.DWDCcorrect, ga.DWDSCcorrect)

% multiplot for each contrast
% figure;
% ft_multiplotER(cfg, contrasts.ortho)
% figure;
% ft_multiplotER(cfg, contrasts.catsem)
% figure;
% ft_multiplotER(cfg, contrasts.subsem)

%% plot a single ERP for each condition
% just to see what's interesting
cfg = [];
cfg.parameter = 'avg';
cfg.baseline = [-0.9 -0.7];
conds = {'DWSCcorrect', 'DWDCcorrect', 'DWDSCcorrect'};
titles = {'DWSC', 'DWDC', 'DWDSC'};

figure;

for iCond = 1:length(conds)

    subplot(3,1,iCond);ft_singleplotER(cfg, ga.(conds{iCond}))
    hold on;
    plot([-0.7 -0.7],[-25 25],'k-') % stim 1 onset
    plot([-0.4 -0.4],[-25 25],'k-') % stim 1 offset
    plot([0 0],[-25 25],'k-') % stim 2 onset
    plot([0.3 0.3],[-25 25],'k-') % stim 2 offset
    plot([-0.9 0.7],[0 0],'k-') % x-axis (potential = 0 uV)
    hold off;
    xAxisLabelsOLD = {'-200', '-100', '0', '100', '200', '300',...
        '400', '500', '600', '700', '800', '900',...
        '1000', '1100', '1200', '1300', '1400'};
    xAxisLabels = {'-900', '-800', '-700', '-600', '-500', '-400',...
        '-300', '-200', '-100', '0', '100', '200',...
        '300', '400', '500', '600', '700'};
    set(gca,'Xtick',-0.9:0.1:0.7)
    set(gca,'XtickLabel',-0.9:0.1:0.7)
    set(gca,'XTickLabel',xAxisLabels)
    title(titles{iCond})

end



% looking only at stim 2
cfg = [];
cfg.parameter = 'avg';
cfg.baseline = [-0.2 0];
cfg.xlim = [-0.2 0.7];
cfg.channel = RUN.channels.anteriorLeft;
conds = {'DWSCcorrect', 'DWDCcorrect', 'DWDSCcorrect'};
titles = {'DWSC', 'DWDC', 'DWDSC'};


figure;

% for iCond = 1:length(conds)

%     subplot(3,1,iCond);
    ft_singleplotER(cfg, ga.(conds{1}), ga.(conds{2}), ga.(conds{3}))
    hold on;
%     plot([-0.7 -0.7],[-25 25],'k-') % stim 1 onset
%     plot([-0.4 -0.4],[-25 25],'k-') % stim 1 offset
    plot([0 0],[-25 25],'k-') % stim 2 onset
    plot([0.3 0.3],[-25 25],'k-') % stim 2 offset
    plot([-0.9 0.7],[0 0],'k-') % x-axis (potential = 0 uV)
    hold off;
    xAxisLabelsOLD = {'-200', '-100', '0', '100', '200', '300',...
        '400', '500', '600', '700', '800', '900',...
        '1000', '1100', '1200', '1300', '1400'};
    xAxisLabels = {'-200', '-100', '0', '100', '200',...
        '300', '400', '500', '600', '700'};
    set(gca,'Xtick',-0.2:0.1:0.7)
    set(gca,'XtickLabel',-0.2:0.1:0.7)
    set(gca,'XTickLabel',xAxisLabels)
    title(titles{iCond})

% end


% plot a single ERP for each contrast
% just to see what's interesting
cfg = [];
cfg.parameter = 'avg';
cfg.baseline = [-0.9 -0.7];
conts = {'ortho', 'catsem', 'subsem'};
titles = {'Orthographic (DWSC - SWSC)', 'Semantic 1 (DWDC - DWSC)',...
    'Semantic 2 (DWDSC - DWSC)'};
figure;

for iCond = 1:length(conts)

    subplot(3,1,iCond);ft_singleplotER(cfg, contrasts.(conts{iCond}))
    hold on;
    plot([-0.7 -0.7],[-25 25],'k-') % stim 1 onset
    plot([-0.4 -0.4],[-25 25],'k-') % stim 1 offset
    plot([0 0],[-25 25],'k-') % stim 2 onset
    plot([0.3 0.3],[-25 25],'k-') % stim 2 offset
    plot([-0.9 0.7],[0 0],'k-') % x-axis (potential = 0 uV)
    hold off;
    xAxisLabelsOLD = {'-200', '-100', '0', '100', '200', '300',...
        '400', '500', '600', '700', '800', '900',...
        '1000', '1100', '1200', '1300', '1400'};
    xAxisLabels = {'-900', '-800', '-700', '-600', '-500', '-400',...
        '-300', '-200', '-100', '0', '100', '200',...
        '300', '400', '500', '600', '700'};
    set(gca,'Xtick',-0.9:0.1:0.7)
    set(gca,'XtickLabel',-0.9:0.1:0.7)
    set(gca,'XTickLabel',xAxisLabels)
    title(titles{iCond})

end


% topos for each condition
cfg = [];
cfg.zlim = [-3 3];
cfg.layout = RUN.template;
cfg.parameter = 'avg';
cfg.interactive='no';
cfg.comment = 'no';
cfg.colorbar = 'no';
cfg.baseline = [-0.2 0];

conds = {'SWSCcorrect', 'DWSCcorrect', 'DWDCcorrect', 'DWDSCcorrect'};
numSteps = 5;
steps = 0.1;
startingTime = 0;

cm_vwfaPlotTopo( cfg, gaCausal, conds, numSteps, steps, startingTime )



% difference wave topos
load([dataPath 'data/contrasts.mat']) % load contrasts
load([dataPath 'data/ga.mat']) % load ga
load([dataPath 'data/RUN.mat']) % load RUN
% actually plotting
cfg = [];
cfg.zlim = [-3 3];
cfg.layout = RUN.template;
cfg.parameter = 'avg';
cfg.interactive='no';
cfg.comment = 'no';
cfg.colorbar = 'no';
cfg.baseline = [-0.2 0];
conds = {'DWSCminusSWSC', 'DWDCminusDWSC', 'DWDSCminusDWSC'};
numSteps = 3;
steps = 0.1;
startingTime = 0.2;

cm_vwfaPlotTopo( cfg, contrasts, conds, numSteps, steps, startingTime )




%% 6. CLUSTERING (AND BEYOND!)

% % first define cfg.neighbours
% cfg = [];
% cfg.method = 'template';
% cfg.template = 'Hydrocel_GSN_128_1.0_TRIM.sfp';
% 
% 
% % actually make neighbours
% cfg.neighbours = ft_prepare_neighbours(cfg, data);



% try this with the MEG demo data first
cfg_neighb        = [];
cfg_neighb.method = 'distance';
neighbours        = ft_prepare_neighbours(cfg_neighb, allsubjFC{1});

cfg = [];
cfg.channel = {'MEG'};
cfg.latency = [0 1];

cfg.method = 'montecarlo';
cfg.statistic = 'depsamplesT';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 2;
cfg.neighbours = neighbours;  % same as defined for the between-trials experiment
cfg.tail = 0;
cfg.clustertail = 0;
cfg.alpha = 0.025;
cfg.numrandomization = 500;


subj = 10;
design = zeros(2,2*subj);
for i = 1:subj
  design(1,i) = i;
end
for i = 1:subj
  design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;

cfg.design = design;
cfg.uvar  = 1;
cfg.ivar  = 2;

[stat] = ft_timelockstatistics(cfg, allsubjFIC{:}, allsubjFC{:})

% plotting the results
% load individual subject data
load('ERF_orig');
% calculate the grand average for each condition
cfg = [];
cfg.channel   = 'all';
cfg.latency   = 'all';
cfg.parameter = 'avg';
GA_FC         = ft_timelockgrandaverage(cfg,allsubjFC{:});  
GA_FIC        = ft_timelockgrandaverage(cfg,allsubjFIC{:});
% "{:}" means to use data from all elements of the variable

cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
GA_FICvsFC = ft_math(cfg,GA_FIC,GA_FC);

figure;  
% define parameters for plotting
timestep = 0.05;      %(in seconds)
sampling_rate = allsubjFC{1}.fsample;
sample_count = length(stat.time);
j = [0:timestep:1];   % Temporal endpoints (in seconds) of the ERP average computed in each subplot
m = [1:timestep*sampling_rate:sample_count];  % temporal endpoints in MEEG samples
% get relevant (significant) values
pos_cluster_pvals = [stat.posclusters(:).prob];


% In case you have downloaded and loaded the data, ensure stat.cfg.alpha exists:
if ~isfield(stat.cfg,'alpha'); stat.cfg.alpha = 0.025; end; % stat.cfg.alpha was moved as the downloaded data was processed by an additional fieldtrip function to anonymize the data.
 
pos_signif_clust = find(pos_cluster_pvals < stat.cfg.alpha);
pos = ismember(stat.posclusterslabelmat, pos_signif_clust);

% plot
for k = 1:20;
     subplot(4,5,k);   
     cfg = [];   
     cfg.xlim=[j(k) j(k+1)];   
     cfg.zlim = [-1.0e-13 1.0e-13];   
     pos_int = all(pos(:, m(k):m(k+1)), 2);
     cfg.highlight = 'on';
     cfg.highlightchannel = find(pos_int);       
     cfg.comment = 'xlim';   
     cfg.commentpos = 'title';   
     cfg.layout = 'CTF151.lay';
     ft_topoplotER(cfg, GA_FICvsFC);
end  




%% doing this for real using timelock data

% try this with the MEG demo data first
cfg_neighb        = [];
cfg_neighb.method = 'distance';
neighbours        = ft_prepare_neighbours(cfg_neighb, timelock.DWDCcorrect{1});

cfg = [];
cfg.channel = RUN.channels.anteriorLeft;
cfg.latency = [0.2 0.3];

cfg.method = 'montecarlo';
cfg.statistic = 'depsamplesT';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 2;
cfg.neighbours = neighbours;  % same as defined for the between-trials experiment
cfg.tail = 0;
cfg.clustertail = 0;
cfg.alpha = 0.05;
cfg.numrandomization = 500;

% cfg.method = 'stats';

% create a design matrix
subj = 11;
design = zeros(2,2*subj);
for i = 1:subj
  design(1,i) = i;
end
for i = 1:subj
  design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;

cfg.design = design;
cfg.uvar  = 1;
cfg.ivar  = 2;

[stat1] = ft_timelockstatistics(cfg, timelock.DWSCcorrect{:}, timelock.SWSCcorrect{:})
[stat2] = ft_timelockstatistics(cfg, timelock.DWDCcorrect{:}, timelock.SWSCcorrect{:})
[stat3] = ft_timelockstatistics(cfg, timelock.DWDSCcorrect{:}, timelock.SWSCcorrect{:})
[stat4] = ft_timelockstatistics(cfg, timelock.DWDCcorrect{:}, timelock.DWSCcorrect{:})
[stat5] = ft_timelockstatistics(cfg, timelock.DWDSCcorrect{:}, timelock.DWSCcorrect{:})
statAll = {stat1, stat2, stat3, stat4, stat5};

% try to compare all of the conditions in one step, doesn't work
subj = 11;
cfg_all = cfg;
design = zeros(2,4*subj);
for i = 1:subj
  design(1,i) = i;
  design(1,subj+i) = i;
  design(1,2*subj+i) = i;
  design(1,3*subj+i) = i;
end

design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;
design(2,2*subj+1:3*subj) = 3;
design(2,3*subj+1:4*subj) = 4;

cfg_all.design = design;
cfg_all.uvar  = 1;
cfg_all.ivar  = 2;

[statNew] = ft_timelockstatistics(cfg_all, timelock.DWDSCcorrect{:},...
    timelock.DWDCcorrect{:}, timelock.DWSCcorrect{:}, timelock.SWSCcorrect{:})

% plotting the results
% load individual subject data
% calculate the grand average for each condition
cfg = [];
cfg.channel   = 'all';
cfg.latency   = 'all';
cfg.parameter = 'avg';
GA_SWSC        = ft_timelockgrandaverage(cfg,timelock.SWSCcorrect{:});
GA_DWSC        = ft_timelockgrandaverage(cfg,timelock.DWSCcorrect{:});
GA_DWDC        = ft_timelockgrandaverage(cfg,timelock.DWDCcorrect{:});
GA_DWDSC       = ft_timelockgrandaverage(cfg,timelock.DWDSCcorrect{:});

cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
GA_DWSCvsSWSC = ft_math(cfg,GA_DWSC,GA_SWSC);
GA_DWDCvsSWSC = ft_math(cfg,GA_DWDC,GA_SWSC);
GA_DWDSCvsSWSC = ft_math(cfg,GA_DWDSC,GA_SWSC);
GA_DWDCvsDWSC = ft_math(cfg,GA_DWDC,GA_DWSC);
GA_DWDsCvsDWSC = ft_math(cfg,GA_DWDSC,GA_DWSC);

figure;  
% define parameters for plotting
timestep = 0.05;      %(in seconds)
sampling_rate = 500;
sample_count = length(stat.time);
j = [0:timestep:1];   % Temporal endpoints (in seconds) of the ERP average computed in each subplot
m = [1:timestep*sampling_rate:sample_count];  % temporal endpoints in MEEG samples

% get relevant (significant) values
stat2check = {'stat1', 'stat2', 'stat3', 'stat4', 'stat5'};
field2check = {'posclusters', 'negclusters'};
pvals2check = {'positive_cluster_pvals', 'negative_cluster_pvals'};

cluster_pvals = {};

for s = 1:length(statAll)
    for f = 1:length(field2check)
        
        if isfield(statAll{s},{field2check{f}}) && ~isempty(statAll{s}.(field2check{f}))
            cluster_pvals{f,s} = [statAll{s}.(field2check{f})(:).prob];
        else
            cluster_pvals{f,s} = [];
        end
    end
    
end

% In case you have downloaded and loaded the data, ensure stat.cfg.alpha exists:
if ~isfield(stat1.cfg,'alpha')
    stat1.cfg.alpha = 0.025; stat2.cfg.alpha = 0.025; stat3.cfg.alpha = 0.025; 
end; % stat.cfg.alpha was moved as the downloaded data was processed by an additional fieldtrip function to anonymize the data.

pos_signif_clust = {};
pos = {};
neg_signif_clust = {};
neg = {};
for s = 1:length(statAll)
    for f = 1:length(field2check)
        if ~isempty(cluster_pvals{f,s}) && strcmp(field2check{f}, 'posclusters')
            pos_signif_clust{s} = find(cluster_pvals{f,s} < stat1.cfg.alpha);
            pos{s} = ismember(statAll{s}.posclusterslabelmat, pos_signif_clust{s});
        elseif ~isempty(cluster_pvals{f,s}) && strcmp(field2check{f}, 'negclusters')
            neg_signif_clust{s} = find(cluster_pvals{f,s} < stat1.cfg.alpha);
            neg{s} = ismember(statAll{s}.negclusterslabelmat, neg_signif_clust{s});
        else
            if strcmp(field2check{f}, 'negclusters')
                neg{s} = [];
            else
                pos{s} = [];
            end
        end
    end
end



% plot
for k = 1:20;
     subplot(4,5,k);   
     cfg = [];   
     cfg.xlim=[j(k) j(k+1)];   
     cfg.zlim = [-1.0e-13 1.0e-13];   
     pos_int = all(pos(:, m(k):m(k+1)), 2);
     cfg.highlight = 'on';
     cfg.highlightchannel = find(pos_int);       
     cfg.comment = 'xlim';   
     cfg.commentpos = 'title';   
%      cfg.layout = 'CTF151.lay';
     ft_topoplotER(cfg, GA_FICvsFC);
end  




%% borrowed from circflank

addpath([RUN.dataPath, 'fieldtrip-20151012/']);
rmpath([RUN.dataPath 'eeglab13_4_4b/']);

load([dataPath, 'data/', 'RUN.mat'])
load([dataPath, 'data/', 'timelock.mat']);

cfg = [];
cfg.channel   = 'all';
cfg.latency   = 'all';
cfg.parameter = 'avg';
GA_SWSC        = ft_timelockgrandaverage(cfg,timelock.SWSCcorrect{:});
GA_DWSC        = ft_timelockgrandaverage(cfg,timelock.DWSCcorrect{:});
GA_DWDC        = ft_timelockgrandaverage(cfg,timelock.DWDCcorrect{:});
GA_DWDSC       = ft_timelockgrandaverage(cfg,timelock.DWDSCcorrect{:});

cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
GA_DWSCvsSWSC = ft_math(cfg,GA_DWSC,GA_SWSC);
GA_DWDCvsSWSC = ft_math(cfg,GA_DWDC,GA_SWSC);
GA_DWDSCvsSWSC = ft_math(cfg,GA_DWDSC,GA_SWSC);
GA_DWDCvsDWSC = ft_math(cfg,GA_DWDC,GA_DWSC);
GA_DWDsCvsDWSC = ft_math(cfg,GA_DWDSC,GA_DWSC);

clusterContrasts.DWSCvsSWSC = GA_DWSCvsSWSC;
clusterContrasts.DWDCvsSWSC = GA_DWDCvsSWSC;
clusterContrasts.DWDSCvsSWSC = GA_DWDSCvsSWSC;
clusterContrasts.DWDCvsDWSC = GA_DWDCvsDWSC;
clusterContrasts.DWDsCvsDWSC = GA_DWDsCvsDWSC;

fn = fieldnames(clusterContrasts);
stat = [];

for iContrast = 1:length(fn)
    
    input = clusterContrasts.(fn{iContrast});
    %input.individual(:,3,:) = mean(input.individual(:,ismember(input.label,channels.occipitalLeft),:),2)
    
    cfg=[];
    cfg.method           = 'triangulation';
    cfg.neighbours       = ft_prepare_neighbours(cfg,RUN.elecLoc);
    cfg.latency          = [0.2 0.3] ; % let's look at the Ninc first
    cfg.channel          = 'all';
    
    cfg.method           = 'montecarlo'; % 'montecarlo' uses pseudorandom sampling based on probabilities/properties to make estimatations
    % the following cfg settings are dependent on the method selected above
    cfg.statistic        = 'indepsamplesT';
    cfg.correctm         = 'cluster';
    cfg.clusterthreshold = 'parametric';
    cfg.clusteralpha     = 0.1;
    cfg.clusterstatistic = 'maxsum';
    cfg.minnbchan        = 3;
    cfg.tail             = 0;
    cfg.clustertail      = 0;
    cfg.alpha            = 0.1;
    cfg.numrandomization = 500;
    cfg.correcttail      = 'prob';
    cfg.avgovertime      = 'no'; % originally 'yes' in circFlank
    cfg.multivariate     = 'no';
    
    % make dummy var filled with zeros for use in statistical tests
    Nsub = 11;
    cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
    cfg.design(2,1:2*Nsub) = [1:Nsub 1:Nsub];
    cfg.ivar                = 1;
%     cfg.uvar = 2;
    dummy = input;
    dummy.avg = zeros(size(input.avg));
    
    
    % actually running stats tests
    [stat.(fn{iContrast})] = ft_timelockstatistics(cfg,input,dummy);
    
end




%% clusters like Clara --> a current work in progress
% not to be confused with moves like Jagger



cfg=[];
cfg.layout    = 'EGI128.lay';
%cfg.neighbourdist=d;
cfg.method='triangulation';
cfg.feedback='no';
neighbours = ft_prepare_neighbours(cfg)




% create neighbours
cfg=[];
cfg.method           = 'triangulation';
neighbours       = ft_prepare_neighbours(cfg,RUN.elecLoc);



alpha = 0.05;
contrast = 'N1';
latency=[0.2 0.3];
time = -200:2:798;
time_range=[find(time==latency(1)*1000):find(time==latency(2)*1000)];
timestep = 0.01;
       
    %N170 any
    
cfg = [];
cfg.channel = {'all'};
cfg.latency = latency;%[.33 .38]%[.170 .220];%[0.58 0.61]; %[.5 1.0];
cfg.neighbours=neighbours;
cfg.method = 'montecarlo';
cfg.statistic = 'depsamplesT';
cfg.correctm = 'cluster';
cfg.clusteralpha = alpha/2;
cfg.clusterstatistic = 'maxsum';
cfg.layout    = 'EGI128.lay';
cfg.minnbchan = 2;
cfg.tail = -1;
cfg.clustertail = -1;
cfg.alpha = alpha;
cfg.numrandomization = 1000;

subj = length(RUN.subNums);
design = zeros(2, 2*subj);
for i=1:subj
    design(1,i) = i;
end
for i=1:subj
    design(1, subj+i) = i;
end


design(2, 1:subj) = 1;
design(2, subj+1:2*subj) = 2;
cfg.design = design;
cfg.uvar = 1;
cfg.ivar = 2;
cfg.channel = RUN.channels.anterior

stat = ft_timelockstatistics(cfg, contrasts.DWDCminusDWSC, contrasts.DWDSCminusDWSC);
% stat = ft_timelockstatistics(cfg, grand_DWSC, grand_DWDC);
subject = [776 794 801 859 880 882 893 895 897 898 900]; 
channels2Use = [];
for iChan = 1:length(RUN.channels.anteriorLeft)
    newChan = RUN.channels.anteriorLeft{iChan}(2:end);
    channels2Use = [channels2Use str2num(newChan)];
end
cm_plotClusterTopoERPBarVWFA('pos', 1, subject, 200:2:300,...
    ga.SWSC, ga.DWSC, ga.DWDC, ga.DWDSC, stat2, 'vwfa_0.1-30Hz', channels2Use)
cm_plotClusterTopoERPBarVWFA('pos', 1, subject, 200:2:300,...
    ga.SWSC, ga.DWSC, ga.DWDC, ga.DWDSC, stat3, 'vwfa_0.1-30Hz', channels2Use)
cm_plotClusterTopoERPBarVWFA('pos', 2, subject, 200:2:300,...
    ga.SWSC, ga.DWSC, ga.DWDC, ga.DWDSC, stat3, 'vwfa_0.1-30Hz', channels2Use)
cm_plotClusterTopoERPBarVWFA('pos', 1, subject, 200:2:300,...
    ga.SWSC, ga.DWSC, ga.DWDC, ga.DWDSC, stat4, 'vwfa_0.1-30Hz', channels2Use)

% cm_plotClusterTopo_50ms('pos', 1, subject, 200:2:300,...
%     ga.SWSC, ga.DWSC, ga.DWDC, ga.DWDSC, stat2, 'vwfa_0.1-30Hz', channels2Use)

cm_plotClusterTopo_50ms_2msSteps('pos', 1, subject, 200:2:300,...
    ga.SWSC, ga.DWSC, ga.DWDC, ga.DWDSC, stat2, 'vwfa_0.1-30Hz', channels2Use)





%P2
plotClusterTopoERPBarVWFA('pos', 1, subject, 451:651, grand_SWSC, grand_DWSC, grand_DWDC, grand_DWDSC, stat, 'vwfa_2-30Hz')

plotClusterTopo_50ms('pos', 1, subject, 451:651, grand_SWSC, grand_DWSC, grand_DWDC, grand_DWDSC, stat, 'vwfa_2-30Hz')

plotClusterTopo_50ms_2msSteps('pos', 1, subject, 451:651, grand_SWSC, grand_DWSC, grand_DWDC, grand_DWDSC, stat, 'vwfa_2-30Hz')


plotClusterTopoERPBarVWFA('neg', 1, subject, 426:451, grand_SWSC, grand_DWSC, grand_DWDC, grand_DWDSC, stat, 'vwfa')


plotClusterTopoERPBarVWFA('pos', 2, subject, 451:651, grand_SWSC, grand_DWSC, grand_DWDC, grand_DWDSC, stat, '9_subj_vwfa')

plotClusterTopoERPBarVWFA('pos', 1, subject, 451:651, grand_SWSC, grand_DWSC, grand_DWDC, grand_DWDSC, stat, '11_subj_vwfa_12122014')
plotClusterTopoERPBarVWFA('neg', 1, subject, 451:651, grand_SWSC, grand_DWSC, grand_DWDC, grand_DWDSC, stat, '11_subj_vwfa_12122014')


plotClusterTopo_50ms_2msSteps('pos', 1, subject, 426:451, grand_SWSC, grand_DWSC, grand_DWDC, grand_DWDSC, stat, 'vwfa_2-30Hz')

plotClusterTopoERPBarVWFA_select('pos', 1, subject, 426:451, grand_SWSC, grand_DWSC, grand_DWDC, grand_DWDSC, stat, 'vwfa_2-30Hz')


load('/Users/clara/Desktop/2012_01_categorization/999_RAPaper/02_N1cluster/n170_50ms_p0.01_minnbchan2.mat')
plotClusterTopoERPBarVWFA_oldLeftN1('neg', 2, subject, 526:551, grand_SWSC, grand_DWSC, grand_DWDC, grand_DWDSC, stat, 'vwfa')

plotClusterTopoERPBarVWFA('neg', 1, subject, 451:551, grand_SWSC, grand_DWSC, grand_DWDC, grand_DWDSC, stat, 'vwfa')

plotClusterTopoERPBarVWFA('neg', 1, subject, 451:651, grand_SWSC, grand_DWSC, grand_DWDC, grand_DWDSC, stat, '11_subj_vwfa_12122014')
plotClusterTopoERPBarVWFA('neg', 2, subject, 451:651, grand_SWSC, grand_DWSC, grand_DWDC, grand_DWDSC, stat, '11_subj_vwfa_12122014')

plotClusterTopoERPBarVWFA('pos', 1, subject, 451:651, grand_SWSC, grand_DWSC, grand_DWDC, grand_DWDSC, stat, '11_subj_vwfa_12122014')
plotClusterTopoERPBarVWFA('pos', 2, subject, 451:651, grand_SWSC, grand_DWSC, grand_DWDC, grand_DWDSC, stat, '11_subj_vwfa_12122014')
plotClusterTopoERPBarVWFA('pos', 3, subject, 451:651, grand_SWSC, grand_DWSC, grand_DWDC, grand_DWDSC, stat, '11_subj_vwfa_12122014')
plotClusterTopoERPBarVWFA('pos', 4, subject, 451:651, grand_SWSC, grand_DWSC, grand_DWDC, grand_DWDSC, stat, '11_subj_vwfa_12122014')
plotClusterTopoERPBarVWFA('pos', 5, subject, 451:651, grand_SWSC, grand_DWSC, grand_DWDC, grand_DWDSC, stat, '11_subj_vwfa_12122014')
plotClusterTopoERPBarVWFA('pos', 6, subject, 451:651, grand_SWSC, grand_DWSC, grand_DWDC, grand_DWDSC, stat, '11_subj_vwfa_12122014')


plotClusterTopoERPBarVWFA('pos', 1, subject, 526:551, grand_SWSC, grand_DWSC, grand_DWDC, grand_DWDSC, stat, '11_subj_2_30Hz_vwfa_12292014')
%plotClusterTopoERPBarVWFA('pos', 2, subject, 526:551, grand_SWSC, grand_DWSC, grand_DWDC, grand_DWDSC, stat, '11_subj_2_30Hz_vwfa_12292014')


plotClusterTopo('pos', 1, subject, 526:551, grand_SWSC, grand_DWSC, grand_DWDC, grand_DWDSC, stat, '11_subj_2_30Hz_vwfa_12292014')

plotClusterTopoERPBarVWFA('neg', 1, subject, 526:551, grand_SWSC, grand_DWSC, grand_DWDC, grand_DWDSC, stat, '11_subj_0.1_30Hz_vwfa_correctOnly_20150113')


plotClusterTopoERPBarVWFA('neg', 1, subject, 451:651, grand_SWSC, grand_DWSC, grand_DWDC, grand_DWDSC, stat, '11_subj_0.1_30Hz_vwfa_correctOnly_20150113')


plotClusterTopo('neg', 1, subject, 451:651, grand_SWSC, grand_DWSC, grand_DWDC, grand_DWDSC, stat, '11_subj_0.1_30Hz_vwfa_12292014')

plotClusterERPBarVWFA('neg', 1, subject, 451:651, grand_SWSC, grand_DWSC, grand_DWDC, grand_DWDSC, stat, '11_subj_0.1_30Hz_vwfa_12292014')













