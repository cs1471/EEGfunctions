function [ my_preproc ] = cm_maxlabPreproc( EEG,currentSub,iSub,RUN )
%Simplifies preprocessing by allowing you to change parameters of RUN
%   RUN.maxlabPreproc contains several ways of preprocessing
%   Output can be saved then loaded for further analysis

global RUN

%% Filter - currently in EEGLAB
%EEG = pop_eegfiltnew(EEG, [], 40, 166, 0, [], 1);
% [EEG, com, b] = pop_eegfiltnew(EEG, locutoff, hicutoff, filtorder,
%                                       revfilt, usefft, plotfreqz, minphase);
% EEGOUT = pop_iirfilt( EEG, locutoff, hicutoff, trans_bw, revfilt, causal);

highpass = RUN.maxlabPreproc.filterSettings{1,1};
lowpass = RUN.maxlabPreproc.filterSettings{1,2};

switch RUN.maxlabPreproc.filterType
    case 'causal'
%         EEG = pop_iirfilt( EEG, highpass, lowpass, [], 0, 1);
%         EEG = pop_eegfiltnew(EEG, highpass, lowpass, [],...
%             0, [], 0, true);

        EEG = pop_eegfiltnew( EEG, highpass, lowpass, [], [0],'minphase',1);

    case 'non-causal'
        EEG = pop_eegfiltnew(EEG, highpass, lowpass, 3300,...
            0, [], 0, false);
end

%% Re-sample!!
% This is done after filtering to avoid aliasing

switch RUN.maxlabPreproc.resample
    case 'yes'
%         rmpath([RUN.dataPath '/eeglab13_4_4b/functions/octavefunc/signal']);
        EEG = pop_resample(EEG,RUN.maxlabPreproc.newSrate);
end


%% Epoch
%Make the DINS even/odd
for ii=1:length(EEG.event)
    order = ii;%EEG.event(i).urevent;
    if mod(order,2) %mod(i,2)
        EEG.event(ii).type = 'DIN1';
    else
        EEG.event(ii).type = 'DIN2';
    end
end

EEG = pop_epoch( EEG, {'DIN2'}, RUN.maxlabPreproc.epochLength);

%% Remove baseline

baselineStart = RUN.maxlabPreproc.baseline(1) * 1000;
baselineEnd = RUN.maxlabPreproc.baseline(2) * 1000;

EEG = pop_rmbase( EEG, [baselineStart    baselineEnd]);


%% Load layout
% useful for interpolating bad channels
EEG = pop_chanedit(EEG, 'load',{[RUN.dataPath, 'Hydrocel_GSN_128_1.0_TRIM.sfp'] 'filetype' 'autodetect'});

%% Interpolate bad channels
% then reject bad channels after removing noisy epochs

% this makes sure there are consistantly 128 total channels in each epoch,
% but each channel in ICA has unique data (not averaged from other channels)

% maybe do something here to actually visualize all channels, then make a
% vector with the indexes of each bad channel...

if length(RUN.maxlabPreproc.interp) < iSub
    pop_eegplot(EEG) % <-- changed to show fewer channels at a time
    keyboard % look for bad channels in the scroll data, save them in 'bad'
    RUN.maxlabPreproc.interp{iSub} = bad;
end

EEG = pop_interp(EEG,RUN.maxlabPreproc.interp{iSub},'spherical');


%% Detect bad epochs
% also make a list of how many trials are detected

% DON'T ACTUALLY REMOVE ANY EPOCHS!!

switch RUN.maxlabPreproc.artDetection
    case 'linear'
%         OUTEEG = pop_rejtrend( INEEG, typerej, elec_comp, ...
%                  winsize, maxslope, minR, superpose, reject,calldisp);
        winsize = EEG.srate * (RUN.preproc.epochLength(2) - RUN.preproc.epochLength(1));
        
        EEG = pop_rejtrend(EEG, 1, [1:65],winsize, RUN.preproc.artCriteria(1), RUN.preproc.artCriteria(2), 1,0,0);
        
        RUN.arfEpochs{iSub} = [];
        RUN.numTrials{iSub} = size(EEG.data,3);
        
    case 'thresh'
%         [EEG Indexes] = pop_eegthresh( INEEG, typerej, elec_comp, lowthresh, ...
%             upthresh, starttime, endtime, superpose, reject);

% reject using all channels and the entire length of the epoch
%         [EEG Indexes] = pop_eegthresh(EEG, 1, [1:128], RUN.preproc.artCriteria(1), RUN.preproc.artCriteria(2), RUN.preproc.epochLength(1), RUN.preproc.epochLength(2), 0,0);

% reject using all channels and only part of the epoch
%         [EEG Indexes] = pop_eegthresh(EEG, 1, [1:128],-90, 90,-.2, 0.4,1,0);
        
% reject using eye channels and only part of the epoch
        [EEG Indexes] = pop_eegthresh(EEG, 1, [14 21 126 127],-80, 80,-.2, 0.4,1,0);
        RUN.maxlabArfEpochs{iSub} = Indexes';
        RUN.maxlabNumTrials{iSub} = size(EEG.data,3);
end


%% Save data
EEG = pop_saveset( EEG, 'filename', [currentSub, '_maxlabPreprocCausal.set'], 'filepath', [RUN.dataPath, 'data\', currentSub, '\']);



end

