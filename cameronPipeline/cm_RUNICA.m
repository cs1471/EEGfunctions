function [ my_preproc ] = cm_RUNICA( EEG,currentSub,iSub,RUN )
%Simplifies preprocessing by allowing you to change parameters of RUN
%   RUN.preproc contains several ways of preprocessing
%   Output can be saved then loaded for further analysis


%% Filter - currently in EEGLAB
%EEG = pop_eegfiltnew(EEG, [], 40, 166, 0, [], 1);
% [EEG, com, b] = pop_eegfiltnew(EEG, locutoff, hicutoff, filtorder,
%                                       revfilt, usefft, plotfreqz, minphase);
% EEGOUT = pop_iirfilt( EEG, locutoff, hicutoff, trans_bw, revfilt, causal);

highpass = RUN.forICA.filterSettings{1,1};
lowpass = RUN.forICA.filterSettings{1,2};

switch RUN.forICA.filterType
    case 'causal'
%         EEG = pop_iirfilt( EEG, highpass, lowpass, [], 0, 1);
        EEG = pop_eegfiltnew(EEG, highpass, lowpass, [],...
            0, [], 0, true);

    case 'non-causal'
        EEG = pop_eegfiltnew(EEG, highpass, lowpass, 3300,...
            0, [], 0, false);
end

%% Re-sample!!
% This is done after filtering to avoid aliasing

switch RUN.forICA.resample
    case 'yes'
        rmpath([RUN.dataPath '\matlab_applications\eeglab13_3_2b\functions\octavefunc\signal']);
        %C:\Users\ccm25\circFlank\matlab_applications\eeglab12_0_2_5b\functions\octavefunc\signal

        EEG = pop_resample(EEG,RUN.forICA.newSrate);
end


%% Epoch

EEG = pop_epoch( EEG, RUN.forICA.stimuli, RUN.forICA.epochLength);

%% Remove baseline

epochStart = RUN.forICA.epochLength(1) * 1000;

EEG = pop_rmbase( EEG, [epochStart    0]);

%% Re-reference
% not done until after ICA

%% Load layout
% useful for interpolating bad channels
% layout has 2 VEOGs and only 1 mastoid (LM)
EEG = pop_chanedit(EEG, 'load',{[RUN.dataPath, 'analysis\mw64_with2VEOGs.ced'] 'filetype' 'autodetect'});

%% Interpolate bad channels
% then reject bad channels after removing noisy epochs

% this makes sure there are consistantly 64 total channels in each epoch,
% but each channel in ICA has unique data (not averaged from other channels)

if RUN.forICA.interp{iSub} ~= 0
    EEG = pop_interp(EEG,RUN.forICA.interp{iSub},'spherical');
end


%% Remove noisy epochs

% does not detect artifacts in the eye channels
[EEG Indexes] = pop_eegthresh(EEG, 1, setdiff(1:64,31:32), RUN.forICA.noiseThresh(1), RUN.forICA.noiseThresh(2),...
    RUN.forICA.epochLength(1), RUN.forICA.epochLength(2), 1, 1);

RUN.badEpochs{iSub} = Indexes';

%% Remove bad channels
% so that each channel has unique info for ICA

EEG = pop_select(EEG, 'nochannel', RUN.forICA.interp{iSub});

%% Run ICA (and save icaweights, icasphere, icawinv)
EEG = pop_runica(EEG, 'icatype', 'runica', 'extended', 1);
EEG = pop_saveset( EEG, 'filename', [currentSub, '_preprocICA.set'], 'filepath', [RUN.dataPath, 'data\', currentSub, '\']);

% used to name the saved version of my function's output
ICAstruct = [];
ICAstruct.icaweights = EEG.icaweights;
ICAstruct.icasphere = EEG.icasphere;
ICAstruct.icawinv = EEG.icawinv;

save([RUN.dataPath, 'data\', currentSub, '\', 'ICAstruct.mat'], 'ICAstruct');


end

