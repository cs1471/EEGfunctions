function [ my_preproc ] = cm_vwfaICA( EEG,currentSub,iSub )
%Simplifies preprocessing by allowing you to change parameters of RUN
%   RUN.preproc contains several ways of preprocessing
%   Output can be saved then loaded for further analysis

global RUN

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
%         rmpath([RUN.dataPath '/eeglab13_4_4b/functions/octavefunc/signal']);
        EEG = pop_resample(EEG,RUN.forICA.newSrate);
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

EEG = pop_epoch( EEG, {'DIN2'}, RUN.forICA.epochLength);

%% Remove baseline

baselineStart = RUN.forICA.baseline(1) * 1000;
baselineEnd = RUN.forICA.baseline(2) * 1000;

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

if length(RUN.forICA.interp) < iSub
    pop_eegplot(EEG) % <-- changed to show fewer channels at a time
    keyboard % look for bad channels in the scroll data, save them in 'bad'
    RUN.forICA.interp{iSub} = bad;
end

EEG = pop_interp(EEG,RUN.forICA.interp{iSub},'spherical');


%% Remove noisy epochs

% does not detect artifacts in the eye channels
[EEG Indexes] = pop_eegthresh(EEG, 1, setdiff(1:128,31:32), RUN.forICA.noiseThresh(1), RUN.forICA.noiseThresh(2),...
    RUN.forICA.epochLength(1), RUN.forICA.epochLength(2), 1, 1);

RUN.badEpochs{iSub} = Indexes';

%% Remove bad channels
% so that each channel has unique info for ICA

EEG = pop_select(EEG, 'nochannel', RUN.forICA.interp{iSub});

%% Run ICA (and save icaweights, icasphere, icawinv)
EEG = pop_runica(EEG, 'icatype', 'runica', 'extended', 1);

filename = strcat(currentSub, '_preprocICA.set');
filepath = strcat(RUN.dataPath, 'data/', currentSub, '/');
EEG = pop_saveset( EEG, 'filename', char(filename), 'filepath', char(filepath));

% used to name the saved version of my function's output
ICAstruct = [];
ICAstruct.icaweights = EEG.icaweights;
ICAstruct.icasphere = EEG.icasphere;
ICAstruct.icawinv = EEG.icawinv;
ICAstruct.icachansind = EEG.icachansind;

save([RUN.dataPath, 'data/', currentSub, '/', 'ICAstruct.mat'], 'ICAstruct');


end

