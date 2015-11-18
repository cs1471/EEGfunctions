function  cm_vwfaPreproc( EEG,currentSub,iSub )
%Simple preprocessing in EEGLAB using parameters in RUN (created Oct '15 by CCM in MAXLAB)
%   RUN.preproc contains several ways of preprocessing
%   Output can be saved then loaded for further analysis

global RUN

%% Filter - currently in EEGLAB
%EEG = pop_eegfiltnew(EEG, [], 40, 166, 0, [], 1);
% [EEG, com, b] = pop_eegfiltnew(EEG, locutoff, hicutoff, filtorder,
%                                       revfilt, usefft, plotfreqz, minphase);
% EEGOUT = pop_iirfilt( EEG, locutoff, hicutoff, trans_bw, revfilt, causal);

highpass = RUN.preproc.filterSettings{1,1};
lowpass = RUN.preproc.filterSettings{1,2};

switch RUN.preproc.filterType
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

switch RUN.preproc.resample
    case 'yes'
        rmpath([RUN.dataPath '\matlab_applications\eeglab13_3_2b\functions\octavefunc\signal']);
        %C:\Users\ccm25\circFlank\matlab_applications\eeglab12_0_2_5b\functions\octavefunc\signal

        EEG = pop_resample(EEG,RUN.preproc.newSrate);
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

EEG = pop_epoch( EEG, {'DIN2'}, RUN.preproc.epochLength);

%% Remove baseline

baselineStart = RUN.preproc.baseline(1) * 1000;
baselineEnd = RUN.preproc.baseline(2) * 1000;

EEG = pop_rmbase( EEG, [baselineStart    baselineEnd]);


%% Load layout

EEG = pop_chanedit(EEG, 'load',{[RUN.dataPath, 'Hydrocel_GSN_128_1.0_TRIM.sfp'] 'filetype' 'autodetect'});

%% Remove bad channels
% copy EEG.chanlocs so we can interp bad channels later
chanlocs = EEG.chanlocs;

% actually removing bad channels
EEG = pop_select(EEG, 'nochannel', RUN.preproc.interp{iSub});

%% Reject components

% copy in ICA stuff (same number of channels/components now)
EEG.icaweights = RUN.icaweights{1,iSub};
EEG.icasphere = RUN.icasphere{1,iSub};
EEG.icawinv = RUN.icawinv{1,iSub};
EEG.icachansind = RUN.icachansind{1,iSub};

switch RUN.preproc.rejComps
    case 'yes_first' % allows for debugging on the first run; use when ICA has changed
        
        pop_topoplot(EEG,0)
        
        % breakpoint for debugging
        keyboard
        
        % make a vector called components
        EEG = pop_subcomp(EEG, components, 1);
        RUN.ICAcomps{iSub} = components;
    
    case 'yes' % no debugging
        
        EEG = pop_subcomp(EEG, RUN.ICAcomps{iSub}, 0);

        
end

%% Interpolate bad channels

if RUN.preproc.interp{iSub} ~= 0
    EEG = pop_interp(EEG, chanlocs,'spherical');
end

%% Re-reference
% nope


%% Detect bad epochs
% also make a list of how many trials are detected

% DON'T ACTUALLY REMOVE ANY EPOCHS!!

switch RUN.preproc.artDetection
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
        RUN.arfEpochs{iSub} = Indexes';
        RUN.numTrials{iSub} = size(EEG.data,3);
end


%% Save data
EEG = pop_saveset( EEG, 'filename', [currentSub, '_preproc.set'], 'filepath', [RUN.dataPath, 'data\', currentSub, '\']);


% next steps will be creating trialStruct and conditions

end

