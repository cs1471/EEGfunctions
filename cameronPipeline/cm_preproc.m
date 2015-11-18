function  cm_preproc( EEG,currentSub,iSub )
%Simple preprocessing in EEGLAB using parameters in RUN (created Oct '15 by CCM in MAXLAB)
%   RUN.preproc contains several ways of preprocessing
%   Output can be saved then loaded for further analysis

%% Flip channels if needed
global RUN

if find(ismember(RUN.preproc.subjectsToFlip,currentSub)) > 0
    EEGheadbox1 = EEG.data(33:64,:);
    EEGheadbox2 = EEG.data(1:32,:);
    
    EEG.data = vertcat(EEGheadbox1, EEGheadbox2);
    
    oldInterp = RUN.preproc.interp{iSub};
    switchInterp = RUN.preproc.interp{iSub};
    
    % switch bad channels
    for h = 1:length(oldInterp)
        if oldInterp(h) >= 33
            switchInterp(h) = oldInterp(h) - 32;
        else
            switchInterp(h) = oldInterp(h) + 32;
        end
    end
    
    RUN.preproc.interp{iSub} = switchInterp;
        
end



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

EEG = pop_epoch( EEG, RUN.preproc.stimuli, RUN.preproc.epochLength);

%% Remove baseline

epochStart = RUN.preproc.epochLength(1) * 1000;

EEG = pop_rmbase( EEG, [epochStart    0]);


%% Load layout

EEG = pop_chanedit(EEG, 'load',{[RUN.dataPath, 'analysis\\mw64_with2VEOGs.ced'] 'filetype' 'autodetect'});

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

% Reref in EEGlab to average mastoid.
% first step is to add a ref channel to the channel location structure (its
% basically an empty "dummy" channel). channel locations are not needed at
% this point
EEG = pop_chanedit(EEG, 'lookup',fullfile(fileparts(which('eeglab.m')),'/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385NEW.elp'),...
    'append',64,'changefield',{65 'labels' 'RM'},'setref',{'1:64' 'RM'});

% rereference to average over all channels (in pop_reref the second input; [] does
% this), 'refloc' with label Ch33 retains the old reference, it basically
% fills in the dummy channel as the old reference channel....
EEG = pop_reref( EEG, [],'refloc',struct('labels',{'RM'},'theta',{[]},'radius',{[]},'X',{[]},'Y',{[]},...
    'Z',{[]},'sph_theta',{[]},'sph_phi',{[]},'sph_radius',{[]},'type',{[]},'ref',{[]},'urchan',{[65]}));

% rereference to average mastoid NOT KEEPING the reference channels
% index channel 32 is the Right mastoid (channel 33 on the cap) 64 is the
% dummy channel we added in earlier and becomes the left mastoid...
%EEG = pop_reref( EEG, [32 64] );
% rereference to average mastoid KEEPING the reference channels

switch RUN.preproc.reref
    case 'average'
        EEG = pop_reref( EEG, 1:65 ,'keepref','on');
        
    case 'average_mastoid'
        EEG = pop_reref( EEG, [33 65] ,'keepref','on');
end

EEG = pop_chanedit(EEG, 'load',{[RUN.dataPath, 'analysis\\mw64_withMastoids_and2VEOGs.ced'] 'filetype' 'autodetect'});



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
        [EEG Indexes] = pop_eegthresh(EEG, 1, [1:65], RUN.preproc.artCriteria(1), RUN.preproc.artCriteria(2), RUN.preproc.epochLength(1), RUN.preproc.epochLength(2), 0,0);
        
        RUN.arfEpochs{iSub} = Indexes';
        RUN.numTrials{iSub} = size(EEG.data,3);
end


%% Save data

EEG = pop_saveset( EEG, 'filename', [currentSub, '_preproc.set'], 'filepath', [RUN.dataPath, 'data\', currentSub, '\']);


% next steps will be creating trialStruct and conditions

end

