function trialStruct = cm_vwfaTrialStruct( RUN, EEG, iSub );
%% cm_vwfaTrialStruct.m is a combination of createTrialStruct and convertToMaxLab_vwfa
% Basically sorting trials and adding relevant info to a dataset

% initialize trialStruct
trialStruct = dataset();

%  Convert to our own format
exptData = load(strcat(RUN.dataPath, 'data/', RUN.subNums{iSub}, '/', RUN.filenames.exptFilename{iSub}));
numSessions = size(exptData.trialOutput,2);
numTrials = size(exptData.trialOutput(1).trials,2);

% % create trial list
trialNum = zeros(1,length(EEG.epoch));
stim = zeros(1,length(EEG.epoch));
condition = zeros(1,length(EEG.epoch)); % fill this in with DWSC, etc. later
responseTime = zeros(1,length(EEG.epoch));
responsePress = zeros(1,length(EEG.epoch));
isCorrect = zeros(1,length(EEG.epoch));
isCorrectREAL = zeros(1,length(EEG.epoch));
correctResp = zeros(1,length(EEG.epoch));

totaltrials = 0;
 
% Go through each session to fill in values like RT, correct, etc.
 for session=1:numSessions
     %  Go through each trial
     numTrials = size(exptData.trialOutput(session).trials,2);
     for trial=1:numTrials
         index = trial+((session-1)*numTrials);
         %  Record the type of stimulus presented (1 = target, 0 = distractor)
         trialNum(index) = index;
         stim(index) = exptData.trialOutput(session).trials(trial).correctResponse;
         responsePress(index) = exptData.trialOutput(session).trials(trial).subjectAnimalResponse;
         responseTime(index) = exptData.trialOutput(session).trials(trial).animalResponseFinished - exptData.trialOutput(session).trials(trial).animalResponseStart;
         isCorrect(index) = exptData.trialOutput(session).trials(trial).subjectHasRightAnswer;
     end
 end

% create exptFilename
exptFilenmae = '';

if (strcmp(RUN.subNums{iSub}, '880') || strcmp(RUN.subNums{iSub}, '882') || strcmp(RUN.subNums{iSub}, '900'))
    exptFilename = ['s' RUN.subNums{iSub} '_vwfa.2408.mat'];
else
    exptFilename = ['s' RUN.subNums{iSub} '.2408.mat'];
end

% Go through the experimental design matrix to get correct responses and
% conditions
if strcmp('s895.2408.mat', exptFilename)
     sameDiffCond = exptData.exptdesign.designMatrix;
     sameDiffCond2 = exptData.exptdesign.designMatrix2;
     
     for i=1:408
         correctResp(i) = cell2mat(sameDiffCond(i,2));
         condition(i) = cell2mat(sameDiffCond(i,1));
         correctResp(i+length(sameDiffCond)) = cell2mat(sameDiffCond2(i,2));
         condition(i+length(sameDiffCond)) = cell2mat(sameDiffCond2(i,1));
     end
     
     word1 = [sameDiffCond(:,3); sameDiffCond2(:,3)];
     word2 = [sameDiffCond(:,4); sameDiffCond2(:,4)];
     
else
     sameDiffCond = exptData.exptdesign.designMatrix;
     for i=1:length(sameDiffCond)
         correctResp(i) = cell2mat(sameDiffCond(i,2));
         condition(i) = cell2mat(sameDiffCond(i,1));
         correctResp(i+length(sameDiffCond)) = cell2mat(sameDiffCond(length(sameDiffCond)+1-i,2));
         condition(i+length(sameDiffCond)) = cell2mat(sameDiffCond(length(sameDiffCond)+1-i,1));
     end
     
     word1 = [sameDiffCond(:,3); fliplr(sameDiffCond(:,3)')'];
     word2 = [sameDiffCond(:,4); fliplr(sameDiffCond(:,4)')'];
end

% actually looking for correct responses
for i=1:length(responsePress)
    if correctResp(i) == responsePress(i)
        isCorrectREAL(i) = 1;
    else
        isCorrectREAL(i) = 0;
    end
end
    

   
% making everything be a part of trialStruct
trialStruct.trialNum = trialNum';
trialStruct.stim = stim';
trialStruct.condition = condition';
trialStruct.responseTime = responseTime';
trialStruct.responsePress = responsePress';
trialStruct.correctResponse = correctResp';
trialStruct.isCorrect = isCorrectREAL';

% add the stimuli too
trialStruct.word1 = word1;
trialStruct.word2 = word2;





