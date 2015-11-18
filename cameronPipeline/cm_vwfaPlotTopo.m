function cm_vwfaPlotTopo( cfg, data, conds, numSteps, steps, startingTime )
%cm_plotTopoVWFA Make a figure with a bunch of topos
%   Can take very general inputs

% make sure the cfg.parameter is present
for iCond = 1:length(conds)

    if ~isfield(data.(conds{iCond}),cfg.parameter) && strcmp(cfg.parameter,'avg')
        data.(conds{iCond}).avg = squeeze(mean(data.(conds{iCond}).individual,1));
    elseif ~isfield(data.(conds{iCond}),cfg.parameter) && ~strcmp(cfg.parameter,'avg')
        disp('The selected cfg.parameter does not exist!')
        return
    end
end


% make a large, empty figure and start filling it in
figure('units', 'normalized', 'outerposition', [0 0 1 1]);

startingTime = startingTime - steps;
numConds = length(conds);

for iTime = 1:numSteps

    cfg.xlim = [startingTime + iTime*steps (startingTime+steps)+iTime*steps];
    
    startTime = num2str((startingTime+iTime*steps)*1000);
    endTime = num2str((startingTime+iTime*steps+steps)*1000);
    
    % making sure 0 is actually printed as '0'
    [dummy digits] = size(endTime);
    if digits > 5
        endTime = '0';
    end
    
    % loop through each condition
    for iCond = 1:numConds
        if iCond == 1
            subplot(numConds,numSteps,iTime);ft_topoplotER(cfg,data.(conds{iCond}))
        else
            subplot(numConds,numSteps,iTime+(numSteps*(iCond-1)));ft_topoplotER(cfg,data.(conds{iCond}))
        end
        
        if iTime == ceil(numSteps/2)
%             topTitle = [conds{iCond} ' topo'];
            theTitle = title({[conds{iCond} ' topo']; [startTime ' to ' endTime ' ms']});
            set(theTitle,'FontSize',14);
        else
            theTitle = title([startTime ' to ' endTime ' ms']);
            set(theTitle,'FontSize',14);
        end
    
    
    end
    
    
    
end


end

