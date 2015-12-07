function cm_plotClusterTopoERPBarVWFA(sign, number, subject, time_range, gaSWSC, gaDWSC, gaDWDC, gaDWDSC,...
    stat, contrast, channels, ga)

%time_range=time_range+200;
qq=number;

% grandm0vsm6 = gaDWDSC;
% grandm0vsm6.individual = gaDWDSC.individual - gaSWSC.individual;

% topoContrast = gaDWDSC;
cfg = [];
cfg.operation = 'subtract';
% cfg.parameter = 'individual';
% topoContrast.individual = ft_math(cfg, ga.(stat.cfg.timelock2compare{1}), ga.(stat.cfg.timelock2compare{2}));
cfg.parameter = 'avg';
topoContrast = ft_math(cfg, ga.(stat.cfg.timelock2compare{1}), ga.(stat.cfg.timelock2compare{2}));

scrsz=[0 0 2560 1440];
figure('Position', [scrsz(1), scrsz(2), scrsz(3)*.75, scrsz(4)*.5]);


if strcmp(sign, 'pos')
    titleText=(['cluster p = ' num2str(stat.posclusters(qq).prob, '%1.3f')]);
    value = ismember(stat.posclusterslabelmat, qq);
    duration = ismember(stat.posclusterslabelmat, qq);
elseif strcmp(sign, 'neg')
    titleText=(['cluster p = ' num2str(stat.negclusters(qq).prob, '%1.3f')]);
    value = ismember(stat.negclusterslabelmat, qq);
    duration = ismember(stat.negclusterslabelmat, qq);
end


% channels=[1:128];

clusterChan = find(any(value(:,:), 2));

durationVector=stat.time(any(duration,1));

cfg = [];
cfg.xlim = [min(durationVector) max(durationVector)];
cfg.highlight = 'on';

cfg.highlightchannel = clusterChan;
cfg.comment = 'xlim';
cfg.commentpos = 'title';
%     cfg.layout    = 'EGI128.lay';
cfg.zlim=[-2 2];
cfg.highlightsize = 6;
cfg.highlightcolor='w';
cfg.colorbar ='yes';
% cfg.colormap = jet;
cfg.parameter = 'avg';

%% plot 1 - contrast topo with cluster channels highlighted
subplot(131)
% subplot(1,3,1)

ft_topoplotER(cfg, topoContrast);
fixFonts_topo
freezeColors
cbfreeze

%     title([titleText ' time [' num2str(min(durationVector)) ' ' num2str(max(durationVector)) ']';...
%         stat.cfg.compare])

if sum(size(stat.cfg.compare)) > 16
    % make sure the title isn't too long
    newTitle = strsplit(stat.cfg.compare, 'vs');
    textLine1 = '';
    textLine2 = '';
    
    if length(newTitle{1}) > length(newTitle{2})
        textLine1 = newTitle{1};
        textLine2 = strcat('vs ', newTitle{2});
    else
        textLine1 = strcat(newTitle{1}, ' vs');
        textLine2 = newTitle{2};
    end
    
    title({[titleText]; ['time: ' num2str(min(durationVector)*1000) '-' num2str(max(durationVector)*1000) ' ms'];...
        [textLine1]; [textLine2]})
else
    title({[titleText]; ['time: ' num2str(min(durationVector)*1000) '-' num2str(max(durationVector)*1000) ' ms'];...
        [stat.cfg.compare]})
end

fixFonts

%% plot 2 - ERPs
% subplot(132)
subplot(1,3,2)

%secondStim=1:400;
%secondStim=201:600;
secondStim=351:750

range=[-200:2:598];
% xlim([-200 600])



        l1 = plot(range,squeeze(mean(gaSWSC.avg(clusterChan,secondStim))), 'b', 'LineWidth', 6);
        hold on;
        l2 = plot(range,squeeze(mean(gaDWSC.avg(clusterChan,secondStim))), 'c--', 'LineWidth', 5);
        l3 = plot(range,squeeze(mean(gaDWDC.avg(clusterChan,secondStim))), 'm--', 'LineWidth', 5);
        l4 = plot(range,squeeze(mean(gaDWDSC.avg(clusterChan,secondStim))), 'r', 'LineWidth', 6);
     

% ERP2plot1 = ga.(stat.cfg.timelock2compare{1});
% ERP2plot2 = ga.(stat.cfg.timelock2compare{2});
% 
% ERPcolor1 = cm_assignColor(ga, stat.cfg.timelock2compare{1});
% ERPcolor2 = cm_assignColor(ga, stat.cfg.timelock2compare{2});
% 
% plot(range,squeeze(mean(ERP2plot1.avg(clusterChan,secondStim))), ERPcolor1{1}, 'LineWidth', ERPcolor1{2})
% hold on;
% plot(range,squeeze(mean(ERP2plot2.avg(clusterChan,secondStim))), ERPcolor2{1}, 'LineWidth', ERPcolor2{2})

%title('P2 cluster ERPs')
ylabel('\muV')

% legend([11 12 l3 l4],{'SWSC', 'DWSC', 'DWDC', 'DWDSC'}, 'Location', 'SouthWest');
legend('SWSC', 'DWSC', 'DWDC', 'DWDSC', 'Location', 'SouthWest')
% legend(stat.cfg.timelock2compare{1}, stat.cfg.timelock2compare{2},...
%     stat.cfg.timelock2compare{3}, stat.cfg.timelock2compare{4}, 'Location', 'SouthWest')
% legend(stat.cfg.timelock2compare{1}, stat.cfg.timelock2compare{2}, 'Location', 'SouthWest')

%ylim([-7.5 2])
xlim([-200 600])

ylims=ylim;
fixFonts(gca, 14, 3)
plot([-200 598], [0 0], 'k--', 'LineWidth', 1)
%plot([-400 -400], [ylims(1) ylims(2)], 'k--', 'LineWidth', 1)
%plot([-200 -200], [ylims(1) ylims(2)], 'k--', 'LineWidth', 1)
plot([0 0], [ylims(1) ylims(2)], 'k--', 'LineWidth', 1)
plot([300 300], [ylims(1) ylims(2)], 'k--', 'LineWidth', 1)

% now to plot the signficant time window
if max(durationVector) == min(durationVector)
    sigLine = plot([max(durationVector)*1000 max(durationVector)*1000], [ylims(1) ylims(2)], 'y', 'LineWidth', 2);
    uistack(sigLine,'bottom')
else
%     sigLine1 = plot([min(durationVector)*1000 min(durationVector)*1000], [ylims(1) ylims(2)], 'y', 'LineWidth', 2);
%     sigLine2 = plot([max(durationVector)*1000 max(durationVector)*1000], [ylims(1) ylims(2)], 'y', 'LineWidth', 2);
%     uistack(sigLine1,'bottom'); uistack(sigLine2,'bottom');
% %     significantPatch =patch([min(durationVector)*1000 max(durationVector)*1000 max(durationVector)*1000 min(durationVector)*1000],...
% %         [ylims(1) ylims(1) ylims(2) ylims(2)],'y');
% %     set(significantPatch,'FaceAlpha',0.5);
    % highlightpatchOne(min(durationVector)*1000, max(durationVector)*1000)
    
    sigPatch = rectangle('Position',...
        [min(durationVector)*1000 ylims(1) (max(durationVector)*1000 - min(durationVector)*1000) abs(ylims(1) - ylims(2))],...
        'FaceColor',[1 1 0]);
    uistack(sigPatch,'bottom')
    
    
    
end

title({titleText; ['time: ' num2str(min(durationVector)*1000) '-' num2str(max(durationVector)*1000) ' ms']})
fixFonts%(gca, 14, 3)
freezeColors

%% plot 3 - error bars

    m0=zeros(length(subject),1);
    m3w=zeros(length(subject),1);
    m3b=zeros(length(subject),1);
    m6=zeros(length(subject),1);

% gaSWSC = ERP2plot1;
% gaDWSC = ERP2plot2;

fn = fieldnames(ga);
allPvals = zeros(length(subject),length(fn));

if strcmp(sign, 'pos')
    
    for iSub = 1:length(subject)
        for iCond = 1:length(fn)
%             if strcmp(fn{iCond}, stat.cfg.timelock2compare{1}) || strcmp(fn{iCond}, stat.cfg.timelock2compare{2})
                allPvals(iSub,iCond) = sum(sum(squeeze(ga.(fn{iCond}).individual(iSub,channels,time_range)).*(stat.posclusterslabelmat==qq)))/sum(sum(stat.posclusterslabelmat==qq));
%             end
        end
    end
    
    
    %         for jjj=1:length(subject)
    %             m0(jjj) = sum(sum(squeeze(gaSWSC.individual(jjj,channels,time_range)).*(stat.posclusterslabelmat==qq)))/sum(sum(stat.posclusterslabelmat==qq));
    %             m3w(jjj) = sum(sum(squeeze(gaDWSC.individual(jjj,channels,time_range)).*(stat.posclusterslabelmat==qq)))/sum(sum(stat.posclusterslabelmat==qq));
    %             m3b(jjj) = sum(sum(squeeze(gaDWDC.individual(jjj,channels,time_range)).*(stat.posclusterslabelmat==qq)))/sum(sum(stat.posclusterslabelmat==qq));
    %             m6(jjj) = sum(sum(squeeze(gaDWDSC.individual(jjj,channels,time_range)).*(stat.posclusterslabelmat==qq)))/sum(sum(stat.posclusterslabelmat==qq));
    %
    %         end
    
    titleText=(['cluster p=' num2str(stat.posclusters(qq).prob, '%1.3f')]);
    
elseif strcmp(sign, 'neg')
    
    for iSub = 1:length(subject)
        for iCond = 1:length(fn)
            allPvals(iSub,iCond) = sum(sum(squeeze(ga.(fn{iCond}).individual(iSub,channels,time_range)).*(stat.negclusterslabelmat==qq)))/sum(sum(stat.negclusterslabelmat==qq));
        end
    end
%         for jjj=1:length(subject)
%             m0(jjj) = sum(sum(squeeze(gaSWSC.individual(jjj,channels,time_range)).*(stat.negclusterslabelmat==qq)))/sum(sum(stat.negclusterslabelmat==qq));
%             m3w(jjj) = sum(sum(squeeze(gaDWSC.individual(jjj,channels,time_range)).*(stat.negclusterslabelmat==qq)))/sum(sum(stat.negclusterslabelmat==qq));
%             m3b(jjj) = sum(sum(squeeze(gaDWDC.individual(jjj,channels,time_range)).*(stat.negclusterslabelmat==qq)))/sum(sum(stat.negclusterslabelmat==qq));
%             m6(jjj) = sum(sum(squeeze(gaDWDSC.individual(jjj,channels,time_range)).*(stat.negclusterslabelmat==qq)))/sum(sum(stat.negclusterslabelmat==qq));
%         end

        titleText=(['cluster p=' num2str(stat.negclusters(qq).prob, '%1.3f')]);

end


% figure out which columns of allPvals to use and what color(s) to use
index1 = find(strcmp(fn, stat.cfg.timelock2compare{1}));
index2 = find(strcmp(fn, stat.cfg.timelock2compare{2}));

cond1 = allPvals(:,index1);
cond2 = allPvals(:,index2);

colors = [0 0 0; 0 0 0];
for iColor = 1:length(stat.cfg.timelock2compare)
    if strfind(stat.cfg.timelock2compare{iColor}, 'SWSC')
        colors(iColor,:) = [0 0 1]; % blue
    elseif strfind(stat.cfg.timelock2compare{iColor}, 'DWSC')
        colors(iColor,:) = [0 1 1]; % cyan
    elseif strfind(stat.cfg.timelock2compare{iColor}, 'DWDC')
        colors(iColor,:) = [1 0 1]; % magenta
    elseif strfind(stat.cfg.timelock2compare{iColor}, 'DWDSC')
        colors(iColor,:) = [1 0 0]; % red
    elseif strfind(stat.cfg.timelock2compare{iColor}, 'DWbothAnimalOrTool')
        colors(iColor,:) = [0 1 0]; % lime green
    elseif strfind(stat.cfg.timelock2compare{iColor}, 'SameAnimalOrTool')
        colors(iColor,:) = [0 0.4 0]; % dark green
    elseif strfind(stat.cfg.timelock2compare{iColor}, 'AcrossCategories')
        colors(iColor,:) = [1 0.6 0]; % organge
    else
        colors(iColor,:) = [1 1 1];
    end
    
end

% remove duplicate colors, so bar graphs aren't both the same color
if sum(size(unique(colors, 'rows'))) == 4
    colors(2,:) = [0 0 0];
end

subplot(133)
% subplot(1,3,3)
cm_writepvalue2condNEW(cond1, cond2, colors)
% writepvalue3cond(m0, m3w, m3b)
%UNCOMMENT ME!! writepvalueNeg(m0, m3w, m3b, m6)
title(titleText)
fixFonts


%   saveFormats('test', ['TopoERPBarCluster' sign num2str(number) 'Contrast' contrast 'clusterAlpha' num2str(stat.cfg.clusteralpha) 'minnbchan' num2str(stat.cfg.minnbchan) ...
%       'latency' num2str(stat.cfg.latency(1)) '_' num2str(stat.cfg.latency(2)) 'numRand' num2str(stat.cfg.numrandomization)], './', 100)

end
