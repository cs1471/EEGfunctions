function ga = cm_fixChanLocs(ga, fn, iConditions)
%cm_fixChanLocs fixes the location of electrodes when ga is first created
%   This makes that anterior is at the top of the figure when doing
%   ft_multiplotER later

chanposX = ga.(fn{iConditions}).elec.chanpos(:,1);
chanposY = ga.(fn{iConditions}).elec.chanpos(:,2);
chanposZ = ga.(fn{iConditions}).elec.chanpos(:,3);
actualChannelPos = [-chanposY, chanposX, chanposZ];

ga.(fn{iConditions}).elec.chanpos = actualChannelPos;
ga.(fn{iConditions}).elec.elecpos = ga.(fn{iConditions}).elec.chanpos;


end

