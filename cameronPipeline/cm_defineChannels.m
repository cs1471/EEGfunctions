function channels = cm_defineChannels(elecLoc, label)
%cm_defineChannels divides channels based on their locations
%   After channels are localized, they are sorted and saved as channels

% make a bunch of categories
% define anterior channels and posterior channels
anterior = elecLoc(:,2) > 0;
posterior = elecLoc(:,2) <= 0;
left = elecLoc(:,1) < 0;
right = elecLoc(:,1) > 0;

% add each category to channels
channels.anterior = label(anterior);
channels.posterior = label(posterior);
channels.left = label(left);
channels.right = label(right);
channels.anteriorLeft = label(anterior & left);
channels.anteriorRight = label(anterior &  right);
channels.posteriorLeft = label(posterior & left);
channels.posteriorRight = label(posterior & right);


end

