function dataAvg = cm_createAvg(data)
%cm_createAvg take a fieldtrip data structure and adds an 'avg' field
%   Requires that the fieldtrip data structure contains 'individual'

fn = fieldnames(data);

for iField = 1:length(fn)
    if isfield(data.(fn{iField}),'individual')
        data.(fn{iField}).avg = squeeze(mean(data.(fn{iField}).individual,1));
    else
        disp('Your data does not contain the required .individual field!')
        return
    end
  
end

dataAvg = data;
    
end

