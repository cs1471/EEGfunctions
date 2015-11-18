function dataBaseline = cm_vwfaBaseline( cfg, data )
%cm_vwfaBaseline will baseline timelock based on values in cfg
%   Exports data that has been baselined so you can make ga and contrasts

fn = fieldnames(data);

for iParam = 1:length(cfg.fields)
    
    param = cfg.fields{iParam};
    
    % loop through each condition (each field in the data)
    for iField = 1:length(fn)
        
        cfg.parameter = param;
        
        
            %loop through each subject
            if strcmp(param,'avg') 
                for iSub = 1:length(data.(fn{iField}))
                    data.(fn{iField}){iSub} = ft_timelockbaseline(cfg, data.(fn{iField}){iSub});
                end
            end
            
      
        
    end
end

dataBaseline = data;

end

