%% get_behaviour_data

% Other m-files required: 
% Subfunctions: 
% MAT-files required: 
% See also:

% Author: Alastair Garner
% email: alastairgarner@outlook.com
% Website: https://github.com/alastairgarner/
% Sep 2019; Last revision: 

%%
function data_array = get_behaviour_data(obj,behaviour,params)
%     behaviour = "rolls"
    
    [n_exp,grp_exp] = obj.get_group_numbers('driver','effector','protocol');
    
    clear data_array
    for jj = 1:10%n_exp
        temp = [];
        
        f = grp_exp == jj;
        dA = obj(f).load_groups();
        dA = dA.filter_by_size(params.area_cutoff);
        
        b = find(strcmp({dA.behaviour.behaviour},behaviour));
        
        temp = dA.behaviour(b);
        
        temp.duration = temp.end - temp.start;
        temp.n_behaviours = accumarray(temp.uniID',~isnan(temp.start),[],@sum)';
        temp.genotype_id = repelem(jj,1,numel(temp.uniID));
        temp.behaviour = behaviour;
        temp.driver = dA.driver;
        temp.n_timestamps = size(unique([[dA.date];[dA.time]]','rows'),1);
        
        [r,c] = find(dA.uniID' == temp.uniID);
        temp.timestamp_index = dA.timestamp_index(r);
        
        [tf,~] = ismember([dA.timeseries.uniID],temp.uniID);
        temp.area = cellfun(@mean,{dA.timeseries(tf).area});
        temp.area_ids = dA.uniID(tf);
        
        data_array(jj) = temp;
    end
end