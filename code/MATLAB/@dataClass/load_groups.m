%% load_groups

% Other m-files required: 
% Subfunctions: 
% MAT-files required: 
% See also:

% Author: Alastair Garner
% email: alastairgarner@outlook.com
% Website: https://github.com/alastairgarner/
% Sep 2019; Last revision: 

%%
function dc_group = load_groups(obj,filter_date_tf)
    % filter by overlapping dates (optional)
    if filter_date_tf
        [~,driv_grp] = obj.get_group_numbers('driver');
        [~,date_grp] = obj.get_group_numbers('date');
        [C,~,~] = unique([driv_grp,date_grp],'rows');
        expsPerDate = accumarray(C(:,2),1,[max(C(:,2)),1]);
        idx = find(expsPerDate == max(driv_grp));

        f = ismember(date_grp,idx);
        if sum(f) ~= 0
            obj = obj(f);
        else
            objFilter = ismember(date_grp,find(expsPerDate > 1));
            if ~sum(objFilter)
                objFilter = driv_grp > 1;
            end
            obj = obj(objFilter);
            fprintf("no control found \n")
            
%             obj = obj(driv_grp > 1);
%             fprintf("no control found \n")
        end
    end


    [n_prot,grp_prot] = obj.get_group_numbers('driver','effector','protocol');
    for xx = 1:n_prot
        f = grp_prot == xx;
        for zz = find(f)'
            obj(zz) = obj(zz).load_data();
        end

        obj(f) = obj(f).update_unique_ids();
        
        f = f & cellfun(@(x,y) numel(x)==numel(y) & numel(x)~=0,{obj.uniID},{obj.track_start})';
        if ~any(f)
            continue
        end
        
        dc_group(xx) = obj(f).merge_dataClasses();
        full_exp = strcat(obj(f).get_full_genotype(),'@',obj(f).get_full_protocol);
        fprintf("...loaded %s \n",full_exp(1));   
    end 
end