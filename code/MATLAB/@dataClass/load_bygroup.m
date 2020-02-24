%% load_bygroupm

% Other m-files required: 
% Subfunctions: 
% MAT-files required: 
% See also:

% Author: Alastair Garner
% email: alastairgarner@outlook.com
% Website: https://github.com/alastairgarner/
% Sep 2019; Last revision: 

%%
function dc_genotype = load_bygroup(obj)
    for ii = 1:numel(obj)
        obj(ii) = obj(ii).load_data();
    end

    %%

    obj = obj.update_unique_ids();

    dc_genotype = obj.merge_dataClasses();
    full_exp = strcat(obj.get_full_genotype,'@',obj.get_full_protocol);
    fprintf("...loaded %s \n",full_exp(1));    
end