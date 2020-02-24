%% load_bygenotype.m

% Other m-files required: 
% Subfunctions: 
% MAT-files required: 
% See also:

% Author: Alastair Garner
% email: alastairgarner@outlook.com
% Website: https://github.com/alastairgarner/
% Sep 2019; Last revision: 

%%
function dc_genotype = load_bygenotype(obj,genotype_number)
    geno = genotype_number;
    
    f = [obj.group_genotype] == geno;
    temp = obj(f);

    for ii = 1:numel(temp)
        temp(ii) = temp(ii).load_data();
    end

    %%

    dca = temp;

    temp = temp.update_unique_ids();

    dc_genotype = temp.merge_dataClasses();
    full_exp = strcat(dca.get_full_genotype,'@',dca.get_full_protocol);
    fprintf("...loaded %s \n",full_exp(1));    
end