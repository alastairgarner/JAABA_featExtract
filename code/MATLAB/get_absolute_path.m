%% get_absolute_path.m

% Other m-files required: 
% Subfunctions: 
% MAT-files required: 
% See also:

% Author: Alastair Garner
% email: alastairgarner@outlook.com
% Website: https://github.com/alastairgarner/
% Sep 2019; Last revision: 

%%
function abs_paths = get_absolute_path(rel_paths)
    rel = string(rel_paths);
    
    current_dir = pwd;
    absol = [];
    for ii = 1:numel(rel)
        cd(rel(ii))
        absol = vertcat(absol,string(pwd));
        cd(current_dir)
    end
    abs_paths = absol;
    
end
