%% blacklist_filter.m

% Other m-files required: 
% Subfunctions: 
% MAT-files required: 
% See also:

% Author: Alastair Garner
% email: alastairgarner@outlook.com
% Website: https://github.com/alastairgarner/
% Sep 2019; Last revision: 

%%
function filepaths = blacklist_filter(filepaths)

    blacklist_file = "./blacklist.txt";
    if ~isfile(blacklist_file)
        lines = [];
    else
        fid = fopen(blacklist_file,'r');
        lines = textscan(fid,'%s\n');
        fclose(fid);
        lines = lines{:};
    end
    
    f = ~contains(filepaths,lines');
    filepaths = filepaths(f);
end

