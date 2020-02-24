%% update_blacklist.m

% Other m-files required: 
% Subfunctions: 
% MAT-files required: 
% See also:

% Author: Alastair Garner
% email: alastairgarner@outlook.com
% Website: https://github.com/alastairgarner/
% Sep 2019; Last revision: 

%%
function update_blacklist(timestamps)

    pr_check = regexp(timestamps,'\d{8}_\d{6}','match','once');

    if any(cellfun(@isempty, pr_check))
        fprintf('bad timestamp \n')
        return
    end

    blacklist_file = "./blacklist.txt";
    if ~isfile(blacklist_file)
        lines = [];
    else
        fid = fopen(blacklist_file,'r');
        lines = textscan(fid,'%s\n');
        fclose(fid);
        lines = lines{:};
    end

    out = [timestamps;lines];
    out = unique(out,'stable');

    fid = fopen(blacklist_file,'w');
    fprintf(fid,'%s\n',out{:});
    fclose(fid);

end

