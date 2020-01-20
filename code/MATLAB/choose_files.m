%% choose_files.m

% Other m-files required: 
% Subfunctions: 
% MAT-files required: 
% See also:

% Author: Alastair Garner
% email: alastairgarner@outlook.com
% Website: https://github.com/alastairgarner/
% Sep 2019; Last revision: 

%%
function filepaths = choose_files(params,filetypes)
    % find all compiled .mat files
    d = dir(fullfile(params.directories.data_compiled,"**","*.mat"));
    % filter out files without timestamp
    fullpaths = fullfile({d.folder},{d.name});
    timestamps = regexp(fullpaths,'(\d{8}_\d{6})', 'tokens', 'once');
    filt = cellfun(@isempty,timestamps);
    timestamps = [timestamps{:}];
    fullpaths = fullpaths(~filt);
    % sort file names by date
    [timestamps,I] = sort(string(timestamps), 'descend');
    fullpaths = fullpaths(I);
    % create string for list dialog
    splits = split(fullpaths,filesep);
    display_names = squeeze(splits(:,:,end));
    data_type = squeeze(splits(:,:,end-1));
    [display_names_unique,~,ic] = unique(display_names,'stable');

    chx = listdlg('ListString',string(display_names_unique));
    % filter file list by selection and data type
    filt_select = ismember(ic,chx);
    filt_type = any(string(data_type') == string(filetypes),2);

    filepaths = fullpaths(filt_select & filt_type)';
    
end