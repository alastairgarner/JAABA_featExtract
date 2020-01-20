%% setup_local_directory.m

% Other m-files required: 
% Subfunctions: 
% MAT-files required: 
% See also:

% Author: Alastair Garner
% email: alastairgarner@outlook.com
% Website: https://github.com/alastairgarner/
% Sep 2019; Last revision: 

%%
function setup_local_directory(params,filepaths)
    % reset temp directory
    if isdir(params.directories.tempdir)
        rmdir(params.directories.tempdir,'s')
    end
    mkdir(params.directories.tempdir)
    % make local directory (for speedier loading of data)
    fullpaths_local = strrep(filepaths,params.directories.data_compiled,params.directories.tempdir);
    tempfolders = unique(cellfun(@fileparts,fullpaths_local, 'UniformOutput', false));
    for ii = 1:numel(tempfolders)
        mkdir(tempfolders{ii})
    end
    % copy files to local folder
    fprintf('setting up local directory...\n')
    for ii = 1:numel(filepaths)
        copyfile(filepaths{ii},fullpaths_local{ii});
    end
    fprintf('local directory ready.\n')
end

