%% main_analysis.m

% Other m-files required: 
% Subfunctions: 
% MAT-files required: 
% See also:

% Author: Alastair Garner
% email: alastairgarner@outlook.com
% Website: https://github.com/alastairgarner/
% Sep 2019; Last revision: 

%% User settings
clear all; clc;
fprintf('\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n');

%% Load Parameters

parameterFile = 'params_AG_external.yaml'; 
% parameterFile = 'params_AG_external.yaml'; 

paramFile = dir(fullfile('.','**',[parameterFile,'*']));
params = yaml.ReadYaml(fullfile(paramFile.folder,paramFile.name));
try
    cd(params.directories.master)
catch
    error('Non-existant directory specified - please update params.yaml file')
end
addpath(fullfile(params.directories.code,'MATLAB'))

params.directories.tempdir = fullfile(params.directories.master,"tempdir");

contentfile = "contents.csv";
pipelines = ["mwt","choreography","salam","jaaba","jb"];
overwrite = false;

%% 

filetypes = {'salam','choreography'};

filepaths = choose_files(params,filetypes);

setup_local_directory(params,filepaths);

deets = dataClass.parse_filepaths(string(filepaths))

obj = dataClass(filepaths(1));

dc = dataClass(filepaths(1:2));

for ii = 1:numel(dc)
    dc(ii) = dc(ii).load_data();
end


%% combine instances

[C,ia,ic] = unique(dc.get_full_timestamp);

for ii = 1:numel(C)
    f = ii == ic;
    fdnames = fieldnames(dc);
    temp = struct();
    for jj = 1:numel(fdnames)
        if strmatch(fdnames{jj},'aniID')
            temp.(fdnames{jj}) = unique([dc(f).(fdnames{jj})],'stable');
            
        elseif strmatch(fdnames{jj},'track_start')
            [~,~,idx] = unique([dc(f).aniID],'stable');
            temp.(fdnames{jj}) = accumarray(idx,[dc(f).(fdnames{jj})]',[],@min)';
            
        elseif strmatch(fdnames{jj},'track_end')
            [~,~,idx] = unique([dc(f).aniID],'stable');
            temp.(fdnames{jj}) = accumarray(idx,[dc(f).(fdnames{jj})]',[],@max)';
            
        elseif strmatch(fdnames{jj},'behaviour')
            cells = {dc(f).(fdnames{jj})};
            goodstructs = cellfun(@(x) numel(fieldnames(x))>0, cells);
            temp.(fdnames{jj}) = vertcat(cells{goodstructs});
        elseif strmatch(fdnames{jj},'timeseries')
            ids = []; cell_arr = [{}];
            for xx = 1:numel(dc)
                ids = [ids,[dc(xx).(fdnames{jj}).id]];
                cell_arr = [cell_arr {struct2cell(dc(xx).(fdnames{jj}))}]
            end
            [un,~,idx] = unique(ids,'stable')
            padarray(struct2cell(dc(1).timeseries)
            cells = {dc.timeseries};
            [cells{:}.id]
            
            goodstructs = arrayfun(@(x) isfield(x.timeseries,'id'),dc);
            cells = arrayfun(@(x) [x.timeseries.id], dc(goodstructs), 'UniformOutput', false);
            [cells{:}]
        elseif strmatch(fdnames{jj},'raw_data')
            
        else
            temp.(fdnames{jj}) = unique([dc(filt).(fdnames{jj})])
        end
    end
end





