%% main.m

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

search_term = "20190531_140148";

%% Load Parameters


parameterFile = 'params_AG.yaml'; 
% parameterFile = 'params_AG.yaml'; 

paramFile = dir(fullfile('.','**',[parameterFile,'*']));
params = yaml.ReadYaml(fullfile(paramFile.folder,paramFile.name));
try
    cd(params.directories.master)
catch
    error('Non-existant directory specified - please update params.yaml file')
end
addpath(fullfile(params.directories.code,'MATLAB'))

%%

fprintf('\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n');

file_structure = dataContainer.get_files(params);
fullpaths = fullfile({file_structure.folder},{file_structure.name});
search_filt = contains(fullpaths,search_term);
fullpaths = fullpaths(search_filt);
timestamps = regexp(fullpaths,'(\d{8}_\d{6})', 'tokens', 'once');
timestamps = string([timestamps{:}]);
[unique_timestamps,~,indicies] = unique(timestamps);

for xx = 1:numel(unique_timestamps)
    idx = indicies == xx;
    dc = dataContainer(fullpaths(idx));
    data = [];
    for ii = 1:numel(dc)
        dc(ii) = dc(ii).load_data();
        comp = dc(ii).compile_data();
        dc(ii).raw_data = [];
        data = vertcat(data,comp);
    end
    if all(arrayfun(@(x) isempty(x.aniID),data))
        continue
    end
    dataContainer.save_data_batch(data,params)
end

txt = jsonencode(data)
txt(1:1000)
saveJSONfile(data,'test.json')

%%
expr1 = ['(?<date>\d{8})[_]'...
            '(?<time>\d{6})'];
expr2 = ['[\@\\\/](?<driver>(?![\d{8}])[\w.]+)[@]'...
            '(?<effector>\w+)[\@\\\/]'];
expr3 = ['[\@\\\/](?<rig>[t]\d{2})[\@\\\/]'];
expr4 = ['[\@\\\/](?<protocol1>[\w\d\_]+)'...
            '[#](?<protocol2>[\w\d\_]+)'...
            '[#](?<protocol3>[\w\d\_]+)'...
            '[#](?<protocol4>[\w\d\_]+)'];
%%

searchterm = "201908"

f = dir("data_compiled/");
{f.name}








