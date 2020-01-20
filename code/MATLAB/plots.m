%% plots.m

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

%% Load Parameters

parameterFile = 'params_AG_external.yaml'; 
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

searchterm = "201908"

f = dir("data_compiled/");
f = f(~[f.isdir]);
filelist = string({f.name})';
[indx,tf] = listdlg('ListString',filelist);

f = f(indx);

get_file_details(filelist)

function file_details = get_file_details(filelist)
    expr{1} = ['(?<date>\d{8})[_]'...
        '(?<time>\d{6})'];
    % expr{2} = ['[\@\\\/](?<driver>(?!20[\d+])[\w.]+)[@]'...
    %     '(?<effector>\w+)[\@\\\/]'];
    expr{2} = ['[\@\\\/](?<driver>(?!\d{8}_d{6})[\w.]+)[@]'...
        '(?<effector>\w+)[\@\\\/]'];
    expr{3} = ['[\@\\\/](?<rig>[t]\d{2})[\@\\\/]'];
    expr{4} = ['[\@\\\/](?<protocol1>[\w\d\_]+)'...
        '[#](?<protocol2>[\w\d\_]+)'...
        '[#](?<protocol3>[\w\d\_]+)'...
        '[#](?<protocol4>[\w\d\_]+)'];

    file_details = struct("date","",...
        "time","",...
        "driver","",...
        "effector","",...
        "protocol1","",...
        "protocol2","",...
        "protocol3","",...
        "protocol4","");
    
    file_details = repelem(file_details,numel(filelist),1);
    for ii = 1:numel(expr)
        matches = regexp(filelist,expr{ii},'names','once');
        filt = ~cellfun(@isempty,matches);
        matches = vertcat(matches{:});
        fnames = fieldnames(matches);
        for jj = 1:numel(fnames)
            field_vals = {matches.(fnames{jj})};
            [file_details(filt).(fnames{jj})] = field_vals{:};
        end
    end
end

fd = file_details;
list_select = strcat([fd.driver],'@',[fd.effector],'@',[fd.protocol1],'@',[fd.protocol2])';
[un_list_select,~,idx] = unique(list_select,'stable');

for ii = 1:numel(un_list_select)
    f = idx == ii;
    file_groups{ii} = list_select(f);
end
idx'


