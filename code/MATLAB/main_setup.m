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

filePopUp = true;

loadLocal = false;

filterByDriver = false;

loadArgs = {'Effector','UAS_Chrimson_attp18_72F11',...
            'Protocol','r_LED30_30s2x15s30s',...
            'Dates',{'20170526','20170822'},...
            'DateFilter','between',...#
            'ExactMatch',{'Effector'}};
%             'Protocol','r*',...


% 'Driver',{'attp2',...
%     'GMR_SS01816',...
%     'GMR_SS01817',...
%     'GMR_SS01816n',...
%     'GMR_SS01817n',...
%     'GMR_SS01321',...
%     'GMR_SS01792',...
%     'GMR_SS02175',...
%     'GMR_SS04189',...
%     'GMR_SS00666',...
%     'GMR_SS00869',...
%     'GMR_SS01750',...
%     'GMR_SS04052',...
%     'GMR_SS04232',...
%     'GMR_SS04248',...
%     'GMR_SS43207'},...

set(0,'DefaultLegendAutoUpdate','off')

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
addpath(genpath(fullfile(params.directories.code,'MATLAB')))

params.directories.tempdir = fullfile(params.directories.master,"tempdir");
params.area_cutoff = 2.5;

contentfile = "contents.csv";
pipelines = ["mwt","choreography","salam","jaaba","jb"];
overwrite = false;

%% 

filetypes = {'salam','choreography','jaaba','jb'};

if ~loadLocal
    if filePopUp
        filepaths = choose_files(params,filetypes);
    else
        filepaths = select_files(params,filetypes,loadArgs{:});
    end

    setup_local_directory(params,filepaths);
end

%%
d = dir(fullfile(params.directories.tempdir,"**","*.mat"));
% filter out files without timestamp
fullpaths = fullfile({d.folder},{d.name});
timestamps = regexp(fullpaths,'(\d{8}_\d{6})', 'tokens', 'once');
filt = cellfun(@isempty,timestamps);
timestamps = [timestamps{:}];
fullpaths = fullpaths(~filt);
% sort file names by date
[timestamps,I] = sort(string(timestamps), 'descend');
filepaths = fullpaths(I);

clear I timestamps d

%% load genotype

% f = contains(filepaths,"attp2@UAS_Chrimson_attp18_72F11@r_LED05_30s2x15s30s");

filepaths = blacklist_filter(filepaths);

dc = dataClass(filepaths);

if filterByDriver
    fid = fopen('screen_genotypes.txt');
    myGenotypes = textscan(fid,'%s\n');
    myGenotypes = string([myGenotypes{:}]);
    fclose(fid);
    filterDriver = ismember(vertcat(dc.driver),myGenotypes);
    dc = dc(filterDriver);
end

dc = dc.update_figure_directories(params);