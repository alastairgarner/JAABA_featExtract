%% main_process.m

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

% search_term = "/2017";

%% Load Parameters


parameterFile = 'params_AG.yaml'; 
% parameterFile = 'params_AG_external.yaml';

paramFile = dir(fullfile('.','**',[parameterFile,'*']));
params = yaml.ReadYaml(fullfile(paramFile.folder,paramFile.name));
try
    cd(params.directories.master)
catch
    error('Non-existant directory specified - please update params.yaml file')
end
addpath(fullfile(params.directories.code,'MATLAB'))

contentfile = "contents.csv";
pipelines = ["mwt","choreography","salam","jaaba","jb"];
overwrite = false;

%%

fprintf('\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n');

file_structure = dataContainer.get_files(params);
fullpaths = fullfile({file_structure.folder},{file_structure.name});
if exist('search_term')
    search_filt = contains(fullpaths,search_term);
    fullpaths = fullpaths(search_filt);
end
if ~overwrite
    fullpaths = dataClass.filterby_contents_file(fullpaths,"contents.csv",pipelines,params);
end
timestamps = regexp(fullpaths,'(\d{8}_\d{6})', 'tokens', 'once');
filt = ~cellfun(@isempty, timestamps);
fullpaths = fullpaths(filt);
timestamps = string([timestamps{:}]);
[unique_timestamps,~,indicies] = unique(timestamps,'stable');

%% Version 1
% 
% for xx = 1:numel(unique_timestamps)
%     idx = indicies == xx;
%     timestamps = regexp(fullpaths(idx),'(\d{8}_\d{6})', 'tokens', 'once');
%     if numel(unique([timestamps{:}])) > 1
%         fprintf("\n too many damn timestamps")
%         break
%     end
%     dc = dataClass(fullpaths(idx));
%     data = [];
%     for ii = 1:numel(dc)
%         dc(ii) = dc(ii).load_data();
%         comp = dc(ii).compile_data();
%         dc(ii).raw_data = [];
%         data = vertcat(data,comp);
%     end
%     
%     if all(arrayfun(@(x) isempty(x.aniID),data))
%         continue
%     end
%     outnames = data.save_data(params);
%     
%     contentfile = "contents.csv";
%     pipelines = ["mwt","choreography","salam","jaaba","jb"];
%     dataClass.update_contents_file(contentfile,outnames,pipelines,params);
% end
    
%% Version 2

% fpath = '/home/alastair/projects/JAABA_featExtract/downloads/20200824_0040/jb-reduced/t93/attp2@UAS_Chrimson_attp18_72F11/r_LED30_30s2x15s30s#n#n#n@100/20181130_170329/trx.mat'
% dc = dataClass({fpath});
% dc = dc.load_data()
% dc.compile_data()

dc = dataClass(fullpaths);

contentfile = "contents.csv";
pipelines = ["mwt","choreography","salam","jaaba","jb"];
for ii = 1:numel(dc)
    fprintf("%d/%d \n",ii,numel(dc));
    
    if isnan(double(dc(ii).date))
        fprintf("bad timestamp \n")
        continue
    end
    
    dc(ii) = dc(ii).load_data();
    comp = dc(ii).compile_data();
    dc(ii).raw_data = [];

    outnames = comp.save_data(params);
    
    dataClass.update_contents_file(contentfile,outnames,pipelines,params);
end


% find(cellfun(@(x) any(contains(x,"jaaba")), {dc.filepath}))



