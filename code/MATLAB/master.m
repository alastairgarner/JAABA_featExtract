%% master.m

% Other m-files required: 
% Subfunctions: 
% MAT-files required: 
% See also:

% Author: Alastair Garner
% email: alastairgarner@outlook.com
% Website: https://github.com/alastairgarner/
% August 2019; Last revision: 

%% ------------- PARAMETER SETUP --------------

params.directories.current =            '/home/alastair/OneDrive/projects/JAABA_featExtract/code/MATLAB';
% params.directories.current =            fileparts(mfilename('fullpath'));

cd(params.directories.current); cd('../..'); 
params.directories.project =            pwd;
params.directories.data =               fullfile(params.directories.project,'data');
params.directories.chore_input =        fullfile(params.directories.project,'data','choreography_input');
params.directories.chore_output =       fullfile(params.directories.project,'data','choreography_results');
params.directories.salam =              fullfile(params.directories.project,'data','salam_feature_extraction');
params.directories.salam =              fullfile(params.directories.project,'data','jaaba_results');
params.directories.jb =              fullfile(params.directories.project,'data','jb_results');

params.filetypes.blobs =                '.blo*';
params.filetypes.choreography_txt =     'compiledChore.txt';
params.filetypes.choreography_dat =     '.dat';
params.filetypes.salam =                'animal_stats*.txt';
params.filetypes.jaaba =                {'trx.mat','scores*short.mat'};
params.filetypes.jb =                   'trx.mat';

addpath(fullfile(params.directories.project,'code','MATLAB'));

%% ------------- SCRIPT STARTS HERE --------------

fprintf('\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n');

files = dataStruct.get_files(params);
data = dataStruct.make_dataStructs(files);

n = 135;
try
    data(n) = data(n).load_blobsdata();
    data(n) = data(n).load_choreographydata();    
    data(n) = data(n).load_salamdata();
    data(n) = data(n).load_jaabadata();
    data(n) = data(n).load_jbdata();
end

pltr = larvaPlot('dataStruct',data(n));

pltr.salam
pltr.jaaba
pltr.jb

%% -------------- COMPILE JB -----------------

x = [vertcat(jb.tss_roll.bStart{:}),vertcat(jb.tss_roll.bEnd{:})]
aniID = [trx.numero_larva_num];
nbBehs = vertcat(trx.nb_action);
nbBehs = [nbBehs{:,6}];
y = repelem(aniID',nbBehs,1);
y = [y, y];

plot(x',y','k')
xlim([0 120])



et = {obj.data_jb.trx.t};
eTime = unique(vertcat(et{:}))';

fields_wanted = {["t"]};

jb = struct();
for ii = 1:length(fields_wanted)
    args = cellstr(fields_wanted{ii});
    len = length(args);
    try
        if len == 1
            jb.(args{end}) = {(obj.data_jb.trx.(args{1}))};
        elseif len == 2
            jb.(args{end}) = {obj.data_jb.trx.(args{1}).(args{2})};
        elseif len == 3
            jb.(args{end}) = {obj.data_jb.trx.(args{1}).(args{2}).(args{3})};
        end
    catch
        fprintf('... No field called jaaba_data.%s \n',strjoin(strcat(args,'.')));
    end
end




%% %------------- COMPILE SALAM ------------------


salam = data(n).compile_salam();


data(n).data_jaaba.allScores.tStart

spdat = sparse(145,3268);
dat = zeros(145,3268);

for ii = 1:145
    rline = jaaba.postprocessed{ii};
    rline(isnan(rline)) = 0;
    bS = jaaba.tStart(ii);
    bE = length(rline);
    dat(ii,1:bE) = rline;
    spdat(ii,1:bE) = rline;
end
rawdat = jaaba.postprocessed

whos spdat dat rawdat

spdat = sparse(145,3268);
dat = zeros(145,3268);
whos dat spdat




%% -------------- COMPILE JAABA -----------------

temp = data(n).data_jaaba

clear jaaba
jaaba.aniID = [temp.trx.id];
jaaba.tStart = temp.allScores.tStart;
jaaba.tEnd = temp.allScores.tEnd;
jaaba.timestamps = temp.timestamps;
jaaba.(temp.behaviorName).bStart = temp.allScores.t0s;
jaaba.(temp.behaviorName).bEnd = temp.allScores.t1s;

get_fields = {["behaviorName"];
        ["timestamps"];
        ["trx","id"];
        ["allScores","tStart"];
        ["allScores","tEnd"];
        ["allScores","t0s"];
        ["allScores","t1s"];
        ["allScores","postprocessed"]};
new_fields = {[""]};

jaaba = struct();
for ii = 1:length(get_fields)
    args = cellstr(get_fields{ii});
    len = length(args);
    try
        if len == 1
            jaaba.(args{end}) = [temp.(args{1})];
        elseif len == 2
            jaaba.(args{end}) = [temp.(args{1}).(args{2})];
        elseif len == 3
            jaaba.(args{end}) = [temp.(args{1}).(args{2}).(args{3})];
        end
    catch
        fprintf('... No field called jaaba_data.%s \n',strjoin(strcat(args,'.')));
    end
end

nBehav = cellfun(@length,jaaba.t0s);
y = repelem(jaaba.id,nBehav);
x = [[jaaba.t0s{:}];[jaaba.t1s{:}]];
x = jaaba.timestamps(x);

plot(x,[y;y],'k')
hold on

beh = 'turns';
nBehav = cellfun(@length,salam.(beh).bStart);
y = repelem([salam.aniID],nBehav);
x = [vertcat(salam.(beh).bStart{:}),vertcat(salam.(beh).bEnd{:})];

plot(x',[y;y]+.5,'m')
hold off




