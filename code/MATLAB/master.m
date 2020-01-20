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

clear all; clc;

parameterFile = 'params_AG.yaml'; 

paramFile = dir(fullfile('.','**',[parameterFile,'*']));
params = yaml.ReadYaml(fullfile(paramFile.folder,paramFile.name))
try
    cd(params.directories.master)
catch
    error('Non-existant directory specified - please update params.yaml file')
end
addpath(fullfile(params.directories.code,'MATLAB'))

% params.filetypes.blobs =                '.blo*';
% params.filetypes.choreography_txt =     'compiledChore.txt';
% params.filetypes.choreography_dat =     '.dat';
% params.filetypes.salam =                'animal_stats*.txt';
% params.filetypes.jaaba =                {'trx.mat','scores*short.mat'};
% params.filetypes.jb =                   'trx.mat';


%% ------------- SCRIPT STARTS HERE --------------

fprintf('\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n');

files = dataStruct.get_files(params);
data = dataStruct.make_dataStructs(files);

doi = "20190423";
toi = "115527";
f = ([data.date] == doi);
% f = ([data.exp_date] == doi) & ([data.exp_time] == toi);

dS = data(f);
dS.load_compile_save(params);

clear dS

files = dir(fullfile(params.directories.dataprocessed,'*.mat'))
fullcontents = [];
for ii = 1:length(files)
    fname = fullfile(files(ii).folder,files(ii).name);
    cont = load(fname,'content');
    fullcontents = [fullcontents; struct2array(cont.content)];
end

imagesc(fullcontents)

%% Load Data

files = dir(fullfile(params.directories.dataprocessed,'*.mat'));
data = [];
for ii = 1:length(files)
    fname = fullfile(files(ii).folder,files(ii).name);
    temp = load(fname,'temp');
    data = [data; temp.temp];
end

control_geno = 'attp2@U*';

data = data.sort_by_genotype(control_geno);

temp = data([data.iscontrol])

% data.get_group(1)

grp = 1

temp = data;

temp = temp([temp.group]==grp);
temp = temp.get_behaviour_data('rolls');
temp = merge_structure_fields(temp);
temp.uniID = get_unique_ids([temp.aniID]);

beh_counts = cellfun(@length, temp.behaviour.bStart);
time_sta = [temp.behaviour.bStart{:}];
time_end = [temp.behaviour.bEnd{:}];
id = [temp.uniID];

ids = repelem(id,beh_counts');

plot([time_sta; time_end],[ids; ids],'k');

temp(1).behaviour

%% curve timeseries

files = dataStruct.get_files(params);
data = dataStruct.make_dataStructs(files);

doi = "20190423";
toi = "115527";
f = ([data.date] == doi);
% f = ([data.exp_date] == doi) & ([data.exp_time] == toi);
dS = data(f);

control_geno = 'attp2@U*';
dS = dS.sort_by_genotype(control_geno);

grp = 4
temp = dS;
temp = temp([temp.group]==grp);

temp = temp(1).load_choredata()
[un_time,ia,ic] = unique(temp.raw_chore.et);
[~,~,ids] = unique(temp(1).raw_chore.id);
[~,xlims] = min(abs((un_time-[29 40])));
[~,stim_start] = min(abs((un_time-[30])));

values_sparse = sparse(ids,ic,[temp.raw_chore.speed],max(ids),max(ic));
hmap = imagesc(values_sparse,[0 60]);
hmap = imagesc(values_sparse,[0 4]);
ax = gca;
ax.XLim = xlims;
hold on
plot([stim_start,stim_start],ax.YLim,'r--')
hold off


cmap = (cmocean('deep'))
cmap(1,:) = [1 1 1];
colormap(cmap)
colorbar

cmap = cmocean('algae',12)

%% Ridgeline Plot


files = dataStruct.get_files(params);
data = dataStruct.make_dataStructs(files);

doi = "20190423";
toi = "115527";
f = ([data.date] == doi);
dS = data(f);

control_geno = 'attp2@U*';
dS = dS.sort_by_genotype(control_geno);

grp = 4
temp = dS;
temp = temp([temp.group]==grp);
temp = temp(1).load_choredata();

mat = [[temp.raw_chore.et]';[temp.raw_chore.id]';[temp.raw_chore.curve]'];
filt = mat(1,:)>=28 & mat(1,:)<= 35;
mat = mat(:,filt);
idx = accumarray(mat(2,:)',1);
idx = idx(find(idx));
curves = mat2cell(mat,3,idx);

num_line = length(idx);
cmap = flipud(viridis(num_line));
hold on
for jj = 1:num_line
    cv = curves{jj}(3,:)./40;
%     cv = [0 diff(curves{jj}(3,:))]./10;
    y = (-jj)+[cv,0,0];
    et = curves{jj}(1,:);
    x = [et,max(et),min(et)];
    patch(x,y,'red','EdgeColor','white','LineWidth',0.1,'FaceColor',cmap(jj,:));
end
hold off
set(gcf,'PaperOrientation','landscape');
pbaspect([1,1,1])

ylim([-num_line,0]);
set(gca,'YTickLabels',[],'YTick',[]);

%%
temp = data(700);
temp = temp.load_blobsdata;
ids = [temp.raw_blobs.aniID];
[C,ia,ic] = unique(ids,'stable');
id_counts = accumarray(ic,1);

blobsData = temp.raw_blobs;
blobsData = structfun(@(x) mat2cell(x,1,id_counts'),blobsData,'UniformOutput',false);
fields = fieldnames(blobsData);

aniData = struct();
for ii = 1:length(fields)
    for jj = 1:length(blobsData.(fields{ii}))
        aniData(jj,1).(fields{ii}) = blobsData.(fields{ii}){jj};
    end
end
duration = cellfun(@range, {aniData.time});
filt = duration >= params.config_choreography.t;
aniData = aniData(filt);

id = 1;
fr = 1;
[x,y] = update_plot_xy(aniData,id,fr);



%%
behaviour = 'roll_tss';

bData = data(1).get_behaviour_data(behaviour);
[tracked,uniIDs] = data(1).get_tracked(behaviour);
[behaved,uniIDs] = data(1).get_behaved(behaviour);

plot(behaved./tracked)

%%

[cb] = cbrewer('qual', 'Set3', 12, 'pchip');
n = 250;
% set a random number generator seed for reproducible results
rng(123)

d{1} = [exprnd(5, 1, n) + 15]';
d{2} = [(randn(1, n) *5) + 20]';

f7 = figure();
h1 = raincloud_plot(d{1}, 'box_on', 1, 'color', cb(1,:), 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
     'box_col_match', 0, 'lwr_bnd', 10, 'cloud_edge_col',cb(1,:));
h2 = raincloud_plot(d{2}, 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35,...
     'box_col_match', 0, 'lwr_bnd', 10, 'cloud_edge_col',cb(4,:));
h1{1}.EdgeAlpha = 0.5;
h2{1}.EdgeAlpha = 0.5;


legend([h1{1} h2{1}], {'Group 1', 'Group 2'});
title(['Figure M7' newline 'A) Dodge Options Example 1']);
set(gca,'XLim', [0 40], 'YLim', [-.075 .15]);
box off



%% plot beh starts

behs = {'salam','rolls';...
    'jaaba','roll_beg_short';...
    'jb','roll_tss'};

dodge = [.15,.35,.55]

[cb] = cbrewer('qual', 'Set3', 12, 'pchip');
cb = cb([1,5,6],:)

xs = {}; clear h
for ii = 1:size(behs,1)
    aniIDs = pltr.(behs{ii,1}).aniID;
    temp = pltr.(behs{ii,1}).(behs{ii,2});
    
%     x = [temp.bStart{:}];
%     filt = x < 75;
%     x = x(filt);
    
    fx = @(x) x(find([x>30 & x<45],1));
    bFirst = cellfun(fx, temp.bStart,'UniformOutput', false);
    x = [bFirst{:}];
%     
    h{ii} = raincloud_plot(x, 'box_on', 1, 'color', cb(ii,:), 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', dodge(ii), 'dot_dodge_amount', dodge(ii),...
     'box_col_match', 0, 'lwr_bnd', 10, 'cloud_edge_col',cb(ii,:));
    h{ii}{1}.EdgeAlpha = 0.5;
    
    xs = [xs {x}];
   
end
set(gca,'XLim', [28 45], 'YLim', [-.25 .60]);
pbaspect([2,1,1])
box off

cl = num2cell(cb(1:3,:),2)';
cl = cb(1:3,:);
h   = rm_raincloud(xs', cb(1,:),1)


%% PLOT ALL
blobs = blobClass.blobClass_batch(data(n).data_blobs);

temp = blobs.get_animals([4]);
% temp = blobs

aniID = [temp.aniID];
uniID = cumsum([1,diff(aniID)~=0]);
frames = [temp.frames];
x_pts = [temp.x_pt];
y_pts = [temp.y_pt];
npts = [temp.npts];
outlines = [temp.outline];

rx = range(x_pts);
ry = range(y_pts);
mx = round(min(x_pts)+[rx/2]);
my = round(min(y_pts)+[ry/2]);
rmax = round(max(rx,ry)/2)*1.1;

p = plot(nan(max(uniID)),'k');

xlim([-rmax rmax]+mx)
ylim([-rmax rmax]+my)
pbaspect([1,1,1])

totframes = [1:max([temp.frames])]';

id_last = uniID;
tic;
for ii = 1:max(totframes)
%     tic;
    f = find(frames == ii);
    ids = uniID(f);
    
%     x = x_pts(f);
%     y = y_pts(f);
%     npt = npts(f);
%     kerr = outlines(f);
%     [xx,yy] = cellfun(@(k,n,x,y) unpackKerr_xy(k,n,x,y), kerr, num2cell(npt), num2cell(x), num2cell(y), 'UniformOutput', false);
    
%     [xx,yy] = cellfun(@(k,n,x,y) unpackKerr_xy(k,n,x,y),...
%         outlines(f),...
%         num2cell(npts(f)),...
%         num2cell(x_pts(f)),...
%         num2cell(y_pts(f)),...
%         'UniformOutput', false);
    
    [xx,yy] = unpackKerr_xy_batch(outlines(f),...
        npts(f),...
        x_pts(f),...
        y_pts(f));
    
    set(p(ids),{'XData'},xx',{'YData'},yy');
    
    dd = setdiff(id_last,ids);
    if any(find(dd))
        set(p(dd),'XData',nan,'YData',nan);
    end
    
    fprintf('\n %.2f', round(ii/30,2))    
    pause(max(1/60-toc,0.0001))
    tic

    id_last = ids;
    
     % optimising
%     fprintf('\n %.2f vs %.4f', round(ii/30,2), toc)
%     tic
%     pause(0.0001)
end


%%

deet = 26;

timestamps = data(n).data_jaaba.timestamps;
xlims = [data(n).data_jaaba.trx.id] == deet ;
timestamps(data(n).data_jaaba.allScores.tStart(xlims))

xlims = [data(n).data_jb.trx.numero_larva_num] == deet;
min(data(n).data_jb.trx(xlims).t)

data(n).data_featextract.peran


pltr = larvaPlot('dataStruct',data(n));

[pltr.salam.rolls.bStart{:}]
[pltr.jaaba.roll_beg_short.bStart{:}]
[pltr.jb.roll_tss.bStart{:}]

types = {'salam','jaaba','jb'};
colcode = {'k','r','b'};
for ii = 1:length(types)
    hold on
    x = [pltr.(types{ii}).tStart;pltr.(types{ii}).tEnd];
    y = repmat([pltr.(types{ii}).aniID],2,1) + ([ii-1]*0.2);
    plot(x,y,colcode{ii});
end


types = {'salam','jaaba','jb'};
behs = {'rolls','roll_beg_short','roll_tss'}
colcode = {'r','k','g'};
for ii = 1:length(types)
    hold on
    y = [pltr.(types{ii}).aniID];
    cnts = cellfun(@length,pltr.(types{ii}).(behs{ii}).bStart);
    y = repelem(y,2,cnts) + ([ii-1]*0.1);
    
    x = [[pltr.(types{ii}).(behs{ii}).bStart{:}]; [pltr.(types{ii}).(behs{ii}).bEnd{:}]];
    
    plot(x,y,colcode{ii},'LineWidth',1);
end

types = {'jb','salam'};
behs = {'roll_tss','rolls'}
colcode = {'r','k','g'};
for ii = 1:length(types)
    hold on
    y = [pltr.(types{ii}).aniID];
    cnts = cellfun(@length,pltr.(types{ii}).(behs{ii}).bStart);
    y = repelem(y,2,cnts) + ([ii-1]*0.1);
    
    x = [[pltr.(types{ii}).(behs{ii}).bStart{:}]; [pltr.(types{ii}).(behs{ii}).bEnd{:}]];
    
    plot(x,y,colcode{ii},'LineWidth',2);
end


%% Sparse Plotter

[[pltr.salam.rolls.bStart{:}]; [pltr.salam.rolls.bEnd{:}]]

beh = 'peran';

tmin = min(pltr.salam.tStart)
tmax = max(pltr.salam.tEnd)
vect = [0:0.001:tmax];

bStarts = [pltr.salam.(beh).bStart{:}];
bEnds = [pltr.salam.(beh).bEnd{:}];
aniIDs = pltr.salam.aniID;
cnts = cellfun(@length,pltr.salam.(beh).bStart);
aniIDs = repelem(aniIDs,1,cnts);

f = ~isnan(bStarts);
iStarts = round(bStarts(f)*1000);
iEnds = round(bEnds(f)*1000);
iIDs = aniIDs(f);

spStarts = sparse(iIDs,iStarts,1,max(iIDs),length(vect));
spEnds = sparse(iIDs,iEnds,-1,max(iIDs),length(vect));
spBoth = spStarts + spEnds;
spBoth = cumsum(spBoth,2);

imagesc(spBoth)
colormap(gray)

whos spStarts spEnds spBoth

empt = sparse(max(iIDs),length(vect));
empt(iIDs',iStarts') = 1;
empt([2,2],[2,3]) = 

matStart = sparse([pltr.salam.rolls.bStart{:}]' == vect);
matEnd = sparse([pltr.salam.rolls.bEnd{:}]' == vect)*-1;
matBoth = matStart + matEnd;
matCumsum = cumsum(matBoth,2);

imagesc(matCumsum)

whos vect mat matStart matEnd matBoth matCumsum

full(matCumsum(1,:))

pltr.salam.rolls


%% -------------- COMPILE JB -----------------






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
timestamps = temp.timestamps;

clear jaaba
jaaba.aniID = [temp.trx.id];
jaaba.tStart = timestamps(temp.allScores.tStart);
jaaba.tEnd = timestamps(temp.allScores.tEnd);
jaaba.(temp.behaviorName).bStart = cellfun(@(x) timestamps(x), temp.allScores.t0s, 'UniformOutput',false);
jaaba.(temp.behaviorName).bEnd = cellfun(@(x) timestamps(x), temp.allScores.t1s, 'UniformOutput',false);
jaaba.timestamps = timestamps;

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




