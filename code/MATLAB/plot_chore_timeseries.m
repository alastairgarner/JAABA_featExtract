%% plot_chore_timeseries.m

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
params = yaml.ReadYaml(fullfile(paramFile.folder,paramFile.name));
try
    cd(params.directories.master)
catch
    error('Non-existant directory specified - please update params.yaml file')
end
addpath(fullfile(params.directories.code,'MATLAB'))


%%

params.directories.data = '/media/alastair/Alastair_EH/Pipeline/choreography_results';

fprintf('\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n');

files = dataStruct.get_files(params);
data = dataStruct.make_dataStructs(files);

tstart = 30;
beh2plot = 'curve';
doi = "201908";
doi = '^(201905[23]|2019042)';
doi = '^2018120[123]';
toi = "115527";
% f = ([data.date] == doi);
f = regexp([data.date],doi);
f = ~cellfun(@isempty, f);
% f = startsWith([data.date],doi);
% f = ([data.exp_date] == doi) & ([data.exp_time] == toi);

data = data(f);

control_geno = 'attp2@U*';
data = data.sort_by_genotype(control_geno);
un_groups = unique([data.group])

behaviour_diffs = {};
for ii = 1:length(un_groups)
    grp = ii;
    temp = data;
    temp = temp([temp.group]==grp);
    
    if all(cellfun(@isempty,{temp.files_chore}))
        behaviour_diffs{ii} = zeros(1,1201);
        continue
    end
    
    for jj = 1:length(temp)
        temp(jj) = temp(jj).load_choredata();
    end
    % temp = data(10).load_choredata();

    bins = [0:1/25:120];
    for jj = 1:length(temp)
        if isempty(temp(jj).raw_chore)
            continue
        end

        times_full = temp(jj).raw_chore.et;
        [times,ia,ic] = unique(times_full);
    %     time_diff = times(1:10)-bins(1:10)
        time_diff = times-bins;
    %     time_diff(time_diff<0) = nan;
        [min_val,idx] = min(abs(time_diff));

        [frame_filter,et_adj] = ismember(ic,idx);
        et_adj(et_adj~=0) = bins(et_adj(et_adj~=0));
        et_adj(et_adj==0) = nan;
        max(min_val);

        temp(jj).raw_chore.et_adj = et_adj;
        temp(jj).raw_chore.frame_filter = ismember(ic,idx);
    end

    merged = merge_structure_fields([temp(:).raw_chore]);
    merged.uniid = get_unique_ids([merged.id]);

    filt = merged.frame_filter;
    [un_time,ia,ic] = unique(merged.et_adj(filt));
    [~,~,ids] = unique(merged.uniid(filt));
    [~,xlims] = min(abs((un_time'-[[-2 +8]+tstart])));
    [~,stim_start] = min(abs((un_time'-[[0]+tstart])));

    values_sparse = sparse(ids,ic,[merged.(beh2plot)(filt)],max(ids),max(ic));
    f = all(full(values_sparse(:,xlims) > 0),2);
    
%     figure()
%     hmap = imagesc(values_sparse(f,:),[0 50]);
% %     hmap = imagesc(values_sparse(f,:),[0 4]);
%     cmap = (cmocean('deep'));
%     cmap = (viridis);
%     cmap(1,:) = [.2 .2 .2];
%     colormap(cmap);
%     colorbar;
% 
%     ax = gca;
%     ax.XLim = xlims;
%     % ax.YLim = [0 50]+.5;
%     hold on
%     plot([stim_start,stim_start],ax.YLim,'r--','LineWidth',2)
%     hold off
%     pbaspect([1,2,1])
    
    [C,ia,ic] = unique(merged.uniid);
    ani_counts = accumarray(ic,1);
    ani_cell = mat2cell(merged.(beh2plot),1,ani_counts);
    ani_cell = cellfun(@(x) [0 diff(x)], ani_cell, 'UniformOutput', false);
    ani_vect = [ani_cell{:}];

    [Y,E] = discretize(merged.et,[0:.1:ceil(max(merged.et))]);
    behaviour_diffs{ii} = accumarray(Y',ani_vect,[],@mean)';
    
    %%
%     savename = strcat(temp(1).get_full_genotype,'@',temp(1).get_full_protocol,'@fig_heatplot');
%     savepath = fullfile('.','figures',savename);
%     print(savepath,'-dpdf','-painters','-fillpage');
    close
end
maxlen = max(cellfun(@length, behaviour_diffs));
cs = cellfun(@(x) padarray(x,[0 maxlen-length(x)],nan,'post'), behaviour_diffs,'UniformOutput',false);
% behaviour_diffs = cellfun(@(x) x(1:maxlen), behaviour_diffs,'UniformOutput',false);
diffs_mat = vertcat(cs{:});


%%

groups = unique([string([data.group])',data.get_full_genotype',data.get_full_protocol',[data.driver]'],'rows','stable')

idx = [7,6,4];
[cb] = cbrewer('qual', 'Dark2', length(idx), 'pchip');
cb = vertcat([0,0,0],cb);

n = 1;
for ii = idx
%     figure()
%     plot([0:maxlen-1]/10,diffs_mat(idx(1),:),'Color',cb(1,:),'LineWidth',3)
%     xlim([28 38])
%     ylim([-.5 1.5])
    hold on
    plot([0:maxlen-1]/10,diffs_mat(ii,:),'Color',cb(n,:),'LineWidth',3)
    n = n+1;
end
ax = gca;
plot([tstart,tstart],ax.YLim,'r--','LineWidth',2)
hold off
legend(groups(idx,4))
xlim([28 38])
% xlim([43 53])
pbaspect([2,1,1])
% ylim([-.5 2])

savename = strcat('fig_curvediff');
savepath = fullfile('.','figures',savename);
print(savepath,'-dpdf','-painters','-fillpage');
close

