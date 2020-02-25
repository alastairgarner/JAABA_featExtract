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

filepaths = choose_files(params,filetypes);

setup_local_directory(params,filepaths);

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

dc = dc.update_figure_directories(params);

%% individual plots

[n_groups,group_num] = dc.get_group_numbers('driver','effector','protocol');

for ii = 1:n_groups
    idx = group_num == ii;
    dct = dc(idx);
    
    %
    dA = dct.load_groups(true);
    dA = dA.filter_by_size(params.area_cutoff);
    %
    [ax,figure_path] = dA.plot_ridgeline_timeseries(params,'crabspeed',25,45,5/12);
    fix_axes(gcf,7,false)
    dataClass.save_figure_catch(ax,figure_path);
    close
    %
    [ax,figure_path] = dA.plot_ridgeline_timeseries(params,'curve',25,45,5/100);
    fix_axes(gcf,7,false)
    dataClass.save_figure_catch(ax,figure_path);
    close
    %
    [ax,figure_path] = dA.plot_ridgeline_timeseries(params,'speed',25,45,5/12);
    fix_axes(gcf,7,false)
    dataClass.save_figure_catch(ax,figure_path);
    close
    %
    [ax,figure_path] = dA.plot_ridgeline_timeseries(params,'midline',25,45,5/5);
    fix_axes(gcf,7,false)
    dataClass.save_figure_catch(ax,figure_path);
    close
    
    %
    behaviours = {"peran","rolls"};
    [ax,figure_path] = dA.plot_ethogram(params,behaviours,25,50,0);
    fix_axes(gcf,7,false)
    dataClass.save_figure_catch(ax,figure_path);
    close
%     %
%     behaviours = {"peran","turns"};
%     [ax,figure_path] = dA.plot_ethogram(params,behaviours,25,50,0);
%     print(figure_path,'-dpdf','-painters','-fillpage');
%     close
%     %
%     behaviours = {"peran","rolls","roll_tss","roll_beg_short"};
%     [ax,figure_path] = dA.plot_ethogram(params,behaviours,25,50,0.2);
%     dataClass.save_figure_catch(ax,figure_path);
%     close
end



%% paired plots

cmap = cbrewer('qual','Dark2',8);
cmap(1,:) = [0.2 0.2 0.2];
% cmap = vertcat([0 0 0],cols);

[n_prot,grp_prot] = dc.get_group_numbers('effector','protocol');

for ii = 1:n_prot
    fprintf("%d/%d \n",ii,n_prot);
    datagrp = [];
    
    f1 = grp_prot == ii;
    dcS = dc(f1);
    [n_exp,grp_exp] = dcS.get_group_numbers('driver','effector','protocol');
    
    fc = strcmp([dcS.driver],"attp2");
    if any(fc)
        controldata = dcS(fc);
        controlnum = unique(grp_exp(fc));
    else
        controldata = [];
        controlnum = 0;
    end
    
    % plot loop
    for jj = 1:n_exp
        if jj == controlnum
            continue
        end
            
        f2 = grp_exp == jj;
        datagrp = vertcat(controldata,dcS(f2));
        
        dataArray = datagrp.load_groups(true);
        dataArray = dataArray.filter_by_size(params.area_cutoff);

        %%
%         figure('Resize','off')
        subplot(2,2,1)
        [ax,figure_path] = dataArray.plot_timeseries(params,"curve",28,40,[0 50],"norm",cmap);
        subplot(2,2,3)
        [ax,figure_path] = dataArray.plot_timeseries(params,"crabspeed",28,40,[0 2.5],"norm",cmap);
        subplot(2,2,2)
        [ax,figure_path] = dataArray.plot_timeseries(params,"curve",28,40,[-1 4],"diff",cmap);
        subplot(2,2,4)
        [ax,figure_path] = dataArray.plot_timeseries(params,"crabspeed",28,40,[-.2 .35],"diff",cmap);
        
        set_legend_relative_position([.75 1.05 .25 .2])
        fix_axes(gcf,7,true)
        figure_path = strrep(figure_path,"crabspeed_diff","by_day");
        dataClass.save_figure_catch(ax,figure_path);
        close
        
        
        %%
        subplot(2,2,1)
        [ax,figure_path] = dataArray.plot_behaviour_timeseries(params,{"rolls"},.2,[25 50],[0 60],cmap);
        subplot(2,2,2)
        [ax,figure_path] = dataArray.plot_behaviour_timeseries(params,{"peran"},.2,[25 50],[0 100],cmap);
        subplot(2,2,3)
        [ax,figure_path] = dataArray.plot_behaviour_timeseries(params,{"turns"},.2,[25 50],[0 100],cmap);
        subplot(2,2,4)
        [ax,figure_path] = dataArray.plot_behaviour_timeseries(params,{"hunches"},.2,[25 50],[0 100],cmap);
        
        set_legend_relative_position([.75 1.05 .25 .2])
        fix_axes(gcf,7,true)        
        figure_path = strrep(figure_path,"hunches","by_day");
        dataClass.save_figure_catch(ax,figure_path);
        close
        
    end
end

%% plot large group

[n_prots,grp_prots] = dc.get_group_numbers('effector','protocol');
for ii = 1:n_prots
    idx = grp_prots == ii;
    
    sum_dur_roll = struct();
    ave_dur_roll = struct();
    ave_lat_roll = struct();
    ave_amp_roll = struct();
    ave_prp_roll = struct();
    ave_num_roll = struct();
    ave_amp_turn = struct();
    ave_dur_turn = struct();
    ave_lat_turn = struct();
    ave_area = struct();
    
    [n_exp,grp_exp] = dc(idx).get_group_numbers('driver','effector','protocol');
    for jj = 1:n_exp
        ix = find(idx)';
        f = grp_exp' == jj;
        ix = ix(f);
        temp = dc(ix).load_groups(false);
%         temp = dc(myfilt).load_groups(false);
        dA = temp.filter_by_size(params.area_cutoff);
        if isempty(dA.behaviour) | isempty(dA.timeseries) | ~isfield(dA.timeseries,'et')
            continue
        end
        
        p = dA.parse_protocol();
        frame = p.start + [0 p.length];
        
        sum_dur_roll = dA.append_2_plotstruct(sum_dur_roll,params,'rolls','duration','all','sum',frame);
        
        ave_dur_roll = dA.append_2_plotstruct(ave_dur_roll,params,'rolls','duration','all','mean',frame);
        
        ave_lat_roll = dA.append_2_plotstruct(ave_lat_roll,params,'rolls','start','first','mean',frame);
        
        ave_amp_roll = dA.append_2_plotstruct(ave_amp_roll,params,'rolls','amplitude','all','mean',frame);
        
        ave_prp_roll = dA.append_2_plotstruct(ave_prp_roll,params,'rolls','duration','all','proportion',frame);
        
        ave_num_roll = dA.append_2_plotstruct(ave_num_roll,params,'rolls','duration','all','count',frame);
        
        ave_amp_turn = dA.append_2_plotstruct(ave_amp_turn,params,'turns','amplitude','all','mean',frame);
        
        ave_dur_turn = dA.append_2_plotstruct(ave_dur_turn,params,'turns','duration','all','mean',frame);
        
        ave_lat_turn = dA.append_2_plotstruct(ave_lat_turn,params,'turns','duration','all','mean',frame);
        
        ave_area = dA.append_2_plotstruct(ave_area,params,'area','duration','all','mean',frame);
        
        if numel(ave_lat_turn.timestamp) ~= numel(ave_lat_turn.dates)
            break
        end
        
        clear temp dA
    end
   
    %%
    highlight_genos = {'GMR_SS01816','GMR_SS04232','GMR_SS04248'};
    
    %% sum_dur_roll
    subplot(2,1,1)
    [ax,figure_path] = dataClass.plot_grouped_beeswarm(sum_dur_roll,'descend',false,{'attp2'},highlight_genos,"Total Roll Duration (s)");
    subplot(2,1,2)
    [ax,figure_path] = dataClass.plot_grouped_beeswarm(sum_dur_roll,'descend',true,{'attp2'},highlight_genos,"\Delta Total Roll Duration (s)");
    fix_axes(gcf,9,true)
%     
    ax = findobj(gcf,'Type','Axes');
    arrayfun(@(x) set(x.XAxis,'FontSize',2),ax);
    arrayfun(@(x) set(x,'TickLength',[0,0]),ax);
    
    figure_path = dc(find(idx,1)).figure_directory;
    figure_path = fullfile(figure_path,"grouped",'sum_dur_roll');
    
    dataClass.save_figure_catch(ax,figure_path);
    close
%     
    %% ave_dur_roll
    subplot(2,1,1)
    [ax,figure_path] = dataClass.plot_grouped_beeswarm(ave_dur_roll,'descend',false,{'attp2'},highlight_genos,"Average Roll Duration (s)");
    subplot(2,1,2)
    [ax,figure_path] = dataClass.plot_grouped_beeswarm(ave_dur_roll,'descend',true,{'attp2'},highlight_genos,"\Delta Average Roll Duration (s)");
    fix_axes(gcf,9,true)
    
    ax = findobj(gcf,'Type','Axes');
    arrayfun(@(x) set(x.XAxis,'FontSize',2),ax);
    arrayfun(@(x) set(x,'TickLength',[0,0]),ax);

    figure_path = dc(find(idx,1)).figure_directory;
    figure_path = fullfile(figure_path,"grouped",'ave_dur_roll');

    dataClass.save_figure_catch(ax,figure_path);
    close
    
    %% ave_lat_roll
    subplot(2,1,1)
    [ax,figure_path] = dataClass.plot_grouped_beeswarm(ave_lat_roll,'descend',false,{'attp2'},highlight_genos,"Average Latency to Roll (s)");
    subplot(2,1,2)
    [ax,figure_path] = dataClass.plot_grouped_beeswarm(ave_lat_roll,'descend',true,{'attp2'},highlight_genos,"\Delta Average Latency to Roll (s)");
    fix_axes(gcf,9,true)
    
    ax = findobj(gcf,'Type','Axes');
    arrayfun(@(x) set(x.XAxis,'FontSize',2),ax);
    arrayfun(@(x) set(x,'TickLength',[0,0]),ax);

    figure_path = dc(find(idx,1)).figure_directory;
    figure_path = fullfile(figure_path,"grouped",'ave_lat_roll');

    dataClass.save_figure_catch(ax,figure_path);
    close
    
    %% ave_amp_roll
    subplot(2,1,1)
    [ax,figure_path] = dataClass.plot_grouped_beeswarm(ave_amp_roll,'descend',false,{'attp2'},highlight_genos,"Average Roll Amplitude");
    subplot(2,1,2)
    [ax,figure_path] = dataClass.plot_grouped_beeswarm(ave_amp_roll,'descend',true,{'attp2'},highlight_genos,"\Delta Average Roll Amplitude");
    fix_axes(gcf,9,true)
    
    ax = findobj(gcf,'Type','Axes');
    arrayfun(@(x) set(x.XAxis,'FontSize',2),ax);
    arrayfun(@(x) set(x,'TickLength',[0,0]),ax);
    
    figure_path = dc(find(idx,1)).figure_directory;
    figure_path = fullfile(figure_path,"grouped",'ave_amp_roll');

    dataClass.save_figure_catch(ax,figure_path);
    close
    
    %% ave_prp_roll
    subplot(2,1,1)
    [ax,figure_path] = dataClass.plot_grouped_beeswarm(ave_prp_roll,'descend',false,{'attp2'},highlight_genos,"Roll Proportion");
    subplot(2,1,2)
    [ax,figure_path] = dataClass.plot_grouped_beeswarm(ave_prp_roll,'descend',true,{'attp2'},highlight_genos,"\Delta Roll Proportion");
    fix_axes(gcf,9,true)
    
    ax = findobj(gcf,'Type','Axes');
    arrayfun(@(x) set(x.XAxis,'FontSize',2),ax);
    arrayfun(@(x) set(x,'TickLength',[0,0]),ax);

    figure_path = dc(find(idx,1)).figure_directory;
    figure_path = fullfile(figure_path,"grouped",'ave_prp_roll');

    dataClass.save_figure_catch(ax,figure_path);
    close
    
    %% ave_num_roll
    subplot(2,1,1)
    [ax,figure_path] = dataClass.plot_grouped_beeswarm(ave_num_roll,'descend',false,{'attp2'},highlight_genos,"Average Number of Rolls");
    fix_axes(gcf,9,true)
    subplot(2,1,2)
    [ax,figure_path] = dataClass.plot_grouped_beeswarm(ave_num_roll,'descend',true,{'attp2'},highlight_genos,"\Delta Average Number of Rolls");
    fix_axes(gcf,9,true)
    
    ax = findobj(gcf,'Type','Axes');
    arrayfun(@(x) set(x.XAxis,'FontSize',2),ax);
    arrayfun(@(x) set(x,'TickLength',[0,0]),ax);
    
    figure_path = dc(find(idx,1)).figure_directory;
    figure_path = fullfile(figure_path,"grouped",'ave_num_roll');

    dataClass.save_figure_catch(ax,figure_path);
    close
    
    %% ave_amp_turn
    subplot(2,1,1)
    [ax,figure_path] = dataClass.plot_grouped_beeswarm(ave_amp_turn,'descend',false,{'attp2'},highlight_genos,"Average Turn Amplitude (s)");
    subplot(2,1,2)
    [ax,figure_path] = dataClass.plot_grouped_beeswarm(ave_amp_turn,'descend',true,{'attp2'},highlight_genos,"\Delta Average Turn Amplitude (s)");
    fix_axes(gcf,9,true)
    
    ax = findobj(gcf,'Type','Axes');
    arrayfun(@(x) set(x.XAxis,'FontSize',2),ax);
    arrayfun(@(x) set(x,'TickLength',[0,0]),ax);

    figure_path = dc(find(idx,1)).figure_directory;
    figure_path = fullfile(figure_path,"grouped",'ave_amp_turn');

    dataClass.save_figure_catch(ax,figure_path);
    close
    
    %% ave_dur_turn
    subplot(2,1,1)
    [ax,figure_path] = dataClass.plot_grouped_beeswarm(ave_dur_turn,'descend',false,{'attp2'},highlight_genos,"Average Turn Duration (s)");
    subplot(2,1,2)
    [ax,figure_path] = dataClass.plot_grouped_beeswarm(ave_dur_turn,'descend',true,{'attp2'},highlight_genos,"\Delta Average Turn Duration (s)");
    fix_axes(gcf,9,true)
    
    ax = findobj(gcf,'Type','Axes');
    arrayfun(@(x) set(x.XAxis,'FontSize',2),ax);
    arrayfun(@(x) set(x,'TickLength',[0,0]),ax);

    figure_path = dc(find(idx,1)).figure_directory;
    figure_path = fullfile(figure_path,"grouped",'ave_dur_turn');

    dataClass.save_figure_catch(ax,figure_path);
    close
    
    %% ave_lat_turn
    subplot(2,1,1)
    [ax,figure_path] = dataClass.plot_grouped_beeswarm(ave_lat_turn,'descend',false,{'attp2'},highlight_genos,"Average Latency to Turn (s)");
    subplot(2,1,2)
    [ax,figure_path] = dataClass.plot_grouped_beeswarm(ave_lat_turn,'descend',true,{'attp2'},highlight_genos,"\Delta Average Latency to Turn (s)");
    fix_axes(gcf,9,true)
    
    ax = findobj(gcf,'Type','Axes');
    arrayfun(@(x) set(x.XAxis,'FontSize',2),ax);
    arrayfun(@(x) set(x,'TickLength',[0,0]),ax);

    figure_path = dc(find(idx,1)).figure_directory;
    figure_path = fullfile(figure_path,"grouped",'ave_lat_turn');

    dataClass.save_figure_catch(ax,figure_path);
    close
    
    %% ave_lat_turn
    subplot(2,1,1)
    [ax,figure_path] = dataClass.plot_grouped_beeswarm(ave_area,'descend',false,{'attp2'},highlight_genos,"Average Area (mm^{2})");
    subplot(2,1,2)
    [ax,figure_path] = dataClass.plot_grouped_beeswarm(ave_area,'descend',true,{'attp2'},highlight_genos,"\Delta Average Area (mm^{2})");
    fix_axes(gcf,9,true)
    
    ax = findobj(gcf,'Type','Axes');
    arrayfun(@(x) set(x.XAxis,'FontSize',2),ax);
    arrayfun(@(x) set(x,'TickLength',[0,0]),ax);
    
    figure_path = dc(find(idx,1)).figure_directory;
    figure_path = fullfile(figure_path,"grouped",'ave_area');

    dataClass.save_figure_catch(ax,figure_path);
    close
end



%%
fig = figure('Visible','off');
figpos = get(fig,'Position');
close(fig);
aspect = [1,2,1];
aspect(1:2) = aspect(1:2)./min(aspect(1:2));
figarea = figpos(3)*figpos(4);
x = sqrt(figarea./[aspect(2)/aspect(1)]);

figpos(3:4) = ceil(aspect(1:2)*x);

[n_groups,group_num] = dc.get_group_numbers('driver','effector','protocol');

for ii = 1:n_groups
    idx = group_num == ii;
    dct = dc(idx);
    %
    dA = dct.load_groups(true);
    dA = dA.filter_by_size(params.area_cutoff);
    
    features = {'curve','crabspeed'};
    for jj = 1:numel(features)
        feature = features{jj};
        fr1 = 25;
        fr2 = 40;

        switch feature
            case 'curve'
                cax = [0 60];
            case 'crabspeed'
                cax = [0 3];
            case 'speed'
                cax = [0 5];
            case 'midline'
                cax = [0 8];
        end

        f = cellfun(@(x) min(x)<fr1 & fr2<max(x), {dA.timeseries.et});
        dA.timeseries = dA.timeseries(f);

        counts = arrayfun(@(x) numel(x.et), dA.timeseries);
        ids = repelem([dA.timeseries.uniID], 1,counts);
        et = [dA.timeseries.et];
        val = [dA.timeseries.(feature)];

        bins = [fr1:.05:fr2];
        Y = discretize(et,bins);
        f = ~isnan(Y);
        [~,~,idf] = unique(ids(f),'stable');

        mat = accumarray([idf,Y(f)'],val(f)',[],@mean);
        
        fig = figure('Position',figpos);
        cmap = colormap(viridis);
        im = imagesc([fr1 fr2],[],mat);
        caxis(cax)
        cb = colorbar;
        cb.Label.String = feature;
        pbaspect([1,2,1])

        ax = set(gca,'YTick',[]);

        dA.plot_stimulation_time(params);
        bx = findobj(gca,'Tag','stimbox');
        ymax = max(get(gca,'YLim'));
        bx.YData = ymax-bx.YData;

        chi = get(gca,'Children');
        chi = [findobj(chi,'Tag','stimbox'); findobj(chi,'-not','Tag','stimbox')];
        set(gca,'Children',chi);

        fix_axes(gcf,11,false)
        
        % Save figure
        figure_dir = dA.figure_directory;
        fig_type = strcat("heatmap_",feature);
        fname = strcat(dA.get_full_genotype);

        fig_path = fullfile(figure_dir,fig_type,fname);
        if ~isdir(fileparts(fig_path))
            mkdir(fileparts(fig_path));
        end

        dataClass.save_figure_catch(gcf,fig_path);
        close
    end
end



%%

dc = dataClass(filepaths(contains(filepaths,'choreography'))');
dc = dc.update_figure_directories(params);

[prot_n,prot_grp] = dc.get_group_numbers('protocol');

for ii = 1:prot_n
    idx = prot_grp == ii;
    dct = dc(idx);
    
    [stamp_n,stamp_grp] = dc(idx).get_group_numbers('timestamp');
    for jj = 89:stamp_n
        ix = stamp_grp == jj;
        dA = dct(ix).load_groups(true);
        
        counts = arrayfun(@(x) numel(x.et), dA.timeseries, 'ErrorHandler', @(a,b) nan);

        et = [dA.timeseries.et];
        ids = [dA.timeseries.uniID];
        ids = repelem(ids,1,counts(counts>0));

        bins = [0:.2:ceil(max(et))];

        metric = [dA.timeseries.crabspeed];
        Y = discretize(et,bins);
        y = accumarray(Y',metric',[numel(bins) 1],@mean);
        n = accumarray(Y',ids',[numel(bins) 1],@(x) numel(unique(x)));
        dev = accumarray(Y',metric',[numel(bins) 1],@std);

        err = dev./sqrt(n);
        err(isnan(err)) = 0;

        err_y = [y+err;flipud(y-err)]';
        err_x = [bins(1:end),fliplr(bins(1:end))];

        titl = dA.get_full_timestamp;
        
        if ~ishandle(ax)        
            title(titl,'Interpreter','none')
            p = plot(bins(1:end),y,'Color',[0 0 0], 'LineWidth', 2);
            ylim([0 1.5])
            pbaspect([2,1,1])
            dA.plot_stimulation_time(params);
            ax = gca;
        else
            set(p,'XData',bins(1:end),'YData',y)
            ob = findobj(gca,'Type','Patch');
            hold on
            dA.plot_stimulation_time(params);
            hold off
        end
        
        res = input('enter value: ','s');
        if isempty(res)
            
        elseif strcmp(res,'n')
            tstamp = dA.get_full_timestamp;
            update_blacklist({char(tstamp)})
        end

    end 
end
