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

filterByDriver = true;

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

%% individual plots

[n_groups,group_num] = dc.get_group_numbers('driver','effector','protocol');

for ii = 1:n_groups
    idx = group_num == ii;
    dct = dc(idx);
    
    prot = dct.parse_protocol();
    xlm = prot(1).start + [-5,prot(1).length+5];
    
    %
    dA = dct.load_groups(true);
    dA = dA.filter_by_size(params.area_cutoff);
    %
%     [ax,figure_path] = dA.plot_ridgeline_timeseries(params,'crabspeed',xlm(1),xlm(2),5/12);
%     fix_axes(gcf,14,false)
%     dataClass.save_figure_catch(ax,figure_path);
%     close
%     %
%     [ax,figure_path] = dA.plot_ridgeline_timeseries(params,'curve',xlm(1),xlm(2),5/100);
%     fix_axes(gcf,14,false)
%     dataClass.save_figure_catch(ax,figure_path);
%     close
%     %
%     [ax,figure_path] = dA.plot_ridgeline_timeseries(params,'speed',xlm(1),xlm(2),5/12);
%     fix_axes(gcf,14,false)
%     dataClass.save_figure_catch(ax,figure_path);
%     close
%     %
%     [ax,figure_path] = dA.plot_ridgeline_timeseries(params,'midline',xlm(1),xlm(2),5/5);
%     fix_axes(gcf,14,false)
%     dataClass.save_figure_catch(ax,figure_path);
%     close
    
    %
    behaviours = {"peran","rolls"};
    [ax,figure_path] = dA.plot_ethogram(params,behaviours,xlm(1),xlm(2),0);
    fix_axes(gcf,14,false)
    dataClass.save_figure_catch(ax,figure_path);
    close
%     %
    behaviours = {"peran","turns"};
    [ax,figure_path] = dA.plot_ethogram(params,behaviours,xlm(1),xlm(2),0);
    fix_axes(gcf,14,false)
    dataClass.save_figure_catch(ax,figure_path);
    close
    
    behaviours = {"run_tss","roll_tss"};
    [ax,figure_path] = dA.plot_ethogram(params,behaviours,xlm(1),xlm(2),0);
    fix_axes(gcf,14,false)
    dataClass.save_figure_catch(ax,figure_path);
    close
    
    behaviours = {"run_tss","stop_tss"};
    [ax,figure_path] = dA.plot_ethogram(params,behaviours,xlm(1),xlm(2),0);
    fix_axes(gcf,14,false)
    dataClass.save_figure_catch(ax,figure_path);
    close
    
    behaviours = {"run_tss","cast_tss"};
    [ax,figure_path] = dA.plot_ethogram(params,behaviours,xlm(1),xlm(2),0);
    fix_axes(gcf,14,false)
    dataClass.save_figure_catch(ax,figure_path);
    close
    
    behaviours = {"run_tss","hunch_tss"};
    [ax,figure_path] = dA.plot_ethogram(params,behaviours,xlm(1),xlm(2),0);
    fix_axes(gcf,14,false)
    dataClass.save_figure_catch(ax,figure_path);
    close
    
    behaviours = {"run_tss","back_tss"};
    [ax,figure_path] = dA.plot_ethogram(params,behaviours,xlm(1),xlm(2),0);
    fix_axes(gcf,14,false)
    dataClass.save_figure_catch(ax,figure_path);
    close
    
    
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
        
        p = dataArray.parse_protocol();
        xlm = p(1).start + [-5 20];

        %%
%         figure('Resize','off')
%         fig = figure('Visible','off');
%         subplot(2,2,1)
        [ax,figure_path] = dataArray.plot_timeseries(params,"curve",xlm(1),xlm(2),[0 50],"norm",cmap);
        fix_axes(gcf,14,true)
        set_legend_relative_position('topleft',[1 1 0 0])
        dataClass.save_figure_catch(ax,figure_path);
        close
%         subplot(2,2,3)
        [ax,figure_path] = dataArray.plot_timeseries(params,"crabspeed",xlm(1),xlm(2),[0 2.5],"norm",cmap);
        fix_axes(gcf,14,true)
        set_legend_relative_position('topleft',[1 1 0 0])
        dataClass.save_figure_catch(ax,figure_path);
        close
%         subplot(2,2,2)
        [ax,figure_path] = dataArray.plot_timeseries(params,"curve",xlm(1),xlm(2),[-1 4],"diff",cmap);
        fix_axes(gcf,14,true)
        set_legend_relative_position('topleft',[1 1 0 0])
        dataClass.save_figure_catch(ax,figure_path);
        close
%         subplot(2,2,4)
        [ax,figure_path] = dataArray.plot_timeseries(params,"crabspeed",xlm(1),xlm(2),[-.2 .35],"diff",cmap);
        fix_axes(gcf,14,true)
        set_legend_relative_position('topleft',[1 1 0 0])
        dataClass.save_figure_catch(ax,figure_path);
        close
%         subplot(2,2,2)
        [ax,figure_path] = dataArray.plot_timeseries(params,"speeddiff",xlm(1),xlm(2)+15,[0 2],"norm",cmap);
        fix_axes(gcf,14,true)
        set_legend_relative_position('topleft',[1 1 0 0])
        dataClass.save_figure_catch(ax,figure_path);
        close
%         subplot(2,2,4)
        [ax,figure_path] = dataArray.plot_timeseries(params,"speeddiff",xlm(1),xlm(2)+15,[-.1 .2],"diff",cmap);
        fix_axes(gcf,14,true)
        set_legend_relative_position('topleft',[1 1 0 0])
        dataClass.save_figure_catch(ax,figure_path);
        close
                
        %%
%         fig = figure('Visible','off');
%         subplot(2,2,1)
        [ax,figure_path] = dataArray.plot_behaviour_timeseries(params,{"rolls"},.2,xlm,[0 60],cmap);
        fix_axes(gcf,14,true)
        set_legend_relative_position('topleft',[1 1 0 0])
        dataClass.save_figure_catch(ax,figure_path);
        close
        
        [ax,figure_path] = dataArray.plot_behaviour_timeseries(params,{"peran"},.2,xlm,[0 100],cmap);
        fix_axes(gcf,14,true)
        set_legend_relative_position('topleft',[1 1 0 0])
        dataClass.save_figure_catch(ax,figure_path);
        close
        
        [ax,figure_path] = dataArray.plot_behaviour_timeseries(params,{"turns"},.2,xlm,[0 100],cmap);
        fix_axes(gcf,14,true)
        set_legend_relative_position('topleft',[1 1 0 0])
        dataClass.save_figure_catch(ax,figure_path);
        close
        
        [ax,figure_path] = dataArray.plot_behaviour_timeseries(params,{"hunches"},.2,xlm,[0 100],cmap);
        fix_axes(gcf,14,true)
        set_legend_relative_position('topleft',[1 1 0 0])
        dataClass.save_figure_catch(ax,figure_path);
        close
        
    end
end

%% plot large group
% 
% make_it_tight = true;
% subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.05], [0.1 0.01], [0.1 0.01]);
% if ~make_it_tight,  clear subplot;  end

cmap = cbrewer('qual','Dark2',8);
cmap = cbrewer('qual','Set2',8);

highlightStruct = struct('attp2',[1 1 2.55]./2.55,...
            'GMR_SS01816',cmap(2,:),...
            'GMR_SS01817',cmap(2,:),...
            'GMR_SS01816n',cmap(2,:),...
            'GMR_SS01817n',cmap(2,:),...
            'GMR_SS01750',cmap(1,:),...
            'GMR_SS04232',cmap(1,:),...
            'GMR_SS00666',cmap(1,:),...
            'GMR_SS04052',cmap(1,:),...
            'GMR_SS01951',cmap(5,:),...
            'GMR_SS04245',cmap(5,:),...
            'GMR_SS04189',cmap(5,:),...
            'GMR_SS04248',cmap(4,:));

% highlightStruct = struct('backup',[0.3 0.3 0.3],...
%             'GMR_SS01750',[207 177 63]/255,...
%             'GMR_SS04232',[231 156 88]/255,...
%             'GMR_SS01816',[96 186 164]/255,...
%             'attp2',[114 136 190]/255,...
%             'GMR_SS01817',[243 214 99]/255,...
%             'GMR_SS02007',[226 127 175]/255,...
%             'GMR_SS04052',[225 114 98]/255,...
%             'GMR_SS04248',[191 214 231]/255);
        
[n_prots,grp_prots] = dc.get_group_numbers('effector','protocol');
for ii = 1:n_prots
    idx = grp_prots == ii;
    
    plotStruct = struct();
    
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
        frame = p(1).start + [0 p(1).length];
        
        plotStruct = dA.append_2_plotstruct(plotStruct,params,'rolls','duration','all','sum',frame);
        plotStruct = dA.append_2_plotstruct(plotStruct,params,'rolls','duration','all','mean',frame);
        plotStruct = dA.append_2_plotstruct(plotStruct,params,'rolls','start','first','mean',frame);
        plotStruct = dA.append_2_plotstruct(plotStruct,params,'rolls','amplitude','all','mean',frame);
        plotStruct = dA.append_2_plotstruct(plotStruct,params,'rolls','duration','all','proportion',frame);
        plotStruct = dA.append_2_plotstruct(plotStruct,params,'rolls','duration','all','count',frame);
        plotStruct = dA.append_2_plotstruct(plotStruct,params,'turns','amplitude','all','mean',frame);
        plotStruct = dA.append_2_plotstruct(plotStruct,params,'turns','duration','all','mean',frame);
        plotStruct = dA.append_2_plotstruct(plotStruct,params,'turns','start','first','mean',frame);
        plotStruct = dA.append_2_plotstruct(plotStruct,params,'area','duration','all','mean',frame);
        
%         if numel(ave_lat_turn.timestamp) ~= numel(ave_lat_turn.dates)
%             break
%         end
        
        clear temp dA
    end
    
%     if numel(fieldnames(ave_area)) == 0
%         continue
%     end
   
    %%    
    figure_path = fullfile(dc(find(idx,1)).figure_directory,...
        'grouped');
    if ~isdir(figure_path)
        mkdir(figure_path)
    end
    
    %% Set up loop
    names = {'sum_dur_roll','ave_dur_roll',...
                'ave_lat_roll','ave_amp_roll',...
                'ave_prp_roll','ave_num_roll',...
                'ave_amp_turn','ave_dur_turn',...
                'ave_lat_turn','ave_area'};
            
    xlabels = {"Total Roll Duration (s)","Average Roll Duration (s)",...
                "Average Latency to Roll (s)","Average Roll Amplitude",...
                "Roll Proportion","Average Number of Rolls",...
                "Average Turn Amplitude (s)","Average Turn Duration (s)",...
                "Average Latency to Turn (s)","Average Area (mm^{2})"};

    %% sum_dur_roll
    for jj = 1:numel(plotStruct)
        figName = fullfile(figure_path,strcat(names{jj}));
        [ax,~,stats] = dataClass.plot_violins(plotStruct(jj),'descend',false,{'attp2'},highlightStruct,xlabels{jj});
        fix_axes(gcf,14,true)
        set(ax.XAxis,'FontSize',4)
        set(ax,'TickLength',[0,0])
        dataClass.save_figure_catch(ax,figName);
        close
        save_stats_struct(figName,stats);

        figName = fullfile(figure_path,strcat(names{jj},'_norm'));
        [ax,~,stats] = dataClass.plot_violins(plotStruct(jj),'descend',true,{'attp2'},highlightStruct,strcat('\Delta ',xlabels{jj}));
        fix_axes(gcf,14,true)
        set(ax.XAxis,'FontSize',4)
        set(ax,'TickLength',[0,0])
        dataClass.save_figure_catch(ax,figName);
        close
        save_stats_struct(figName,stats);
    end
    
    % Get n of genotype and n of control by same day
%     groups = plotStruct(1).group;
%     dates = plotStruct(1).dates;
%     
%     f = groups == 1;
%     controlDates = [dates(f)]';
%     expGroups = groups(~f);
%     expDates = dates(~f);
%     
%     nControl = accumarray(expGroups',expDates',[],@(x) sum(ismember(controlDates,x)));
%     nExp = accumarray(expGroups',expDates',[],@(x) sum(ismember(x,controlDates)));
%     
%     printCell = [cellstr(plotStruct(1).labels)', num2cell([nExp,nControl])];
%     [~,I] = sort(plotStruct(1).labels);
%     printCell = printCell(I,:)';
%     
%     fid = fopen('./screen_n_values.txt','w');
%     fprintf(fid,'%s %d %d \n',printCell{:});
%     fclose(fid);
    
end

%%
% highlight_genos = {'GMR_SS01816','GMR_SS04232','GMR_SS04248'};

cmap = cbrewer('qual','Dark2',8);
cmap = cbrewer('qual','Set2',8);

highlightStruct = struct('attp2',[1 1 2.55]./2.55,...
            'GMR_SS01816',cmap(2,:),...
            'GMR_SS01817',cmap(2,:),...
            'GMR_SS01816n',cmap(2,:),...
            'GMR_SS01817n',cmap(2,:),...
            'GMR_SS01750',cmap(1,:),...
            'GMR_SS04232',cmap(1,:),...
            'GMR_SS00666',cmap(1,:),...
            'GMR_SS04052',cmap(1,:),...
            'GMR_SS01951',cmap(5,:),...
            'GMR_SS04245',cmap(5,:),...
            'GMR_SS04189',cmap(5,:),...
            'GMR_SS04248',cmap(4,:));
        
% highlightStruct = struct('backup',[0.3 0.3 0.3],...
%             'GMR_SS01750',[207 177 63]/255,...
%             'GMR_SS04232',[231 156 88]/255,...
%             'GMR_SS01816',[96 186 164]/255,...
%             'attp2',[114 136 190]/255,...
%             'GMR_SS01817',[243 214 99]/255,...
%             'GMR_SS02007',[226 127 175]/255,...
%             'GMR_SS04052',[225 114 98]/255,...
%             'GMR_SS04248',[191 214 231]/255);

[n_prots,grp_prots] = dc.get_group_numbers('effector','protocol');
for ii = 1:n_prots
    plotStruct = struct();
    
    idx = grp_prots == ii;
        
    figure_path = fullfile(dc(find(idx,1)).figure_directory,...
        'grouped');
    if ~isdir(figure_path)
        mkdir(figure_path)
    end
    
    [n_exp,grp_exp] = dc(idx).get_group_numbers('driver','effector','protocol');
    for jj = 1:n_exp
        ix = find(idx)';
        f = grp_exp' == jj;
        ix = ix(f);
        temp = dc(ix).load_groups(false);
%         temp = dc(myfilt).load_groups(false);
        dA = temp.filter_by_size(params.area_cutoff);
        if isempty(dA.timeseries) | ~isfield(dA.timeseries,'et')
            continue
        end
        
        p = dA.parse_protocol();
        frame = p(1).start + [0 p(1).length];
        
        plotStruct = dA.append_choreography_metric(plotStruct,'speed','mean',[15 30],'all');
        plotStruct = dA.append_choreography_metric(plotStruct,'curve','mean',[15 30],'all');
        plotStruct = dA.append_choreography_metric(plotStruct,'crabspeed','mean',[15 30],'all');
        
        plotStruct = dA.append_choreography_metric(plotStruct,'curve','mean',frame,'all');
        plotStruct = dA.append_choreography_metric(plotStruct,'crabspeed','mean',frame,'all');  
    end
    
    for jj = 1:numel(plotStruct)
        figName = fullfile(figure_path,plotStruct(jj).name);
        label = strsplit(plotStruct(jj).name,'_');
        [ax,~,stats] = dataClass.plot_violins(plotStruct(jj),'descend',false,{'attp2'},highlightStruct, strcat(label(1),'_',label(3)) );
        fix_axes(gcf,14,true)
        set(ax.XAxis,'FontSize',4)
        set(ax,'TickLength',[0,0])
        dataClass.save_figure_catch(ax,figName);
        close
        save_stats_struct(figName,stats);

        figName = fullfile(figure_path,strcat(plotStruct(jj).name,'_norm'));
        [ax,~,stats] = dataClass.plot_violins(plotStruct(jj),'descend',true,{'attp2'},highlightStruct, strcat('\Delta ',label(1),'_',label(3)) );
        fix_axes(gcf,14,true)
        set(ax.XAxis,'FontSize',4)
        set(ax,'TickLength',[0,0])
        dataClass.save_figure_catch(ax,figName);
        close
        save_stats_struct(figName,stats);
    end
end

%%
% fig = figure('Visible','off');
% figpos = get(fig,'Position');
% close(fig);
% aspect = [1,2,1];
% aspect(1:2) = aspect(1:2)./min(aspect(1:2));
% figarea = figpos(3)*figpos(4);
% x = sqrt(figarea./[aspect(2)/aspect(1)]);
% 
% figpos(3:4) = ceil(aspect(1:2)*x);

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
        
        p = dA.parse_protocol();
        frame = p(1).start + [0 p(1).length];
        fr1 = frame(1)-5;
        fr2 = fr1+15;

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
        
        if isempty(fieldnames(dA.timeseries))
            continue
        end

        counts = arrayfun(@(x) numel(x.et), dA.timeseries);
        ids = repelem([dA.timeseries.uniID], 1,counts);
        et = [dA.timeseries.et];
        val = [dA.timeseries.(feature)];

        bins = [fr1:.05:fr2];
        Y = discretize(et,bins);
        f = ~isnan(Y);
        [~,~,idf] = unique(ids(f),'stable');

        mat = accumarray([idf,Y(f)'],val(f)',[],@mean);
        
%         fig = figure('Position',figpos);
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
        
        title(dA.driver,'Interpreter','none')

        fix_axes(gcf,14,false)
        
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

[n_prots,grp_prots] = dc.get_group_numbers('effector','protocol');
for ii = 1:n_prots
    idx = grp_prots == ii;
    
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
        frame = p(1).start + [0 p(1).length];
        
        
        Behaviour = 'peran';
        Metric = 'frequency';
        Instance = 'first';
        Method = 'mean';
        obj = dA;
        
        f = strcmp({obj.behaviour.behaviour},Behaviour);
        
        pS = obj.behaviour(f);
        [~,id_filt] = ismember(pS.uniID,obj.uniID);
        pS.timestamp = obj.timestamp_index(id_filt);
        
        filt = pS.track_start < frame(1) & pS.track_end > frame(2);
        
        if any(strcmp(Instance,{'first','last'}))
            switch Instance
                case 'first'
                    func = @(x) min([x(x>frame(1) & x<frame(2));nan]);
                case 'last'
                    func = @(x) max([x(x>frame(1) & x<p.frame(2));nan]);
            end
            ref = [pS.uniID;pS.start]';
            [dds,~,uds] = unique(pS.uniID,'stable');
            val = accumarray(uds,pS.start',[],func);
            nf = ismember(ref,[dds',val],'rows');
            filt = filt & nf';
        end
        
        structFilter = structfun(@(x) numel(x)==numel(filt), pS);
        pS = structfun(@(x) x(filt), pS,'UniformOutput', false, 'ErrorHandler',@(x,y) NaN);
        pS.behaviour = obj.behaviour(f).behaviour;
        
        [~,~,id] = unique(pS.uniID,'stable');
        pS.uniID = id';
        ts = pS.timestamp;
        
        % Get metric
        switch Metric
            case 'duration'
                val = pS.end - pS.start;
            case 'start'
                val = pS.start;
            case 'end'
                val = pS.end;
            case 'frequency'
                val = pS.frequency;
            case 'amplitude'
                val = pS.amplitude;
        end
        
        % metric by animal
        tsIndex = accumarray(id,ts',[],@mean);
        switch Method
            case 'all'
                tot = val';
            case 'mean'
                tot = accumarray(id,val',[],@(x) mean(x,'omitnan'));
            case 'median'
                tot = accumarray(id,val',[],@(x) median(x,'omitnan'));
            case 'sum'
                tot = accumarray(id,val',[],@(x) sum(x,'omitnan'));
                tot(tot==0) = nan;
            case 'count'
                tot = accumarray(id,val',[],@(x) sum(~isnan(x)));
                tot(tot==0) = nan;
            case 'proportion'
                tot = accumarray(id,val',[],@(x) any(~isnan(x)));
        end
        
    end
end

%%

dc = dataClass(filepaths(contains(filepaths,'choreography'))');
dc = dc.update_figure_directories(params);

clear ax
[prot_n,prot_grp] = dc.get_group_numbers('protocol');
for ii = 1:prot_n
    idx = prot_grp == ii;
    dct = dc(idx);
    
    [stamp_n,stamp_grp] = dc(idx).get_group_numbers('timestamp');
    for jj = 1:stamp_n
        ix = stamp_grp == jj;
        try
            dA = dct(ix).load_groups(true);
        catch
            tstamp = dct(ix).get_full_timestamp;
            update_blacklist({char(tstamp)})
            continue
        end
            
        counts = arrayfun(@(x) numel(x.et), dA.timeseries, 'ErrorHandler', @(a,b) nan);

        et = [dA.timeseries.et];
        ids = [dA.timeseries.uniID];
        ids = repelem(ids,1,counts(counts>0));

        bins = [0:.2:ceil(max(et))];

        metric = [dA.timeseries.curve];
        Y = discretize(et,bins);
        y = accumarray(Y',metric',[numel(bins) 1],@mean);
        n = accumarray(Y',ids',[numel(bins) 1],@(x) numel(unique(x)));
        dev = accumarray(Y',metric',[numel(bins) 1],@std);

        err = dev./sqrt(n);
        err(isnan(err)) = 0;

        err_y = [y+err;flipud(y-err)]';
        err_x = [bins(1:end),fliplr(bins(1:end))];

        titl = dA.get_full_timestamp;
        
        fprintf(strcat(titl,'\n'))
        
        if ~exist('ax')     
            title(titl,'Interpreter','none')
            p = plot(bins(1:end),y,'Color',[0 0 0], 'LineWidth', 2);
            ylim([0 max(y*3)])
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
