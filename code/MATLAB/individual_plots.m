%% individual_plots

selectionDialog = true;

PlotGenos = {'attp2',...
        'GMR_SS01816',...
        'GMR_SS01817',...
        'GMR_SS01816n',...
        'GMR_SS01817n',...
        'GMR_SS01321',...
        'GMR_SS01792',...
        'GMR_SS04189',...
        'GMR_SS01750',...
        'GMR_SS00666',...
        'GMR_SS00869',...
        'GMR_SS04052',...
        'GMR_SS04232',...
        'GMR_SS04248'};
    
PlotGenos = {'attp2',...
        'GMR_SS01816',...
        'GMR_SS01817'}
    
    
PantonePastel = [[207 177 63];...
            [231 156 88];...
            [96 186 164];...
            [114 136 190];...
            [243 214 99];...
            [226 127 175];...
            [225 114 98];...
            [191 214 231] ] / 255;
        
colourRef = struct('attp2',[0.3 0.3 0.3],...
            'GMR_SS01750',[207 177 63]/255,...
            'GMR_SS04232',[231 156 88]/255,...
            'GMR_SS04189',[231 156 88]/255,...
            'GMR_SS04052',[96 186 164]/255,...
            'GMR_SS02175',[96 186 164]/255,...
            'GMR_SS00666',[114 136 190]/255,...
            'GMR_SS01792',[114 136 190]/255,...
            'GMR_SS01817',[243 214 99]/255,...
            'GMR_SS01817n',[243 214 99]/255,...
            'GMR_SS04248',[226 127 175]/255,...
            'GMR_SS01321',[226 127 175]/255,...
            'GMR_SS01816',[225 114 98]/255,...
            'GMR_SS01816n',[225 114 98]/255,...
            'tbd',[191 214 231]/255);
    
cmap = PantonePastel;

% cmap = cbrewer('qual','Dark2',8);
% 
% colourRef = struct('attp2',cmap(8,:),...
%             'GMR_SS01750',[207 177 63]/255,...
%             'GMR_SS04232',cmap(2,:),...
%             'GMR_SS01816',cmap(1,:),...
%             'GMR_SS04248',cmap(3,:),...
%             'GMR_SS01817',[243 214 99]/255,...
%             'GMR_SS00666',[226 127 175]/255,...
%             'GMR_SS04052',[225 114 98]/255,...
%             'tbd',[191 214 231]/255);

% 
% cmap = cbrewer('qual','Dark2',8);
% % cmap(1,:) = [0.2 0.2 0.2];
% 
% x = repelem([1:0.1:10],8,1);
% y = sin(x - 0.5*[0:7]');
% 
% Lines = line(x',y','LineWidth',3);
% ylim([-3 3])
% pbaspect([2,1,1])
% 
% cmapCell = num2cell(cmap,2);
% [Lines(:).Color] = cmapCell{:};
% fix_axes(gcf,14,true)
% 
% set(gcf,'InvertHardCopy','off')
% print(gcf,'colortest','-dsvg','-painters')

%% 

if selectionDialog
    [listDisp,~,listIdx] = unique([dc.driver]);
    chx = listdlg('ListString',listDisp);
    genoFilter = ismember(listIdx,chx);
else
    genoFilter = ismember([dc.driver],PlotGenos);
end
dcS = dc(genoFilter);

%% paired plots
[n_prot,grp_prot] = dcS.get_group_numbers('effector','protocol');
% [n_prot,grp_prot] = dcS.get_group_numbers('protocol');
for ii = 1:n_prot
    fprintf("%d/%d \n",ii,n_prot);
    datagrp = [];
    
    f1 = grp_prot == ii;
    dA = dcS(f1).load_groups(true);
    if numel(dA(1).effector) > 1
        effectors = cellfun(@(x) x(end), {dA.effector}, 'UniformOutput', false);
        [dA.effector] = effectors{:};
    end
    [~,I] = sort([dA.driver],'ascend');
    dA = dA(I);
    ff = [dA.driver]=='attp2';    
    dA = [dA(ff),dA(~ff)];
    dA = dA.filter_by_size(params.area_cutoff);

    prot = dA.parse_protocol();
    xlm = prot(1).start + [-5,20];

    [uniqueDrivers,~,idxDrivers] = unique([dA.driver],'stable');
    counts = accumarray(idxDrivers,1);
    repeats = arrayfun(@(x) [1:x], counts, 'UniformOutput', false);
    repeats = fliplr([repeats{:}]);
%     [~,I] = ismember(unique(uniqueDrivers,'stable'),fieldnames(colourRef));
    [~,I] = ismember([dA.driver],fieldnames(colourRef));
    cmapCell = struct2cell(colourRef);
    cmapFound = vertcat(cmapCell{I(I~=0)});
    cmapFound = cmapFound - [cmapFound./2].*[repeats(I~=0)'-1];
%     cmapFound = cmapFound + (0.2.*[1 1 1].*[repeats'-1]);
    cmapFound(cmapFound>1) = 1;
    cmapFound(cmapFound<0) = 0;
    cmapFilt = ~ismember(cmap,cmapFound,'rows');
    cmapOrdered = vertcat(cmapFound,cmap(cmapFilt,:));

    figPrefix = cellfun(@(x) x(end-3:end),cellstr(unique([dcS.driver])),'UniformOutput', false);
    figPrefix(1:end-1) = strcat(figPrefix(1:end-1),'-');
    figPrefix = strcat(figPrefix{:});
    figurePath = fullfile(dA(1).figure_directory,figPrefix);
    
    if numel(fieldnames(vertcat(dA.timeseries))) == 0 | ~numel(vertcat(dA.timeseries))
        ax = [];
        fig_path = [];
        continue
    end

    %%
    [ax,figure_path] = dA.plot_timeseries(params,"curve",xlm(1),xlm(2),[0 50],"norm",cmapOrdered);
    fix_axes(gcf,14,true);
    set_legend_relative_position('topleft',[1 1 0 0])
    dataClass.save_figure_catch(ax,strrep(figure_path,fileparts(figure_path),figurePath));
    close

    [ax,figure_path] = dA.plot_timeseries(params,"crabspeed",xlm(1),xlm(2),[0 2.5],"norm",cmapOrdered);
    fix_axes(gcf,14,true);
    set_legend_relative_position('topleft',[1 1 0 0])
    dataClass.save_figure_catch(ax,strrep(figure_path,fileparts(figure_path),figurePath));
    close

    [ax,figure_path] = dA.plot_timeseries(params,"curve",xlm(1),xlm(2),[-1 4],"diff",cmapOrdered);
    fix_axes(gcf,14,true);
    set_legend_relative_position('topleft',[1 1 0 0])
    dataClass.save_figure_catch(ax,strrep(figure_path,fileparts(figure_path),figurePath));
    close

    [ax,figure_path] = dA.plot_timeseries(params,"crabspeed",xlm(1),xlm(2),[-.2 .35],"diff",cmapOrdered);
    fix_axes(gcf,14,true);
    set_legend_relative_position('topleft',[1 1 0 0])
    dataClass.save_figure_catch(ax,strrep(figure_path,fileparts(figure_path),figurePath));
    close
    
    [ax,figure_path] = dA.plot_timeseries(params,"speeddiff",xlm(1),xlm(2)+15,[0 1.2],"norm",cmapOrdered);
    fix_axes(gcf,14,true);
    set_legend_relative_position('topleft',[1 1 0 0])
    dataClass.save_figure_catch(ax,strrep(figure_path,fileparts(figure_path),figurePath));
    close
    
    [ax,figure_path] = dA.plot_timeseries(params,"speeddiff",xlm(1),xlm(2)+15,[-.15 .2],"diff",cmapOrdered);
    fix_axes(gcf,14,true);
    set_legend_relative_position('topleft',[1 1 0 0])
    dataClass.save_figure_catch(ax,strrep(figure_path,fileparts(figure_path),figurePath));
    close


    %%
    [ax,figure_path] = dA.plot_behaviour_timeseries(params,{"rolls"},.2,xlm,[0 60],cmapOrdered);
    fix_axes(gcf,14,true)
    set_legend_relative_position('topleft',[1 1 0 0])
    dataClass.save_figure_catch(ax,strrep(figure_path,fileparts(figure_path),figurePath));
    close

    [ax,figure_path] = dA.plot_behaviour_timeseries(params,{"peran"},.2,xlm,[0 100],cmapOrdered);
    fix_axes(gcf,14,true)
    set_legend_relative_position('topleft',[1 1 0 0])
    dataClass.save_figure_catch(ax,strrep(figure_path,fileparts(figure_path),figurePath));
    close

    [ax,figure_path] = dA.plot_behaviour_timeseries(params,{"turns"},.2,xlm,[0 100],cmapOrdered);
    fix_axes(gcf,14,true)
    set_legend_relative_position('topleft',[1 1 0 0])
    dataClass.save_figure_catch(ax,strrep(figure_path,fileparts(figure_path),figurePath));
    close

    [ax,figure_path] = dA.plot_behaviour_timeseries(params,{"hunches"},.2,xlm,[0 100],cmapOrdered);
    fix_axes(gcf,14,true)
    set_legend_relative_position('topleft',[1 1 0 0])
    dataClass.save_figure_catch(ax,strrep(figure_path,fileparts(figure_path),figurePath));
    close
    
end


%% 
%%%%

cmap = cbrewer('qual','Dark2',8);

highlightStruct = struct('attp2',[1 1 2.55]./2.55,...
            'GMR_SS01816',cmap(1,:),...
            'GMR_SS01817',cmap(1,:),...
            'GMR_SS01750',cmap(2,:),...
            'GMR_SS04232',cmap(2,:),...
            'GMR_SS00666',cmap(2,:),...
            'GMR_SS04052',cmap(2,:),...
            'GMR_SS01951',cmap(4,:),...
            'GMR_SS04245',cmap(4,:),...
            'GMR_SS04189',cmap(4,:),...
            'GMR_SS04248',cmap(5,:));

[n_prots,grp_prots] = dcS.get_group_numbers('effector','protocol');
% [n_prots,grp_prots] = dcS.get_group_numbers('protocol');
for ii = 1:n_prots
    idx = grp_prots == ii;
    
    plotStruct = struct();
    
    [n_exp,grp_exp] = dcS(idx).get_group_numbers('driver','effector','protocol');
    for jj = 1:n_exp
        ix = find(idx)';
        f = grp_exp' == jj;
        ix = ix(f);
        temp = dcS(ix).load_groups(false);
%         temp = dc(myfilt).load_groups(false);
        dA = temp.filter_by_size(params.area_cutoff);
        if isempty(dA.behaviour) | isempty(dA.timeseries) | ~isfield(dA.timeseries,'et')
            continue
        end
        
        p = dA.parse_protocol();
        frame = p(1).start + [0 15];
        frameFirst5 = frame(1)+[0 5];
        frameBaseline = frame(1)+[-15 0];
        frame5After = frame(2) + [0 5];
        frame15After = frame(2) + [0 15];
        
        
        %%%% Baseline metrics
        plotStruct = dA.append_2_plotstruct(plotStruct,params,'turns','amplitude','all','mean',frameBaseline); % Turn amplitude
        plotStruct = dA.append_2_plotstruct(plotStruct,params,'turns','duration','all','mean',frameBaseline); % Turn duration
        plotStruct = dA.append_2_plotstruct(plotStruct,params,'turns','amplitude','all','count',frameBaseline); % Turn number
        
        plotStruct = dA.append_2_plotstruct(plotStruct,params,'peran','amplitude','all','mean',frameBaseline); % Crawl amplitude
        plotStruct = dA.append_2_plotstruct(plotStruct,params,'peran','frequency','all','mean',frameBaseline); % Crawl frequency
        
        plotStruct = dA.append_choreography_metric(plotStruct,'speeddiff','mean',frameBaseline,'all'); % Mean para speed
        plotStruct = dA.append_choreography_metric(plotStruct,'curve','mean',frameBaseline,'all');  % Mean curve
        plotStruct = dA.append_choreography_metric(plotStruct,'crabspeed','mean',frameBaseline,'all'); % mean crabspeed
        
        %%%% First 5 second stim
        plotStruct = dA.append_2_plotstruct(plotStruct,params,'rolls','start','first','mean',frame); % latency to first roll
        plotStruct = dA.append_2_plotstruct(plotStruct,params,'rolls','amplitude','all','mean',frameFirst5); % Roll amplitude
        plotStruct = dA.append_2_plotstruct(plotStruct,params,'turns','amplitude','all','mean',frameFirst5); % Turn amplitude
        
        plotStruct = dA.append_choreography_metric(plotStruct,'speeddiff','mean',frameFirst5,'all'); % Mean para speed
        plotStruct = dA.append_choreography_metric(plotStruct,'curve','mean',frameFirst5,'all'); % Mean curve
        plotStruct = dA.append_choreography_metric(plotStruct,'crabspeed','mean',frameFirst5,'all');  % mean crabspeed
        
        %%%% WHole stim metrics
        plotStruct = dA.append_2_plotstruct(plotStruct,params,'rolls','amplitude','all','mean',frame); % Average roll amplitude
        plotStruct = dA.append_2_plotstruct(plotStruct,params,'rolls','amplitude','all','count',frame); % Number of rolls
        plotStruct = dA.append_2_plotstruct(plotStruct,params,'rolls','duration','all','mean',frame); % Average roll duration
        plotStruct = dA.append_2_plotstruct(plotStruct,params,'rolls','duration','all','sum',frame); % Average roll duration
        plotStruct = dA.append_2_plotstruct(plotStruct,params,'rolls','duration','all','proportion',frame); % Proportion of animals that rolled
        plotStruct = dA.append_2_plotstruct(plotStruct,params,'turns','amplitude','all','mean',frame); % Average turn amplitude
        plotStruct = dA.append_2_plotstruct(plotStruct,params,'peran','amplitude','all','mean',frame); % Average crawl amplitude
        plotStruct = dA.append_2_plotstruct(plotStruct,params,'peran','frequency','all','mean',frame); % Crawl amplitude
        
        plotStruct = dA.append_choreography_metric(plotStruct,'speeddiff','mean',frame,'all');  % Mean para speed
        plotStruct = dA.append_choreography_metric(plotStruct,'curve','mean',frame,'all'); % Mean curve
        plotStruct = dA.append_choreography_metric(plotStruct,'crabspeed','mean',frame,'all'); % mean crabspeed
        
        plotStruct = dA.append_choreography_metric(plotStruct,'curve','firstpeak',frame,'all'); % First curve peak
        plotStruct = dA.append_choreography_metric(plotStruct,'crabspeed','firstpeak',frame,'all'); % first crabspeed peak
        
        %%%% Comparison
%         plotStruct = dA.append_2_plotstruct(plotStruct,params,'peran','amplitude','normalise','mean',frame5After); % Normlise crawl amp after stim
%         plotStruct = dA.append_2_plotstruct(plotStruct,params,'peran','frequency','all','mean',frame5After); % Peristalsis freq after stim
%         plotStruct = dA.append_choreography_metric(plotStruct,'speeddiff','mean',frame5After,'normalise'); % Normalised speed diff after stim
        
        plotStruct = dA.append_2_plotstruct(plotStruct,params,'peran','amplitude','normalise','mean',frame15After); % Normlise crawl amp after stim
        plotStruct = dA.append_2_plotstruct(plotStruct,params,'peran','frequency','normalise','mean',frame15After);
%         plotStruct = dA.append_2_plotstruct(plotStruct,params,'peran','amplitude','all','mean',frame15After); % Peristalsis freq after stim
%         plotStruct = dA.append_2_plotstruct(plotStruct,params,'peran','frequency','all','mean',frame15After); % Peristalsis freq after stim
        plotStruct = dA.append_choreography_metric(plotStruct,'speeddiff','mean',frame15After,'normalise'); % Normalised speed diff after stim
        
        %%%% Other
        plotStruct = dA.append_2_plotstruct(plotStruct,params,'area','duration','all','mean',frame); % Animal area (mm^2)

        clear temp dA
    end
    
    if numel(unique([plotStruct.labels])) < n_exp
        effectors = [dcS(idx).effector];
        [~,ia,~] = unique(grp_exp);
        newLabels = arrayfun(@(x) strcat(x.labels," > ",effectors(ia)), plotStruct, 'UniformOutput', false);
        [plotStruct.labels] = newLabels{:};
    end
   
    %%%%% Remove non-overlapping dates
    for jj = 1:numel(plotStruct)
        pos = find(strcmp(plotStruct(jj).labels,'attp2'));
        if isempty(pos) | [numel(pos) > 1]
            continue
        end
        f = plotStruct(jj).group == pos;
        dateFilt = ismember(plotStruct(jj).dates,plotStruct(jj).dates(~f));
        plotStruct(jj).y = plotStruct(jj).y(dateFilt);
        plotStruct(jj).group = plotStruct(jj).group(dateFilt);
        plotStruct(jj).timestamp = plotStruct(jj).timestamp(dateFilt);
        plotStruct(jj).dates = plotStruct(jj).dates(dateFilt);
        try
            plotStruct(jj).id = plotStruct(jj).id(dateFilt);
            plotStruct(jj).uniID = plotStruct(jj).uniID(dateFilt);
        end
    end
    %%%%    
    figPrefix = cellfun(@(x) x(end-3:end),cellstr(unique([dcS.driver])),'UniformOutput', false);
    figPrefix(1:end-1) = strcat(figPrefix(1:end-1),'-');
    figPrefix = strcat(figPrefix{:});
    figurePath = fullfile(dcS(find(idx,1)).figure_directory,figPrefix,'grouped');
    if ~isdir(figurePath)
        mkdir(figurePath)
    end
    
    %%%% Set up loop

    lims = struct("area_",[0 6],...
        "crabspeed_m",[0 2],...
        "curve_m",[0 40],...
        "speed_",[0 2.2],...
        "speeddiff_",[0 3],...
        "peran_amplitude", [0 3],...
        "peran_frequency", [0 2.5]);

    %% sum_dur_roll
    for jj = 1:numel(plotStruct)
        if isempty(plotStruct(jj).y)
            continue
        end
%         label = strsplit(plotStruct(jj).name,'_');
        filename = fullfile(figurePath,plotStruct(jj).name);
        [ax,~,stats] = dataClass.plot_violins(plotStruct(jj),'label',false,{'attp2'},colourRef,plotStruct(jj).name);
        if isempty(ax)
            close all
            continue
        end
        metric = plotStruct(jj).name;
        if startsWith(metric,fieldnames(lims))
            flds = fieldnames(lims);
            f = cellfun(@(x) startsWith(metric,x),flds);
            ylim([lims.(flds{f})]);
        end
        fix_axes(gcf,14,true)
        set(ax.XAxis,'FontSize',8)
        set(ax,'TickLength',[0,0])
        dataClass.save_figure_catch(ax,filename);
        close
        
        save_stats_struct(filename,stats)
        
        %%%%% Normalised Figures (by same-effector control)
        filename = fullfile(figurePath,strcat(plotStruct(jj).name,'_diff'));
        [ax,~,stats] = dataClass.plot_violins(plotStruct(jj),'label',true,{'attp2'},colourRef,plotStruct(jj).name);
        if isempty(ax)
            close all
            continue
        end
        get(gca,'YLim');
        fix_axes(gcf,14,true)
        set(ax.XAxis,'FontSize',8)
        set(ax,'TickLength',[0,0])
        dataClass.save_figure_catch(ax,filename);
        close
        
        save_stats_struct(filename,stats)
        
    end
end



% 
% %%
% 
% cmap = cbrewer('qual','Dark2',8);
% 
% highlightStruct = struct('attp2',[1 1 2.55]./2.55,...
%             'GMR_SS01816',cmap(1,:),...
%             'GMR_SS01817',cmap(1,:),...
%             'GMR_SS01750',cmap(2,:),...
%             'GMR_SS04232',cmap(2,:),...
%             'GMR_SS00666',cmap(2,:),...
%             'GMR_SS04052',cmap(2,:),...
%             'GMR_SS01951',cmap(4,:),...
%             'GMR_SS04245',cmap(4,:),...
%             'GMR_SS04189',cmap(4,:),...
%             'GMR_SS04248',cmap(5,:));
% 
% 
% genoFilter = ismember([dc.driver],PlotGenos);
% dcS = dc(genoFilter);
% 
% [n_prots,grp_prots] = dcS.get_group_numbers('effector','protocol');
% for ii = 1:n_prots
%     plotStruct = struct();
%     
%     idx = grp_prots == ii;
%         
%     [n_exp,grp_exp] = dcS(idx).get_group_numbers('driver','effector','protocol');
%     for jj = 1:n_exp
%         ix = find(idx)';
%         f = grp_exp' == jj;
%         ix = ix(f);
%         temp = dcS(ix).load_groups(false);
% %         temp = dc(myfilt).load_groups(false);
%         dA = temp.filter_by_size(params.area_cutoff);
%         if isempty(dA.timeseries) | ~isfield(dA.timeseries,'et')
%             continue
%         end
%         
%         p = dA.parse_protocol();
%         frame = p(1).start + [0 p(1).length];
%         frameBef = [-15 0] + p(1).start;
%         frameDur = p(1).start + [0 p(1).length];
%         frameBtw = p(1).start + p(1).length + [0 p(1).interval];
%         
%         plotStruct = dA.append_choreography_metric(plotStruct,'speed','mean',frameBef,'all');
%         plotStruct = dA.append_choreography_metric(plotStruct,'curve','mean',frameBef,'all');
%         plotStruct = dA.append_choreography_metric(plotStruct,'crabspeed','mean',frameBef,'all');
%         
%         plotStruct = dA.append_choreography_metric(plotStruct,'speed','mean',frameDur,'all');
%         plotStruct = dA.append_choreography_metric(plotStruct,'curve','mean',frameDur,'all');
%         plotStruct = dA.append_choreography_metric(plotStruct,'crabspeed','mean',frameDur,'all');  
%         
%         plotStruct = dA.append_choreography_metric(plotStruct,'speed','mean',frameBtw,'all');
%         plotStruct = dA.append_choreography_metric(plotStruct,'curve','mean',frameBtw,'all');
%         plotStruct = dA.append_choreography_metric(plotStruct,'crabspeed','mean',frameBtw,'all');  
%     end
%     
%     figPrefix = cellfun(@(x) x(end-3:end),PlotGenos(2:end),'UniformOutput', false);
%     figPrefix(1:end-1) = strcat(figPrefix(1:end-1),'-');
%     figPrefix = strcat(figPrefix{:});
%     figurePath = fullfile(dA.figure_directory,figPrefix,'grouped');
%     if ~isdir(figurePath)
%         mkdir(figurePath)
%     end
%     
%     for jj = 1:numel(plotStruct)
%         label = strsplit(plotStruct(jj).name,'_');
%         [ax,~] = dataClass.plot_violins(plotStruct(jj),'descend',false,{'attp2'},highlightStruct, strcat(label(1),'_',label(3)) );
%         fix_axes(gcf,14,true)
%         set(ax.XAxis,'FontSize',4)
%         set(ax,'TickLength',[0,0])
%         dataClass.save_figure_catch(ax,...
%             fullfile(figurePath,strcat(plotStruct(jj).name) ));
%         close
% 
%         [ax,~] = dataClass.plot_violins(plotStruct(jj),'descend',true,{'attp2'},highlightStruct, strcat(label(1),'_',label(3),'_norm') );
%         fix_axes(gcf,14,true)
%         set(ax.XAxis,'FontSize',4)
%         set(ax,'TickLength',[0,0])
%         dataClass.save_figure_catch(ax,...
%             fullfile(figurePath,strcat(plotStruct(jj).name,'_norm') ));
%         close
%     end
% end
%     
% 
% %%
% 
% 
% 
% cmap = cbrewer('qual','Dark2',8);
% 
% highlightStruct = struct('attp2',[1 1 2.55]./2.55,...
%             'GMR_SS01816',cmap(1,:),...
%             'GMR_SS01817',cmap(1,:),...
%             'GMR_SS01750',cmap(2,:),...
%             'GMR_SS04232',cmap(2,:),...
%             'GMR_SS00666',cmap(2,:),...
%             'GMR_SS04052',cmap(2,:),...
%             'GMR_SS01951',cmap(4,:),...
%             'GMR_SS04245',cmap(4,:),...
%             'GMR_SS04189',cmap(4,:),...
%             'GMR_SS04248',cmap(5,:));
% 
% genoFilter = ismember([dc.driver],PlotGenos);
% dcS = dc(genoFilter);
% 
% [n_prots,grp_prots] = dcS.get_group_numbers('effector','protocol');
% for ii = 1:n_prots
%     idx = grp_prots == ii;
%     
%     plotStruct = struct();
%     
%     [n_exp,grp_exp] = dcS(idx).get_group_numbers('driver','effector','protocol');
%     for jj = 1:n_exp
%         ix = find(idx)';
%         f = grp_exp' == jj;
%         ix = ix(f);
%         temp = dcS(ix).load_groups(false);
% %         temp = dc(myfilt).load_groups(false);
%         dA = temp.filter_by_size(params.area_cutoff);
%         if isempty(dA.behaviour) | isempty(dA.timeseries) | ~isfield(dA.timeseries,'et')
%             continue
%         end
%         
%         p = dA.parse_protocol();
%         frame = p(1).start + [0 p(1).length];
%         
%         plotStruct = dA.append_2_plotstruct(plotStruct,params,'rolls','duration','all','sum',frame);
%         plotStruct = dA.append_2_plotstruct(plotStruct,params,'rolls','duration','all','mean',frame);
%         plotStruct = dA.append_2_plotstruct(plotStruct,params,'rolls','start','first','mean',frame);
%         plotStruct = dA.append_2_plotstruct(plotStruct,params,'rolls','amplitude','all','mean',frame);
%         plotStruct = dA.append_2_plotstruct(plotStruct,params,'rolls','duration','all','proportion',frame);
%         plotStruct = dA.append_2_plotstruct(plotStruct,params,'rolls','duration','all','count',frame);
%         plotStruct = dA.append_2_plotstruct(plotStruct,params,'turns','amplitude','all','mean',frame);
%         plotStruct = dA.append_2_plotstruct(plotStruct,params,'turns','duration','all','mean',frame);
%         plotStruct = dA.append_2_plotstruct(plotStruct,params,'turns','start','first','mean',frame);
%         plotStruct = dA.append_2_plotstruct(plotStruct,params,'area','duration','all','mean',frame);
%         
% 
%         clear temp dA
%     end
% 
%    
%     %%%%    
%     figPrefix = cellfun(@(x) x(end-3:end),PlotGenos(2:end),'UniformOutput', false);
%     figPrefix(1:end-1) = strcat(figPrefix(1:end-1),'-');
%     figPrefix = strcat(figPrefix{:});
%     figurePath = fullfile(dc(find(idx,1)).figure_directory,figPrefix,'grouped');
%     if ~isdir(figurePath)
%         mkdir(figurePath)
%     end
%     
%     %%%% Set up loop
%     names = {'sum_dur_roll','ave_dur_roll',...
%                 'ave_lat_roll','ave_amp_roll',...
%                 'ave_prp_roll','ave_num_roll',...
%                 'ave_amp_turn','ave_dur_turn',...
%                 'ave_lat_turn','ave_area'};
%             
%     xlabels = {"Total Roll Duration (s)","Average Roll Duration (s)",...
%                 "Average Latency to Roll (s)","Average Roll Amplitude",...
%                 "Roll Proportion","Average Number of Rolls",...
%                 "Average Turn Amplitude (s)","Average Turn Duration (s)",...
%                 "Average Latency to Turn (s)","Average Area (mm^{2})"};
% 
%     %% sum_dur_roll
%     for jj = 1:numel(plotStruct)
%         [ax,~] = dataClass.plot_violins(plotStruct(jj),'descend',false,{'attp2'},highlightStruct,xlabels{jj});
%         fix_axes(gcf,14,true)
%         set(ax.XAxis,'FontSize',4)
%         set(ax,'TickLength',[0,0])
%         dataClass.save_figure_catch(ax,...
%             fullfile(figurePath,strcat(names{jj})));
%         close
% 
%         [ax,~] = dataClass.plot_violins(plotStruct(jj),'descend',true,{'attp2'},highlightStruct,strcat('\Delta ',xlabels{jj}));
%         fix_axes(gcf,14,true)
%         set(ax.XAxis,'FontSize',4)
%         set(ax,'TickLength',[0,0])
%         dataClass.save_figure_catch(ax,...
%             fullfile(figurePath,strcat(names{jj},'_norm')));
%         close
%     end
%     
% end
    


%%


    
    