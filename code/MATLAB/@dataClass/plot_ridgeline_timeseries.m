%% plot_ridgeline_timeseries

% Other m-files required: 
% Subfunctions: 
% MAT-files required: 
% See also:

% Author: Alastair Garner
% email: alastairgarner@outlook.com
% Website: https://github.com/alastairgarner/
% Sep 2019; Last revision: 

%%
function [ax,fig_path] = plot_ridgeline_timeseries(obj,params,feature,frame_start,frame_end,y_scaler)
    if numel(unique(obj.get_full_genotype)) > 1
        fprintf("too many genotypes specified... please filter the data\n")
        return
    end
    if numel(fieldnames(obj.timeseries)) == 0
        ax = []; fig_path = [];
        fprintf("ridgeline - no timeseries data \n")
        return
    end
    fr1 = frame_start;
    fr2 = frame_end;
    
    counts = arrayfun(@(x) numel(x.et), obj.timeseries);

    et = [obj.timeseries.et];
    metric = [obj.timeseries.(string(feature))];
    ids = [obj.timeseries.uniID];
    ids = repelem(ids,1,counts(counts>0));

    func1 = @(x) max([false, (min(x)<fr1 & max(x)>fr2)]);
    func2 = @(x) x>fr1 & x<fr2;

    f1 = cellfun(func1,{obj.timeseries.et}, 'UniformOutput', false);
    f2 = cellfun(@(x,y) x & func2(y), f1, {obj.timeseries.et}, 'UniformOutput', false);
    f2 = [f2{:}];

    et = et(f2);
    metric = metric(f2);
    [uids,~,ids] = unique(ids(f2));
    
    if ~numel(uids)
        ax = []; fig_path = [];
        fprintf("ridgeline - no good animals \n")
        return
    end
    
    numels = numel(unique(ids));
    if numels > 250
        samp_ids = randsample([1:numels],250);
        ismem = ismember(ids,samp_ids);
        
        metric = metric(ismem);
        et = et(ismem);
        [uids,~,ids] = unique(ids(ismem));
    end  
    
    cols = viridis(numel(uids));
    hold on
    for kk = 1:numel(uids)
        f = ids == kk;
%         plot(et(f),(metric(f)/mean(metric))+ids(f)','Color',cols(ii,:),'LineWidth',2);
%         plot(et(f),(metric(f)*scl)+ids(f)','Color','white','LineWidth',2);
        plot(et(f),(metric(f)*y_scaler)+ids(f)','Color',cols(kk,:),'LineWidth',1);
    end
    hold off
    
    max_y = max((metric(f)*y_scaler)+ids(f)') +1;
    
    % reorder lines
    child = get(gca,"Children");
    set(gca,'YTickLabel',[],...
        "Children",flipud(child))

    % add scale bar
    ylm = get(gca,'YLim');
    lowerlim = ylm(1) + (.05*range(ylm));
    xlm = get(gca,'XLim');
    leftlim = xlm(1) + (.1*range(xlm));
    hold on
    text(leftlim,lowerlim,strcat(" ",string(5/y_scaler)),...
        'VerticalAlignment','bottom','HorizontalAlignment','left',...
        'FontSize',8, 'BackgroundColor', 'white');
    plot([leftlim leftlim],[lowerlim lowerlim+5],'k','LineWidth', 2);
    hold off
    %
    ylim([0 max_y*1.04])
    ylabel(feature)
    xlabel('Time (s)')
    pbaspect([1,2,1])
    
    %
    obj(1).plot_stimulation_time(params);
    
    % Save figure
    figure_dir = obj.figure_directory;
    fig_type = strcat("ridge_",feature);
%     fname = strcat(obj.get_full_genotype,".pdf");
    fname = strcat(obj.get_full_genotype);
    
    fig_path = fullfile(figure_dir,fig_type,fname);
    if ~isdir(fileparts(fig_path))
        mkdir(fileparts(fig_path));
    end
    
    % Define axis as output
    ax = gca;
    
end