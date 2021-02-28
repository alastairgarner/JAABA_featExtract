%% plot_violins

% Other m-files required: 
% Subfunctions: 
% MAT-files required: 
% See also:

% Author: Alastair Garner
% email: alastairgarner@outlook.com
% Website: https://github.com/alastairgarner/
% Sep 2019; Last revision: 

function [ax,figure_path,stats] = plot_violins(plot_struct,sortOrder,normalise_control_tf,control_driver,highlightStruct,y_label)

    args.Bandwidth = [];
    args.Width = 0.4;
    args.BoxWidth = 0.15;
    args.PlotPoints = false;
    args.PlotViolins = true;
    args.AverageType = 'median';
    args.StatsTest = 'wilcoxon-stringent';
    args.ScatterSize = 1;
    args.TypeHighlight = 'violin';
    args.widthRatio = 1/20;
    highlightColor = [1 0 0];
    args.AverageMarkerSize = 3;
    args.AstericksRotation = 90;
    args.AstericksHorzAlign = 'left';
    args.AstericksVertAlign = 'middle';
    args.AstericksAdjust = 0.5;
    args.Annotation = false;
    args.plotCI = true;
    
    if isempty(plot_struct.y)
        ax = [];
        figure_path = [];
        stats = [];
        return
    end
    
    %%%%%%% Run Statistics
%     if strcmp(args.StatsTest,'fisher') | islogical(plot_struct.y)
%         stats = run_stats_test(plot_struct,'Test','fisher-bonferroni');
%         args.StatsTest = 'fisher';
%     elseif strcmp(args.StatsTest,'wilcoxon-stringent')
%         stats = run_stats_test(plot_struct,'Test','wilcoxon-bonferroni');
%     elseif sum(strcmp(plot_struct.labels,'attp2')) ~= 1
%         stats = run_stats_test(plot_struct,'Test','kruskal-dunn');
%     end

    goodGroups = unique(plot_struct.group(~isnan(plot_struct.y)),'stable');
    goodControl = ismember(goodGroups,find(contains(plot_struct.labels,'attp2')));
    goodControl = any(goodControl);
    
    %%%%%%% Run Statistics
    if [strcmp(args.StatsTest,'fisher') | islogical(plot_struct.y)]
        stats = run_stats_test(plot_struct,'Test','fisher-bonferroni');
        args.StatsTest = 'fisher';
    end
    
    %%%%% Switch to Fishers Exact, if logical array is presented
    if islogical(plot_struct.y)
        args.PlotPoints = true;
        args.PlotViolins = false;
        args.StatsTest = 'fisher';
        args.AverageType = 'mean';
        args.ScatterSize = 8;
        args.TypeHighlight = 'scatter';
        
        data = plot_struct.y;
        grp = plot_struct.group;
        tstamps = plot_struct.timestamp;        
        [C,ia,ic] = unique([grp',tstamps'],'stable','rows');
        
%         [~,~,idxDate] = unique([grp;plot_struct.dates]','stable','rows');
%         plot_struct.yByDay = accumarray(idxDate,data',[max(idxDate) 1],@(x) mean(x,'omitnan'))';
%         plot_struct.yByDay(idxDate)
        
        plot_struct.group = C(:,1)';
        plot_struct.timestamp = C(:,2)';
        plot_struct.dates = plot_struct.dates(ia);
        plot_struct.n = accumarray(ic,1)';
        plot_struct.y = accumarray(ic,data',[max(ic) 1],@(x) mean(x,'omitnan'))';
    end
    
    %%%%% Normalise data by median/mean of control
    if normalise_control_tf & goodControl
        plot_struct.yRaw = plot_struct.y;
        Data = nan(1,numel(plot_struct.y)); % set up blank array to fill with normalised values
        %%% Check if there are >1 controls to normalise by
        str = split([plot_struct.labels],'>');
        if size(str,3) > 1
            [~,~,ic] = unique(str(:,:,end),'stable');
        else
            ic = ones(numel(str(:,:,1)),1);
        end
        %%% Normalise by attp2 for each effector
        for xx = 1:size(str,3)
            labelFilt = ic==xx;
            dataFilt = ismember(plot_struct.group,find(labelFilt));
            ctrl = find(contains(plot_struct.labels,control_driver)' & labelFilt);
            controlFilt = ismember(plot_struct.group,ctrl);
            
            dates = plot_struct.dates;
            [~,~,idxAll] = unique(dates,'stable');
            idxCont = idxAll(controlFilt);
            
            if strcmp(args.StatsTest,'fisher')
                ctrl_val = accumarray(idxCont,plot_struct.y(controlFilt)',[max(idxAll) 1],@(x) mean(x,'omitnan'),NaN);
            else
                ctrl_val = accumarray(idxCont,plot_struct.y(controlFilt)',[max(idxAll) 1],@(x) median(x,'omitnan'),NaN);
            end
            
            Data(dataFilt) = [plot_struct.y(dataFilt)-ctrl_val(idxAll(dataFilt))'];

            %%% Update the stats struct with the adjusted values
            if strcmp(args.StatsTest,'fisher')
                df = dataFilt;
%                 cellFilters = arrayfun(@(x) ismember(plot_struct.dates(df),plot_struct.dates(df & plot_struct.group==x)) & controlFilt(df),...
%                     [find(labelFilt)'], 'UniformOutput', false);
                cellFilters = arrayfun(@(x) ismember(plot_struct.dates,plot_struct.dates(plot_struct.group==x)) & controlFilt,...
                    [find(labelFilt)'], 'UniformOutput', false);
                ctrlByDay = cellfun(@(x) sum(plot_struct.y(x).*plot_struct.n(x)) ./ sum(plot_struct.n(x)) ,cellFilters);
                                
                statsMat = [stats(labelFilt).Average; stats(labelFilt).CILower; stats(labelFilt).CIUpper]';
                statsMatAdj = num2cell(statsMat-stats(ctrl).Average);
                statsMatAdj = num2cell(statsMat-ctrlByDay');
                [stats(labelFilt).Average] = statsMatAdj{:,1};
                [stats(labelFilt).CILower] = statsMatAdj{:,2};
                [stats(labelFilt).CIUpper] = statsMatAdj{:,3};
            end
        end
        
    else
        Data = plot_struct.y;

    end
    plot_struct.y = Data;
    
    %%%%%%% Run Statistics
    if [sum(contains(plot_struct.labels,'attp2'))] ~= 1 & [~strcmp(args.StatsTest,'fisher')]
        stats = run_stats_test(plot_struct,'Test','kruskal-dunn');
    elseif strcmp(args.StatsTest,'wilcoxon-stringent')
        stats = run_stats_test(plot_struct,'Test','wilcoxon-bonferroni');
    end
    
    %%%%%%% Check there is good data to plot
    Groups = [plot_struct.group];
    [uniGroups,~,idxGroups] = unique(Groups,'stable');
    Average = [stats.Average];
    if all(isnan(Average))
        ax = [];
        figure_path = [];
        stats = [];
        return
    end
    
    %%%%%%% Set the plot order
    if isnumeric(sortOrder)
        sortIndex = sortOrder;
    elseif ismember(sortOrder,{'descend','ascend'})
        [~,sortIndex] = sort(Average,sortOrder);
    else
        [~,sortIndex] = sort(lower(plot_struct.labels(unique(plot_struct.group,'stable'))));
    end
    
    %%%%%%%%%% Determine appropriate Y limits
    yMax = accumarray(Groups',Data',[numel(uniGroups) 1],@(x) ceil(quantile(x,.95)));
    yMin = accumarray(Groups',Data',[numel(uniGroups) 1],@(x) floor(quantile(x,.05)));
    yLims = [min(yMin),max(yMax)];
    if yLims(1)==yLims(2)
        if yLims(1)==0
            yLims = [-0.5 0.5];
        else
            yLims(1) = 0;
        end
    end
    
    %%%%%%%%%%%
    nGroups = numel(uniGroups);
    args.widthRatio = nGroups*args.widthRatio;

    %%%%%%%% Set cutom plot size for groups with less than 8 genotypes
    if numel(unique(plot_struct.group)) <= 7
        args.Width =0.4;
        args.BoxWidth = 0.05;
        args.widthRatio = 0.5;
        args.AverageMarkerSize = 5;
        args.AstericksRotation = 0;
        args.AstericksHorzAlign = 'center';
        args.AstericksVertAlign = 'middle';
        args.AstericksAdjust = 0;
        args.Annotation = true;
    end
    
    %%%%%%% Set plot size
    ax = gca;
    pb = get(ax,'Position');
    boxPoints = pb+[0,0,pb(1:2)];
    midPoints = mean(reshape(boxPoints',2,[]),2)';
    newPoints = abs([midPoints,0,0] - [0.5 0.5 1 1].*repmat(0.2*[args.widthRatio 1],1,2));
    set(ax, 'Position',newPoints)

    %%%%%%%% Sort Stats struct to match plot order
    statsIdx = arrayfun(@(x) find([stats.index]==x,1),[1:numel(uniGroups)]);
    statsOrdered = stats(statsIdx);
    statsOrdered = statsOrdered(sortIndex);
        

    hold('on');
    for ii = 1:numel(uniGroups)

        f = Groups == sortIndex(ii);
        data = Data(f);
        groups = Groups(f);
        timestamps = plot_struct.timestamp(f);
    %     pos = mean(Groups(f));
        pos = ii;
        
        ave = statsOrdered(pos).Average;
        if args.plotCI  
            boxY = repelem([statsOrdered(pos).CILower,statsOrdered(pos).CIUpper],1,2); % 95% Confidence Intervals
        else 
            boxY = statsOrdered(pos).Average + [statsOrdered(pos).STD.*[-1 -1 1 1]]; % Standard Deviation
        end

        if all(isnan(data)) | isnan(ave)
%             "Hello"
            continue
        end

        % calculate kernel density estimation for the violin
        [density, value] = ksdensity(data, 'bandwidth', args.Bandwidth);
        density = density(value >= min(data) & value <= max(data));
        value = value(value >= min(data) & value <= max(data));
        value(1) = min(data);
        value(end) = max(data);

        % all data is identical
        if min(data) == max(data)
            density = 1;
        end
        width = args.Width/max(density);
        
        if args.PlotPoints
            % plot the data points within the violin area
            if length(density) > 1
                jitterstrength = interp1(value, density*width, data);
            else % all data is identical:
                jitterstrength = density*width;
            end
            jitter = 2*(rand(size(data))-0.5);
            ScatterPlot = ...
                line(pos + jitter.*jitterstrength, data,...
                    'LineStyle','none', 'Marker','.',...
                    'MarkerSize',args.ScatterSize,'Color',[100 100 100]/255,...
                    'DisplayName','scatter');
            
        end
        
        if args.PlotViolins
            % plot the violin
            ViolinPlot =  ... % plot color will be overwritten later
                patch([pos+density*width pos-density(end:-1:1)*width], ...
                     [value value(end:-1:1)],...
                     [1,0,0],'EdgeColor','none',...
                     'FaceColor',[100 100 100]/255,'FaceAlpha',.6,...
                     'DisplayName','violin');

        end
        
        % plot boxes (Confidence Int vs Interquartile Range vs S.E.M)
        BoxPlot = ... % plot color will be overwritten later
            patch(pos+[-1,1,1,-1]*args.BoxWidth,...
                 boxY,...
                 [1 1 1],'EdgeColor','none','FaceColor',[0,0,0],...
                 'DisplayName','box');
                 
        % plot the data median
        MedianPlot = line(pos, ave,...
            'LineStyle','none', 'Marker','o','MarkerSize',args.AverageMarkerSize,...
            'MarkerEdgeColor','k','MarkerFaceColor','w',...
            'DisplayName','median');

    end
    hold('off');

    Labels = plot_struct.labels(sortIndex);
    controlPosition = contains(Labels,control_driver);
    highlightPosition = cellfun(@(x) contains(Labels,x), fieldnames(highlightStruct),'UniformOutput', false);
    highlightPosition = vertcat(highlightPosition{:}).*[1:numel(highlightPosition)]';
    highlightPosition = max(highlightPosition,[],1);
        
    highlightColors = struct2cell(highlightStruct);
    highlightOrder = highlightPosition(find(highlightPosition));
    highlightColors = highlightColors(highlightOrder);
    
    cnt = accumarray(highlightOrder',1);
    repeats = arrayfun(@(x) [1:x], cnt, 'UniformOutput', false);
    repeats = fliplr([repeats{:}]);

    cmapFound = vertcat(highlightColors{:});
    highlightColors = cmapFound + [cmapFound./2].*[repeats'-1];
    highlightColors(highlightColors>1) = 1;
    highlightColors(highlightColors<0) = 0;
    highlightColors = num2cell(highlightColors,2);
    
    Labels = strrep(Labels,'GMR_','');

    ax = gca;
    set(ax,'XLim',[0 max(pos)+1],...
        'XTick',1:numel(uniGroups),...
        'XTickLabel',Labels,...
        'XTickLabelRotation',45,...
        'YLim',yLims,...
        'TickLabelInterpreter','none');

    ax.YLabel.Interpreter = 'none';
    ax.YLabel.String = y_label;

    violins = flipud(findobj(ax,'DisplayName',args.TypeHighlight));
    yPosition = arrayfun(@(x) round(mean(x.XData,'omitnan')),violins);
    [~,controlIndex] = ismember(yPosition,find(controlPosition));
    [~,highlightIndex] = ismember(yPosition,find(highlightPosition));
    
    if numel(violins)~=0
        try
            set(violins(logical(controlIndex)),'FaceColor',[1 1 2.55]./2.55,'FaceAlpha',.8);
            set(violins(logical(highlightIndex)),'FaceAlpha',1);
            [violins(logical(highlightIndex)).FaceColor] = highlightColors{:};        
        catch
            set(violins(logical(controlIndex)),'Color',[1 1 2.55]./2.55);
            [violins(logical(highlightIndex)).Color] = highlightColors{:}; 
        end
    end
    
    if strcmp(args.StatsTest,'kruskal-dunn')
        ax = gca;
        figure_path = [];
        return
    end

    %%%%%% Plot significance stars (not Kruskal though)
    stars = sum([statsOrdered.pval] < [.05, .01, .001]',1);
    xvals = find(stars);
    stars = stars(xvals);
    stars = arrayfun(@(x) repelem('*',x),stars,'UniformOutput', false)';

    if args.Annotation & any(xvals)
        pbpos = plotboxpos;
        xlims = get(gca,'XLim');
        ylims = get(gca,'YLim');

        textBoxY = pbpos(2)+pbpos(4);
        textBoxWidth = pbpos(3)./range(xlims)./2;
        positions = pbpos(1) + [[xvals./range(xlims)].*pbpos(3)];
        dim = positions'.*[1 0 0 0] + [-textBoxWidth, textBoxY, 2*textBoxWidth, 0.01];
        annontateFunc = @(x,y) annotation('textbox',x,'String',string(y),...
            'HorizontalAlignment','center',...
            'VerticalAlignment','bottom',...
            'LineStyle','none','Margin',0,...
            'FontSize',12);
        t = cellfun(annontateFunc, num2cell(dim,2),stars);
    else        
        yval = max(get(gca,'YLim')) + [range(get(gca,'YLim'))*.04];
        yval = repelem(yval,numel(xvals));

        text(xvals+args.AstericksAdjust,yval,string(stars),'HorizontalAlignment',args.AstericksHorzAlign,'FontSize',12,...
            'Rotation',args.AstericksRotation,'VerticalAlignment',args.AstericksVertAlign);
    end
    
    ax = gca;
    figure_path = [];

end
