%% plot_grouped_beeswarm_v2
% Subfunctions: 
% MAT-files required: 
% See also:

% Author: Alastair Garner
% email: alastairgarner@outlook.com
% Website: https://github.com/alastairgarner/
% Sep 2019; Last revision: 

%%
function [ax,figure_path] = plot_grouped_beeswarm_v2(plot_struct,sorted,normalise_control_tf,control_driver,highlight_driver,y_label)

%     plot_struct = sum_dur_roll;
%     normalise_control_tf = true;
%     sorted = 'descend'
%     control_driver = {'attp2'}
%     highlight_driver = {'GMR_SS01816'}
%     y_label = "\Delta Total Roll Duration (s)";

    y = plot_struct.y;
    grp = plot_struct.group;
    labels = plot_struct.labels;
    
    if normalise_control_tf & any(plot_struct.labels == control_driver)
        temp = plot_struct;
        ctrl = find(temp.labels == control_driver);
        ff = temp.group == ctrl;

        dates = temp.dates;
        [~,~,ix_a] = unique(dates,'stable');
        ix_c = ix_a(ff);

    %     accumarray(ix_a,temp.y',[max(ix_a) 1],@mean,NaN)
        ctrl_val = accumarray(ix_c,temp.y(ff)',[max(ix_a) 1],@mean,NaN);
        y = [temp.y-ctrl_val(ix_a)'];
    end
    
    control_col = [1 1 2.55]./2.55;
    highlight_col = [1 0 0];
    labels = strrep(labels,"GMR_","");
    
    dash_width = .6;
    
    counts = accumarray(grp',1);
    means = accumarray(grp',y',[],@(x) mean(x,'omitnan'));
    err = accumarray(grp',y',[],@(x) std(x,'omitnan'))./sqrt(counts);
    switch sorted
        case 'none'
            order = [1:numel(unique(plot_struct.group))];
        case 'ascend'
            [~,order] = sort(means,'ascend');
        case 'descend'
            [~,order] = sort(means,'descend');
    end
    
    cnts = accumarray([plot_struct.group]',1);    
    yy = mat2cell(y,1,counts);
    
    ave_tf = false;
    if ~ave_tf
        yy = yy(order);
    end
    
    y_mean = repelem(means(order),1,2);
    x_mean = [1:numel(order)]' + ([-1 +1]*dash_width/2);
    y_err = y_mean + ([-1 +1].*err(order));
    x_err = [1:numel(order)] .* [1;1];
    
    hold on
    p = plotSpread(yy,'distributionColors',[100 100 100]/255);
    e = plot(x_err,y_err','k','LineWidth',.5);
    m = plot(x_mean',y_mean','k','LineWidth',1.5);
    hold off
    
    if all(isnan(p{1}))
        ax = gca;
        figure_path = [];
        return
    end        
    
    obs = findobj(gca,'Marker','.');
    plotted = fliplr(double(string({obs.DisplayName})));
    [~,control_pos] = find(fliplr(plot_struct.labels(order(plotted))) == string(control_driver)');
    [~,highlight_pos] = find(fliplr(plot_struct.labels(order(plotted))) == string(highlight_driver)');
%     [~,control_pos] = find(fliplr(plot_struct.labels(order)) == string(control_driver)');
%     [~,highlight_pos] = find(fliplr(plot_struct.labels(order)) == string(highlight_driver)');
    control_col = repelem({control_col},1,numel(control_pos));
    highlight_col = repelem({highlight_col},1,numel(highlight_pos));
    [obs(control_pos).Color] = control_col{:};
    [obs(highlight_pos).Color] = highlight_col{:};
    
    xlim([0 max(grp)+1])
    pbaspect([4,1,1])
    
    set(gca,'XTickLabel',labels(order),'XTickLabelRotation',45,'LineWidth',1,'TickLabelInterpreter','none');
    obs = findobj(gca,'Marker','.');
    set(obs,'MarkerSize',8);
    ylabel(y_label)
    
    to_top = findobj(gca,'LineStyle','-');
    to_bot = findobj(gca,'-not','LineStyle','-','-not','Type','Axes');
    set(gca,'Children',[to_top;to_bot])
    
    [pval,hval,stats] = mult_comp_wilcox(plot_struct);
    stars = sum(pval(order) < [.05, .01, .005]',1);
    xvals = find(stars);
    stars = stars(xvals);
    stars = arrayfun(@(x) repelem('*',x),stars,'UniformOutput', false)';
    
    stepsize = range(get(gca,'XLim'))/numel(pval);
%     xvals = xvals*stepsize;
    yval = max(get(gca,'YLim')) + [range(get(gca,'YLim'))*.02];
    yval = repelem(yval,numel(xvals));
    
    text(xvals,yval,string(stars),'HorizontalAlignment','center','FontSize',15)
    
    ax = gca;
    figure_path = [];
end