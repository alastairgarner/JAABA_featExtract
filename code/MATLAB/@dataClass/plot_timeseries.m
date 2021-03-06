%% plot_timeseries

% Other m-files required: 
% Subfunctions: 
% MAT-files required: 
% See also:

% Author: Alastair Garner
% email: alastairgarner@outlook.com
% Website: https://github.com/alastairgarner/
% Sep 2019; Last revision: 

%%
function [ax,fig_path] = plot_timeseries(obj,params,feature,frame_start,frame_end,y_limits,calculation,colormap)
    fr1 = frame_start;
    fr2 = frame_end;
    calc = calculation;
    cmap = colormap;
    
    linestyle = repelem({'-'},size(cmap,1));
    [~,~,ic] = unique(cmap, 'rows', 'stable');
    linestyle([false; diff(ic)==0]') = {':'};
        
    for xx = 1:numel(obj)
        data = obj(xx);
        if numel(fieldnames(data.timeseries)) == 0 | ~numel(data.timeseries)
            ax = [];
            fig_path = [];
            continue
        end
        counts = arrayfun(@(x) numel(x.et), data.timeseries, 'ErrorHandler', @(a,b) nan);
        
        cellDiff = cellfun(@(x,y) x-y,{data.timeseries.speed},{data.timeseries.crabspeed}, 'UniformOutput', false);
        [data.timeseries.speeddiff] = cellDiff{:};
        
        switch calc
            case "norm"
                metric = [data.timeseries.(feature)];
            case "diff"
                filt = arrayfun(@(x) ~isempty(x.et), data.timeseries);
                metric = arrayfun(@(x) [0 diff(x.(feature))], data.timeseries(filt), 'UniformOutput', false);
                metric = [metric{:}];
            otherwise
                metric = [data.timeseries.curve];
        end

        et = [data.timeseries.et];
        ids = [data.timeseries.uniID];
        ids = repelem(ids,1,counts(counts>0));

        bins = [0:.2:ceil(max(et))];

        Y = discretize(et,bins);
        y = accumarray(Y',metric',[numel(bins) 1],@mean);
        n = accumarray(Y',ids',[numel(bins) 1],@(x) numel(unique(x)));
        dev = accumarray(Y',metric',[numel(bins) 1],@std);

        err = dev./sqrt(n);
        err(isnan(err)) = 0;

        err_y = [y+err;flipud(y-err)]';
        err_x = [bins(1:end),fliplr(bins(1:end))];

        hold on
        plot(bins(1:end),y,linestyle{xx},'Color',cmap(xx,:), 'LineWidth', 3, 'Tag', string(xx))
        patch(err_x,err_y,cmap(xx,:),'EdgeColor','none','FaceAlpha',.3);
        hold off
        
    end
    
    if isempty(findobj('Type','Figure'))
        ax = [];
        fig_path = [];
        return
    end

    % fix Z order
    lines = findobj(gca,'Type','Line');
    other = findobj(gca,'-not','Type','Line','-not','Type','Axes');
    [~,I] = sort(double(string({lines.Tag})),'descend');
    set(gca,'Children',[lines(I);other]);

    % clean up axes
    xlim([fr1 fr2])
    ylim(y_limits)
    ylabel(strcat(feature," ",calc));
    xlabel("Time (s)")
    pbaspect([3,1,1])
    
    %
    obj(1).plot_stimulation_time(params);
    %
    disp_names = {obj.driver};
    if numel(unique(cellstr(disp_names))) < numel(disp_names)
        disp_names = strcat(string({obj.driver}),' > ',string({obj.effector}));
    end
    disp_names = cellfun(@(x) strrep(x,'GMR_',''), disp_names,'UniformOutput', false);
    
    relative_position = [.7 1.05 .3 .2];
    plts = findobj(gca,'Type','Line');
%     flexi_legend(flipud(plts),disp_names,'RelativePosition',relative_position,'Interpreter', 'none','Box','off');
    legend(gca,flipud(plts),disp_names,'Interpreter', 'none','Box','off');
%     set(gca,'LegendColorbarListeners',[]);
    
    % set up save name and path 
    figure_dir = obj.figure_directory;
    fig_type = strcat("timeseries_",feature,"_",calc);
    fname = strcat(obj(end).get_full_genotype());

    fig_path = fullfile(figure_dir,'byDay',fname,fig_type);
    if ~isdir(fileparts(fig_path))
        mkdir(fileparts(fig_path));
    end

    ax = gca;
    
end