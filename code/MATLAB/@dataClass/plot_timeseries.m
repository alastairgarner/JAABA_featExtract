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
    
    for xx = 1:numel(obj)
        data = obj(xx);
        if numel(fieldnames(data.timeseries)) == 0
            ax = [];
            figure_path = [];
            continue
        end
        counts = arrayfun(@(x) numel(x.et), data.timeseries, 'ErrorHandler', @(a,b) nan);
        
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
        plot(bins(1:end),y,'Color',cmap(xx,:), 'LineWidth', 2, 'Tag', string(xx))
        patch(err_x,err_y,cmap(xx,:),'EdgeColor','none','FaceAlpha',.3);
        hold off
        
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
    pbaspect([2,1,1])
    
    %
    obj(1).plot_stimulation_time(params);
    %
    disp_names = {obj.driver};
    disp_names = cellfun(@(x) strrep(x,'GMR_',''), disp_names,'UniformOutput', false);

    relative_position = [.7 1.05 .3 .2];
    plts = findobj(gca,'Type','Line');
%     flexi_legend(flipud(plts),disp_names,'RelativePosition',relative_position,'Interpreter', 'none','Box','off');
    legend(gca,flipud(plts),disp_names,'Interpreter', 'none','Box','off');
%     set(gca,'LegendColorbarListeners',[]);
    
    % set up save name and path 
    figure_dir = obj.figure_directory;
    fig_type = strcat("timeseries_",feature,"_",calc);
%     fname = strcat(obj(end).get_full_genotype,".pdf");
    fname = strcat(obj(end).get_full_genotype);


    fig_path = fullfile(figure_dir,fig_type,fname);

    ax = gca;
    addlistener(ax,'Position','PostSet',@(~,~) fprintf('Change detected'));
    
end