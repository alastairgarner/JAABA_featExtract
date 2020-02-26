%% plot_ethogram

% Other m-files required: 
% Subfunctions: 
% MAT-files required: 
% See also:

% Author: Alastair Garner
% email: alastairgarner@outlook.com
% Website: https://github.com/alastairgarner/
% Sep 2019; Last revision: 

%%params,{"rolls"},bin_width,[25 50],[-.1 .2],cmap
function [ax,fig_path] = plot_behaviour_timeseries(obj,params,behaviour,bin_width,x_lims,y_lims,cmap)
    fr1 = x_lims(1);
    fr2 = x_lims(2);
    
    tmax = ceil(max(cellfun(@(x) max([x.track_end]),{obj.behaviour}, 'ErrorHandler', @(a,b) nan)));
    bins = [0:bin_width:tmax];
    
    x = []; y = []; err_x = []; err_y = {};
    n = 1;
    for ii = 1:numel(obj)
        data = obj(ii);
        
        if isempty(data.behaviour)
            continue
        end

        [f,~] = find(string({data.behaviour.behaviour})' == string(behaviour));
        for jj = 1:numel(f)
            beh = data.behaviour(f(jj));

            if numel(beh) == 0
                continue
            end

            indicies = data.timestamp_index;
            [B,I] = ismember([beh.uniID],[data.uniID]);
            indicies = indicies(I);
            idx_cell = num2cell(unique(indicies,'stable'));

            cells = cellfun(@(x) get_behaved(beh,bins,[indicies==x])', idx_cell, 'UniformOutput', false);
            behaved_stack = vertcat(cells{:});

            cells = cellfun(@(x) get_tracked(beh,bins,[indicies==x])', idx_cell, 'UniformOutput', false);
            tracked_stack = vertcat(cells{:});

            percentge_stack = behaved_stack./tracked_stack;
            sem = std(percentge_stack,1)./sqrt(numel(idx_cell));

            tracked = get_tracked(beh,bins,[])';
            behaved = get_behaved(beh,bins,[])';

            y{n} = behaved./tracked;
            x{n} = bins(1:end);

            err = [y{n}+sem,fliplr(y{n}-sem)];
            err(isnan(err)) = 0;
            err_x{n} = [bins(1:end),fliplr(bins(1:end))];
            err_y{n} = err;
            n = n+1;
        end
    end
    
    if isempty(x)
        ax = [];
        fig_path = [];
        return
    end
    
    % vert concatonation of data across timestamps
    x = vertcat(x{:});
    y = vertcat(y{:}) *100;
    err_x = vertcat(err_x{:});
    err_y = vertcat(err_y{:}) *100;
    
    hold on
    for ii = 1:size(x,1)
        plot(x(ii,:),y(ii,:),'Color',cmap(ii,:), 'LineWidth', 2)
        patch(err_x(ii,:),err_y(ii,:),cmap(ii,:),'EdgeColor','none','FaceAlpha',.3);
    end
    hold off
    
    % fix Z order
    lines = findobj(gca,'Type','Line');
    other = findobj(gca,'-not','Type','Line','-not','Type','Axes');
    [~,I] = sort(double(string({lines.Tag})),'descend');
    set(gca,'Children',[lines(I);other]);
    
    % clean up axes
    ylabel(strcat(behaviour{:}, " %"));
    xlim(x_lims)
    if isempty(y_lims)
        ylim([0 max(max(err_y))]);
    else
        ylim(y_lims)
    end
    xlabel("Time (s)")
    pbaspect([3,1,1])
    
    %
    obj(1).plot_stimulation_time(params);
    %
    disp_names = {obj.driver};
    disp_names = cellfun(@(x) strrep(x,'GMR_',''), disp_names,'UniformOutput', false);

    relative_position = [.7 1.05 .3 .2];
    plts = findobj(gca,'Type','Line');
%     flexi_legend(flipud(plts),disp_names,'RelativePosition',relative_position,'Interpreter', 'none','Box','off');
    leg = legend(gca,flipud(plts),disp_names,'Interpreter', 'none','Box','off');
    
    % set up save name and path 
    figure_dir = obj.figure_directory;
    fig_type = strcat("beh_timeseries_",behaviour{:});
%     fname = strcat(obj(end).get_full_genotype,".pdf");
    fname = strcat(obj(end).get_full_genotype);

    fig_path = fullfile(figure_dir,fig_type,fname);

    ax = gca;
    
    
    %% Nested Functions
    % get number of animals tracked
    function tracked = get_tracked(beh,bins,filter)
        if isempty(filter)
            filter = repelem(true,1,numel(beh.start));
        end
        [~,idx,~] = unique(beh.uniID);
        idx = intersect(find(filter),idx');
        
        tracked = discretize([beh.track_start(idx);...
            beh.track_end(idx)], bins);
        tracked = tracked(:,~any(isnan(tracked),1));
        tracked = cumsum( accumarray(tracked(1,:)',1,[numel(bins),1])...
            - accumarray(tracked(2,:)',1,[numel(bins),1]) );
    end

    % get number of animals performing the behaviour per bin
    function behaved = get_behaved(beh,bins,filter)
        if isempty(filter)
            filter = repelem(true,1,numel(beh.start));
        end
        [~,idx,~] = unique([beh.uniID;beh.start;beh.end]','rows');
        idx = intersect(find(filter),idx');
        
        behaved = discretize([beh.start(idx);...
            beh.end(idx)], bins);
        behaved = behaved(:,~any(isnan(behaved),1));
        behaved = cumsum( accumarray(behaved(1,:)',1,[numel(bins),1])...
            - accumarray(behaved(2,:)',1,[numel(bins),1]));
    end
end