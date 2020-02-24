%% plot_ethogram

% Other m-files required: 
% Subfunctions: 
% MAT-files required: 
% See also:

% Author: Alastair Garner
% email: alastairgarner@outlook.com
% Website: https://github.com/alastairgarner/
% Sep 2019; Last revision: 

%%
function [ax,fig_path] = plot_ethogram(obj,params,behaviours_cell,frame_start,frame_end,offset)
    fr1 = frame_start;
    fr2 = frame_end;
    behaviours = behaviours_cell;
    
    if isempty(obj.behaviour)
        ax = []; fig_path = [];
        fprintf("ethogram - no behaviour data \n")
        return
    end
    
    f = string({obj.behaviour.behaviour})' == string(behaviours);
    [idx,~] = find(f);
    if ~any(f)
        ax = []; fig_path = [];
        fprintf("ethogram - specified behaviours not present \n")
        return
    end
    behs = obj.behaviour(idx);

    counts = arrayfun(@(x) numel(x.id), behs);
    grp = repelem([1:numel(counts)],1,counts);
    
    tS = [behs.track_start];
    tE = [behs.track_end];
    starts = [behs.start];
    ends = [behs.end];
    [~,~,ids] = unique([behs.uniID],'stable');

    if ~isempty(fr1) & ~isempty(fr2)
        f1 = (tS<fr1 & tE >fr2);
        f2 = (starts>fr1 & starts<fr2) | (starts<fr1 & ends>fr1) | (starts<fr2 & ends>fr2);
        f = f1 & f2;

        grp = grp(f);
        starts = starts(f);
        ends = ends(f);
        [~,~,ids] = unique(ids(f));
    else
        fr1 = floor(min(tS));
        fr2 = ceil(max(tE));
    end
    
    if isempty(starts)
        fprintf("NONE of the specified behaviours are present \n")
        ax = [];
        fig_path = [];
        return
    end
    
    numels = numel(unique(ids));
    if numels > 500
        samp_ids = randsample([1:numels],500);
        ismem = ismember(ids,samp_ids);
        
        grp = grp(ismem);
        starts = starts(ismem);
        ends = ends(ismem);
        [~,~,ids] = unique(ids(ismem));
    end
    
    cols = cbrewer('qual','Dark2',8);
    cols = vertcat([0.7 0.7 0.7],cols);
    
    for ii = 1:max(grp)
        f = grp == ii;
        if isempty(sum(f))
            continue
        end
        
        ofs = offset*(ii-1);
        y = repmat(ids(f)'+ofs,2,1);
        x = [starts(f); ends(f)];

        plot(x,y,'Color', cols(ii,:),'LineWidth',2)
        hold on
    end
    hold off
    
    pbaspect([1,2,1])
    xlim([fr1 fr2])
    ylim([0 max(ids)*1.04])
    set(gca,'YTickLabel',[],'TickLength',[0,0])
    
    %
    obj(1).plot_stimulation_time(params);
    
    % Save figure
    figure_dir = obj.figure_directory;
    formatspec = [repmat('%s_',1,numel(behaviours)-1),'%s'];
    fig_type = strcat("ethogram_",sprintf(formatspec,behaviours{:}));
%     fname = strcat(obj.get_full_genotype,".pdf");
    fname = strcat(obj.get_full_genotype);
    
    fig_path = fullfile(figure_dir,fig_type,fname);
    if ~isdir(fileparts(fig_path))
        mkdir(fileparts(fig_path));
    end
    
    % Define axis as output
    ax = gca;
    
end