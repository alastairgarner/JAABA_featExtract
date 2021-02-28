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
function plot_stimulation_time(obj,params)

    x_lim = get(gca,'XLim');
    y_lim = get(gca,'YLim');
    asp = get(gca,'PlotBoxAspectRatio');

    col = {[255 79 0];[170 220 215]};
    col = {[225 114 98];[170 220 215]};

    fnames = fieldnames(obj);
    fnames = fnames(startsWith(fnames,"protocol"));
    protocols = cellfun(@(x) obj.(x), fnames);

    exp = '[_](?<sta>\d+)s(?<n>\d+)x(?<len>\d+)s(?<int>\d+)';
    prot = regexp(protocols,exp,'names');
    prot = vertcat(prot{:});

    p = arrayfun(@(x) structfun(@double,x,'UniformOutput',false), prot);

    func = @(x) [x.sta x.sta+x.len] + ([[1:x.n]-1] * (x.int+x.len))';
    bounds = arrayfun(func, p, 'UniformOutput', false);
    bounds = vertcat(bounds{:});
    x = repelem(bounds,1,2);
    % stim bar - color
    stimbar_height = diff(y_lim)/diff(x_lim)*[asp(1)/asp(2)]*0.2;
    ymax = max(y_lim);
    y = [ymax*.98] + [-stimbar_height 0 0 -stimbar_height];
    y = repmat(y,size(bounds,1),1);
    patch(x',y',col{1}/255, 'EdgeColor', 'none','Tag','stimbox');
    % stim box - grey
%         y = repmat([y_lim,fliplr(y_lim)],size(bounds,1),1);        
%         patch(x',y',[.9 .9 .9], 'EdgeColor', 'none','Tag','stimbox');

    children = get(gca,"Children");
    f = arrayfun(@(x) strcmp(x.Tag,"stimbox"), children);
    % send to bottom
    set(gca,"Children",[children(~f);children(f)]);
    % bring to top
%     set(gca,"Children",[children(f);children(~f)]);
    set(gca,'XLim',x_lim,'YLim',y_lim);

end