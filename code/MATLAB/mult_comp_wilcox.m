%% mult_comp_wilcox.m

% Other m-files required: 
% Subfunctions: 
% MAT-files required: 
% See also:

% Author: Alastair Garner
% email: alastairgarner@outlook.com
% Website: https://github.com/alastairgarner/
% Sep 2019; Last revision: 

%%
function [p,h,stats] = mult_comp_wilcox(plot_struct)
    val = plot_struct;

    counts = accumarray(val.group',1);
    means = accumarray(val.group',val.y',[],@(x) mean(x,'omitnan'));
    cells = mat2cell(val.y,1,counts);
    v = strcmpi(val.labels,'attp2');
    means = means(~v);
    labels = val.labels(~v);

    clear p h stats
    for ii = 1:numel(cells)
        if ~all(isnan(cells{ii}))
    %         [p,h,stats] = cellfun(@(x) ranksum(cells{v},x),cells);
            [p(ii),h(ii),stats(ii)] = ranksum(cells{v},cells{ii});
        else
            p(ii) = 1;
            h(ii) = false;
            stats(ii) = struct('zval',NaN,'ranksum',NaN);
        end
    end
    
    p = p.*sum(~v);
    p(p>1) = 1;

end
