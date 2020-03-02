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
function [p,h,stats] = mult_comp_wilcox_byday(plot_struct)
    val = plot_struct;

    counts = accumarray(val.group',1);
    means = accumarray(val.group',val.y',[],@(x) mean(x,'omitnan'));
    cells = mat2cell(val.y,1,counts);
    dates = mat2cell(val.dates,1,counts);
     
    v = strcmpi(val.labels,'attp2');
    means = means(~v);
    labels = val.labels(~v);
    
    f = val.group == find(v);
    [C,ia,ic] = unique(val.dates,'stable');
    Y = accumarray(ic,val.group',[],@(x) numel(unique(x)));
    comps = Y(ic)-1;
%     comps(comps==0) = 1;
    
    clear p h stats
    for ii = 1:numel(cells)
        overlap = ismember(dates{v},dates{ii});
        if ~all(isnan(cells{ii})) & ~all(isnan([cells{v}(overlap), nan]))
    %         [p,h,stats] = cellfun(@(x) ranksum(cells{v},x),cells);
%             [p(ii),h(ii),stats(ii)] = ranksum(cells{v},cells{ii});
%             overlap = ismember(dates{v},dates{ii});
            [p(ii),h(ii),stats(ii)] = ranksum(cells{v}(overlap),cells{ii},'method','approximate');
            cval(ii) = max(comps(overlap));
        else
            p(ii) = 1;
            h(ii) = false;
            stats(ii) = struct('zval',NaN,'ranksum',NaN);
            cval(ii) = NaN;
        end
    end
    
    p = p.*cval;
    p(p>1) = 1;

end
