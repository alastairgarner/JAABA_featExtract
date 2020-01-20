%% get_tracked

% Other m-files required: 
% Subfunctions: 
% MAT-files required: 
% See also:

% Author: Alastair Garner
% email: alastairgarner@outlook.com
% Website: https://github.com/alastairgarner/
% August 2019; Last revision:

function [tracked,uniIDs] = get_tracked(obj,behaviour)
    
    temp = obj.get_behaviour_data(behaviour);

    bins = 0:.2:max(ceil(temp.tEnd));
    flt = ~isnan(temp.tStart);

    aniIDs = temp.aniID(flt);
    uniIDs = 1:length(aniIDs);
    [binStart,~] = discretize(temp.tStart(flt),bins);
    [binEnd,~] = discretize(temp.tEnd(flt),bins);

    spStart = sparse(uniIDs',binStart',[1],max(uniIDs),length(bins));
    spEnd = sparse(uniIDs',binEnd',[-1],max(uniIDs),length(bins));
    spBoth = cumsum((spStart+spEnd),2);
    tracked = full(sum(spBoth,1));
end