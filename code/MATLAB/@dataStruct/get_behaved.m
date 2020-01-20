%% get_behaved

% Other m-files required: 
% Subfunctions: 
% MAT-files required: 
% See also:

% Author: Alastair Garner
% email: alastairgarner@outlook.com
% Website: https://github.com/alastairgarner/
% August 2019; Last revision:

function [behaved,uniIDs] = get_behaved(obj,behaviour)

    temp = obj.get_behaviour_data(behaviour);

    bins = 0:.2:max(ceil(temp.tEnd));
    flt = ~isnan(temp.tStart);

    aniIDs = temp.aniID(flt);
    uniIDs = 1:length(aniIDs);

    temp = temp.(behaviour);
    counts = cellfun(@length,temp.bStart);
    uniIDs = repelem(uniIDs,counts(flt));
    bStart = [temp.bStart{flt}];
    bEnd = [temp.bEnd{flt}];

    flt = ~isnan(bStart);    
    uniIDs = uniIDs(flt);
    [binStart,~] = discretize(bStart(flt),bins);
    [binEnd,~] = discretize(bEnd(flt),bins);

    spStart = sparse(uniIDs',binStart',[1],max(uniIDs),length(bins));
    spEnd = sparse(uniIDs',binEnd',[-1],max(uniIDs),length(bins));
    spBoth = cumsum((spStart+spEnd),2);
    behaved = full(sum(spBoth,1));
end