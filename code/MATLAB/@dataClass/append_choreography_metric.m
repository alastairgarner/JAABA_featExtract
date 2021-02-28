%% append_choreography_metric

% Other m-files required: 
% Subfunctions: 
% MAT-files required: 
% See also:

% Author: Alastair Garner
% email: alastairgarner@outlook.com
% Website: https://github.com/alastairgarner/
% Sep 2019; Last revision: 

%%
function plotStructUpdated = append_choreography_metric(obj,plotStruct,feature,metric,frame,instance)
% 
% obj = dA;
% feature = 'curve';
% frame = [15 30];
% metric = 'mean';
% instance = 'all';

tmpStruct = plotStruct;
ts = obj.timeseries;

if isempty(ts) | ~any(strcmp(fieldnames(ts),feature))
    plotStructUpdated = tmpStruct;
end

if strcmp(feature,'speeddiff')
    diffs = arrayfun(@(x) [x.speed]-[x.crabspeed], ts, 'UniformOutput', false);
    [ts.speeddiff] = diffs{:};
end

if strcmp(instance,'normalise')
    p = obj.parse_protocol();
    frameBaseline = [5 p(1).start];
    anonFunc = @(x,y) x ./ mean(x([frameBaseline(1) <= y] & [frameBaseline(2) >= y]));
    valNormalised = cellfun(anonFunc,{ts.(feature)},{ts.et}, 'UniformOutput', false);
    [ts.(feature)] = valNormalised{:};
end

[etMin,etMax] = cellfun(@(x) bounds(x), {ts.et});
f1 = [frame(1) >= etMin] & [frame(2) <= etMax];
ts = ts(f1);

anonFunc = @(x,y) x([frame(1) <= y] & [frame(2) >= y]);
valCell = cellfun(anonFunc,{ts.(feature)},{ts.et},'UniformOutput', false)';        
etCell = cellfun(anonFunc,{ts.et},{ts.et},'UniformOutput', false)';        
switch feature
    case 'curve'
        peakProminence = 2;
    case 'crabspeed'
        peakProminence = 0.4;
end

switch metric
    case 'mean'
        val = cellfun(@(x) mean(x,'omitnan'), valCell)';
    case 'median'
        val = cellfun(@(x) median(x,'omitnan'), valCell)';
    case 'firstpeak'
        peakFunc = @(x) findpeaks( movmean([0 diff(x)],7),...
            'MinPeakProminence',peakProminence,'NPeaks',1);
        [~,loc] = cellfun(peakFunc,valCell,...
            'UniformOutput', false);
        val = cellfun(@(x,y) x(y), etCell, loc, 'UniformOutput', false)';
        filt = cellfun(@isempty,val);
        val(filt) = {nan};
        val = [val{:}];
        
end

if numel(val) == 0
    plotStructUpdated = tmpStruct;
    return
end

[tsFilt,I] = ismember([obj.uniID],[ts.uniID]);
[~,~,tsFilt] = unique(obj.timestamp_index(tsFilt),'stable');

tsDates = double([obj.date(tsFilt')]);

valName = strcat(feature,'_',metric,'_',string(frame(1)),'_',string(frame(2)));

if numel(fieldnames(tmpStruct))==0
    tmpStruct.name = valName;
    tmpStruct.id = [ts.id];
    tmpStruct.uniID = [ts.uniID];
    tmpStruct.y = [val];
    tmpStruct.dates = tsDates;
    tmpStruct.labels = obj.driver;
    tmpStruct.timestamp = tsFilt';
    tmpStruct.group = repelem(1,1,numel(val));
    
elseif ~any(ismember(string({tmpStruct.name}),valName))
    len = numel(tmpStruct)+1;
    tmpStruct(len).name = valName;
    tmpStruct(len).id = [ts.id];
    tmpStruct(len).uniID = [ts.uniID];
    tmpStruct(len).y = [val];
    tmpStruct(len).dates = tsDates;
    tmpStruct(len).labels = obj.driver;
    tmpStruct(len).timestamp = tsFilt';
    tmpStruct(len).group = repelem(1,1,numel(val));
    
else
    ix = ismember(string({tmpStruct.name}),valName);
    tmpStruct(ix).id = [tmpStruct(ix).id ts.id];
    tmpStruct(ix).uniID = [tmpStruct(ix).uniID ts.uniID];
    tmpStruct(ix).y = [tmpStruct(ix).y val];
    tmpStruct(ix).dates = [tmpStruct(ix).dates tsDates];
    tmpStruct(ix).labels = [tmpStruct(ix).labels obj.driver];
    tmpStruct(ix).timestamp = [tmpStruct(ix).timestamp tsFilt'];
    tmpStruct(ix).group = [tmpStruct(ix).group ,...
                            repelem(max([tmpStruct(ix).group,0])+1,1,numel(val))];
    
end
    
plotStructUpdated = tmpStruct;

end


    