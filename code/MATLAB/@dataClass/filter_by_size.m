%% fix_axes

% Other m-files required: 
% Subfunctions: 
% MAT-files required: 
% See also:

% Author: Alastair Garner
% email: alastairgarner@outlook.com
% Website: https://github.com/alastairgarner/
% Sep 2019; Last revision: 

%%
function obj = filter_by_size(obj,area_threshold)
    for xx = 1:numel(obj)
        if isempty(obj(xx).timeseries) | isempty(obj(xx).behaviour) | ~isfield(obj(xx).timeseries,'et')
            fprintf('Timestamp lacking data... skipping\n')
            continue
        end
        f = cellfun(@(x) mean(x)>=area_threshold, {obj(xx).timeseries.area});
        gd_ids = [obj(xx).timeseries(f).uniID];

        obj(xx).timeseries = obj(xx).timeseries(f);

        for yy = 1:numel(obj(xx).behaviour)
            filt = ismember(obj(xx).behaviour(yy).uniID,gd_ids);
            fnames = fieldnames(obj(xx).behaviour(yy));
            temp = cellfun(@(x) obj(xx).behaviour(yy).(x)(filt), fnames,...
                'UniformOutput', false, 'ErrorHandler', @(a,b) []);
            temp{1} = obj(xx).behaviour(yy).behaviour;
            fields = [fnames,temp]';
        %                 temp = struct(fields{:})
            obj(xx).behaviour(yy) = struct(fields{:});
        end
    end
end