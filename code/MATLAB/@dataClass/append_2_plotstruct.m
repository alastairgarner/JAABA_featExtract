%% append_2_plotstruct

% Other m-files required: 
% Subfunctions: 
% MAT-files required: 
% See also:

% Author: Alastair Garner
% email: alastairgarner@outlook.com
% Website: https://github.com/alastairgarner/
% Sep 2019; Last revision: 

%%
function plot_struct = append_2_plotstruct(obj,plot_struct,params,behaviour,metric,instance,method,frame)
%     obj = dA;
%     
%     behaviour = 'rolls';
%     instance = 'all';
%     metric = 'duration';

    if numel(fieldnames(plot_struct)) == 0
        fnames = {'x','y','c','group','timestamp','dates','labels'};
        fnames = [fnames;repelem({[]},1,numel(fnames))];
        plot_struct = struct(fnames{:});
        group = 0;
    else
        group = max(plot_struct.group);
    end
    
    x = []; y = []; c = []; labels = []; timestamps = []; dates =[];
    for xx = 1:numel(obj)
    	f = strcmp({obj(xx).behaviour.behaviour},behaviour);
        if any(f)
            temp = obj(xx).behaviour(f);
            
            f3 = temp.start > frame(1) & temp.start < frame(2);
        elseif strcmp('area',behaviour)
            metric = 'amplitude';
            temp = obj(xx).behaviour(end);
            temp = structfun(@(x) [], temp, 'UniformOutput', false);
            temp.behaviour = 'area';
            [C,ia,ic] = unique([obj(xx).behaviour.uniID]);
            temp.uniID = C;
            tstart = [obj(xx).behaviour.track_start];
            tend = [obj(xx).behaviour.track_end];
            temp.track_start = tstart(ia);
            temp.track_end = tend(ia);
            temp.start = tstart(ia);
            temp.end = tend(ia);
            
            gd_id = ismember([obj(xx).timeseries.uniID],C);
            temp.amplitude = cellfun(@mean,{obj(xx).timeseries(gd_id).area});
            
            f3 = repelem(true,1,numel(temp.uniID));
        end
        [~,id_filt] = ismember(temp.uniID,obj(xx).uniID);
        temp.timestamp = obj(xx).timestamp_index(id_filt);
        un_dates = double(obj(xx).date(unique(temp.timestamp)));
        [~,~,temp.timestamp] = unique(temp.timestamp,'stable');
        
        %
        f1 = ~isnan(temp.start);
        f2 = temp.track_start < frame(1) & temp.track_end > frame(2);
        filt = f1 & f2 & f3;
        
        timestamp_n = accumarray(temp.timestamp(f2),temp.uniID(f2)',...
            [max(temp.timestamp(f2)) 1],@(x) numel(unique(x)));
%         un_dates = double(obj(xx).date(unique(temp.timestamp(f2))));
        
        switch instance
            case 'all'
                
            case 'first'
                ref = [temp.uniID;temp.start]';
                [dds,~,uds] = unique(temp.uniID,'stable');
                cnts = accumarray(uds,1);
                
                func = @(x) min([x(x>frame(1) & x<frame(2));nan]);
                val = accumarray(uds,temp.start',[],func);
                nf = ismember(ref,[dds',val],'rows');
                filt = filt & nf';
            case 'last'
                ref = [temp.uniID;temp.start]';
                [dds,~,uds] = unique(temp.uniID,'stable');
                cnts = accumarray(uds,1);
                
                func = @(x) max([x(x>frame(1) & x<p.frame(2));nan]);
                val = accumarray(uds,temp.start',[],func);
                nf = ismember(ref,[dds',val],'rows');
                filt = filt & nf';
        end
        
        switch metric
            case 'duration'
                val = temp.end(filt) - temp.start(filt);
            case 'start'
                val = temp.start(filt);
            case 'end'
                val = temp.end(filt);
            case 'amplitude'
                val = temp.amplitude(filt);            
        end
        
        [~,~,id] = unique(temp.uniID(filt),'stable');
%         [~,~,ts] = unique(temp.timestamp(filt),'stable');
        ts = temp.timestamp(filt);
        
        % metric by animal
        ts_i = accumarray(id,ts',[],@mean);
        switch method
            case 'mean'
                tot = accumarray(id,val',[],@mean);
            case 'sum'
                tot = accumarray(id,val',[],@sum);
            case 'count'
                tot = accumarray(id,1);
            case 'proportion'
                counts = accumarray(ts_i,1,[numel(timestamp_n),1]);
                tot = counts./timestamp_n;
                tot = repelem(tot',1,counts)';                
        end
                
        % metric by timestamp
        ave_tot = accumarray(ts_i,tot,[numel(timestamp_n),1],@mean,NaN)';
        grp = repelem([xx+group],1,numel(ave_tot));
        
        x = [x, grp];
        y = [y, ave_tot];
        labels = [labels obj(xx).driver];
        timestamps = [timestamps 1:numel(timestamp_n)];
        dates = [dates un_dates];
    end
        
%     plot_struct.x = [plot_struct.x, x];
    plot_struct.y = [plot_struct.y, y];
    plot_struct.group = [plot_struct.group, x];
    plot_struct.labels = [plot_struct.labels, labels];
    plot_struct.timestamp = [plot_struct.timestamp, timestamps];
    plot_struct.dates = [plot_struct.dates, dates];
    
end