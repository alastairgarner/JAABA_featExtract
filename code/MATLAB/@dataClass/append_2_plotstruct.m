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
function tmpStruct = append_2_plotstruct(obj,plotStruct,params,behaviour,metric,instance,method,frame)
%     
%     behaviour = 'rolls';
%     instance = 'all';
%     metric = 'duration';

    tmpStruct = plotStruct;
    
    valName = strcat(behaviour,'_',metric,'_',method,'_',string(frame(1)),'_',string(frame(2)));
    if numel(fieldnames(tmpStruct))~=0
        mtch = ismember(string({tmpStruct.name}),valName);
       if any(mtch)
           group = max([tmpStruct(mtch).group]);
       else
           group = max([tmpStruct.group])-1;
       end
    else
        group = 0;
    end
    
    %%
    %%%%
    x = []; y = []; c = []; labels = []; timestamps = []; dates = [];
    for xx = 1:numel(obj)
    	f = strcmp({obj(xx).behaviour.behaviour},behaviour);
        if any(f)
            temp = obj(xx).behaviour(f);
            
            if strcmp(instance,'normalise')
                prot = obj.parse_protocol();
                fr = [5 prot(1).start];

                filt = temp.start > fr(1) & temp.start < fr(2);
                [ids,ia,ic] = unique([temp.uniID]);

                baseline = accumarray(ic(filt),temp.(metric)(filt)',[max(ic),1],@(x) mean(x,'omitnan'));
                baseline = baseline(ic)';
                baseline(baseline==0) = nan;

                temp.(metric) = temp.(metric)./baseline;
            end
            
            f3 = temp.start > frame(1) & temp.start < frame(2);
            f3 = f3 | isnan(temp.start);
%             f3 = repelem(true,1,numel(temp.id));
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
        filt = f2 & f3;
        
        if ~any(f2)
            return
        end
        
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
            case 'frequency'
                val = temp.frequency(filt);  
        end
        
        [~,~,id] = unique(temp.uniID(filt),'stable');
%         [~,~,ts] = unique(temp.timestamp(filt),'stable');
        ts = temp.timestamp(filt);
        
        % metric by animal
        ts_i = accumarray(id,ts',[],@mean);
        switch method
            case 'mean'
%                 tot = accumarray(id,val',[],@mean);
                tot = accumarray(id,val',[],@(x) mean(x,'omitnan'));
            case 'sum'
%                 tot = accumarray(id,val',[],@sum);
                tot = accumarray(id,val',[],@(x) sum(x,'omitnan'));
                tot(tot==0) = nan;
            case 'count'
%                 tot = accumarray(id,1);
                tot = accumarray(id,val',[],@(x) sum(~isnan(x)));
                tot(tot==0) = nan;
            case 'proportion'
%                 counts = accumarray(ts_i,1,[numel(timestamp_n),1]);
%                 tot = counts./timestamp_n;
%                 tot = repelem(tot',1,counts)';
                tot = accumarray(id,val',[],@(x) any(~isnan(x)));
        end
        
        % metric by timestamp
        if false
            ave_tot = accumarray(ts_i,tot,[numel(timestamp_n),1],@mean,NaN)';
            grp = repelem([xx+group],1,numel(ave_tot));
            tstamps = 1:numel(timestamp_n);
        else
            ave_tot = tot';
%             mtch = strcmp(tmpStruct.labels,obj.driver);
%             if plot_struct.group(mtch)
%                 repelem(mean(plot_struct.group(mtch)),1,numel(ave_tot))
%             else
            grp = repelem([xx+group],1,numel(ave_tot));
%             end
            tstamps = ts_i';
            un_dates = un_dates(ts_i');
        end
        
        x = [x, grp];
        y = [y, ave_tot];
        labels = [labels obj(xx).driver];
        timestamps = [timestamps tstamps];
        dates = [dates un_dates];
    end
    
    
    %%
    valName = strcat(behaviour,'_',metric,'_',method,'_',string(frame(1)),'_',string(frame(2)));
    
    if numel(fieldnames(tmpStruct))==0
        tmpStruct.name = valName;
        tmpStruct.y = y;
        tmpStruct.group = x;
        tmpStruct.labels = labels;
        tmpStruct.timestamp = timestamps;
        tmpStruct.dates = dates;

    elseif ~any(ismember(string({tmpStruct.name}),valName))
        len = numel(tmpStruct)+1;
        tmpStruct(len).name = valName;
        tmpStruct(len).y = y;
        tmpStruct(len).group = x;
        tmpStruct(len).labels = labels;
        tmpStruct(len).timestamp = timestamps;
        tmpStruct(len).dates = dates;

    else
        ix = ismember(string({tmpStruct.name}),valName);
        tmpStruct(ix).y = [tmpStruct(ix).y, y];
        tmpStruct(ix).group = [tmpStruct(ix).group, x];
        tmpStruct(ix).labels = [tmpStruct(ix).labels, labels];
        tmpStruct(ix).timestamp = [tmpStruct(ix).timestamp, timestamps];
        tmpStruct(ix).dates = [tmpStruct(ix).dates, dates];

    end
    
    if strcmp(method,'proportion')
        try
            tmpStruct(len).y = logical(tmpStruct(len).y);
        catch
            tmpStruct(ix).y = logical(tmpStruct(ix).y);
    end
    
end