%% merge_dataClasses

% Other m-files required: 
% Subfunctions: 
% MAT-files required: 
% See also:

% Author: Alastair Garner
% email: alastairgarner@outlook.com
% Website: https://github.com/alastairgarner/
% Sep 2019; Last revision: 

%%
function one_dataClass = merge_dataClasses(obj)
    dca = obj;
    % merge identical behaviours across timestamps
    pipelines = [dca.pipeline];
    [pipes,~,idx] = unique(pipelines);
    dca_new = dataClass;
    for ii = 1:numel(pipes)
        clear temp
        f = idx == ii;

        fnames = fieldnames(dca(f));
        for jj = 1:numel(fnames)
            if isstring(dca(find(f,1)).(fnames{jj}))
                dca_new(ii,1).(fnames{jj}) = vertcat(dca(f).(fnames{jj}))';
            elseif isnumeric(dca(find(f,1)).(fnames{jj}))
                dca_new(ii,1).(fnames{jj}) = [dca(f).(fnames{jj})];
            end
        end

        % merge choreography data
        if strcmp(pipes{ii},'choreography')
            f = find(f); % EDIT 20200302
            ff = cellfun(@(x) numel(fieldnames(x)),{dca(f).timeseries}) ~= 0; % EDIT 20200302
            temp = vertcat(dca(f(ff)).timeseries);
            dca_new(ii,1).timeseries = temp;
        % merge salam/jaaba/jb data
        elseif any(strcmp(pipes{ii},{'salam','jaaba','jb'}))
            f = find(f); % EDIT 20200302
            ff = cellfun(@(x) numel(fieldnames(x)),{dca(f).behaviour}) ~= 0; % EDIT 20200302
%             behs = vertcat(dca(f).behaviour); % EDIT 20200302
            behs = vertcat(dca(f(ff)).behaviour);
            if isempty(behs)
                continue
            end
            
            [C,~,ic] = unique({behs.behaviour});
            % loop through behaviours
            for jj = 1:numel(C)
                f2 = ic == jj;
                fnames = fieldnames(behs(f2));
                fields = cellfun(@(x) [behs(f2).(x)], fnames, 'UniformOutput', false);
                fields(1) = unique({behs(f2).behaviour});
                temp(jj,1) = cell2struct(fields,fnames);
            end
            dca_new(ii,1).behaviour = temp;
        end 
    end

    dca = dca_new;
    %%
    [C,ia,ic] = unique(dca.get_full_timestamp);

    fdnames = fieldnames(dca);
    temp = dataClass();
    for jj = 1:numel(fdnames)
        if any(strcmp(fdnames{jj},{'aniID','timestamp_index'}))
            [C,ia,ic] = unique([dca.uniID],'stable');
            ids = [dca.(fdnames{jj})];
            temp.(fdnames{jj}) = ids(ia);

        elseif strcmp(fdnames{jj},'uniID')
            temp.(fdnames{jj}) = unique([dca.uniID],'stable');
            
        elseif any(strcmp(fdnames{jj},["date","time"]))
            dts = unique([[dca.date];[dca.time]]','rows','stable');
            temp.date = dts(:,1)';
            temp.time = dts(:,2)';

        elseif strcmp(fdnames{jj},'track_start')
            [~,~,idx] = unique([dca.uniID],'stable');
            temp.(fdnames{jj}) = accumarray(idx,[dca.(fdnames{jj})]',[],@min)';

        elseif strcmp(fdnames{jj},'track_end')
            [~,~,idx] = unique([dca.uniID],'stable');
            temp.(fdnames{jj}) = accumarray(idx,[dca.(fdnames{jj})]',[],@max)';

        elseif strcmp(fdnames{jj},'behaviour')
            cells = {dca.(fdnames{jj})};
            goodstructs = cellfun(@(x) numel(fieldnames(x))>0, cells);
            temp.(fdnames{jj}) = vertcat(cells{goodstructs});

        elseif strcmp(fdnames{jj},'timeseries')
            notjaaba = ~cellfun(@(x) all(contains(x,'jaaba')), {dca.pipeline});
            dca_gd = dca(notjaaba);
            fnames = []; ids = [];
            for xx = 1:numel(dca_gd)
                fnames = vertcat(fnames, fieldnames(dca_gd(xx).(fdnames{jj})));
            end
            fnames = unique(fnames,'stable');
            un_ids = unique([dca_gd.uniID]);
            nobjs = length(un_ids);
            f_mat = [fnames, repelem({[]}, numel(fnames),1)]';
            templat = struct(f_mat{:});
            templat = repelem(templat,nobjs,1);

            for xx = 1:numel(dca_gd)
                try
                    logica = ismember(un_ids,[dca_gd(xx).timeseries.uniID]);
                catch
                    continue
                end
                id_cell = num2cell(un_ids);
                [templat.id] = id_cell{:};
                for yy = 1:numel(fnames)
                    [templat(logica).(fnames{yy})] = dca_gd(xx).timeseries.(fnames{yy});
                end
            end
            temp.(fdnames{jj}) = templat;

        elseif strcmp(fdnames{jj},'raw_data')
            temp.(fdnames{jj}) = [];
        else
            temp.(fdnames{jj}) = unique([dca.(fdnames{jj})]);
        end
    end

    one_dataClass = temp;

%% Old Version
%     dca = obj;
%     
%     [C,ia,ic] = unique(dca.get_full_timestamp);
% 
%     fdnames = fieldnames(dca);
%     temp = dataClass();
%     for jj = 1:numel(fdnames)
%         if strcmp(fdnames{jj},'aniID')
%             temp.(fdnames{jj}) = unique([dca.(fdnames{jj})],'stable');
% 
%         elseif strcmp(fdnames{jj},'track_start')
%             [~,~,idx] = unique([dca.aniID],'stable');
%             temp.(fdnames{jj}) = accumarray(idx,[dca.(fdnames{jj})]',[],@min)';
% 
%         elseif strcmp(fdnames{jj},'track_end')
%             [~,~,idx] = unique([dca.aniID],'stable');
%             temp.(fdnames{jj}) = accumarray(idx,[dca.(fdnames{jj})]',[],@max)';
% 
%         elseif strcmp(fdnames{jj},'behaviour')
%             cells = {dca.(fdnames{jj})};
%             goodstructs = cellfun(@(x) numel(fieldnames(x))>0, cells);
%             temp.(fdnames{jj}) = vertcat(cells{goodstructs});
% 
%         elseif strcmp(fdnames{jj},'timeseries')
%             notjaaba = [dca.pipeline]' ~= 'jaaba';
%             dca_gd = dca(notjaaba);
%             fnames = []; ids = [];
%             for xx = 1:numel(dca_gd)
%                 fnames = vertcat(fnames, fieldnames(dca_gd(xx).(fdnames{jj})));
%             end
%             fnames = unique(fnames,'stable');
%             un_ids = unique([dca_gd.aniID]);
%             nobjs = length(un_ids);
%             f_mat = [fnames, repelem({[]}, numel(fnames),1)]';
%             templat = struct(f_mat{:});
%             templat = repelem(templat,nobjs,1);
% 
%             for xx = 1:numel(dca_gd)
%                 try
%                     logica = ismember(un_ids,[dca_gd(xx).timeseries.id]);
%                 catch
%                     continue
%                 end
%                 id_cell = num2cell(un_ids);
%                 [templat.id] = id_cell{:};
%                 for yy = 1:numel(fnames)
%                     [templat(logica).(fnames{yy})] = dca_gd(xx).timeseries.(fnames{yy});
%                 end
%             end
%             temp.(fdnames{jj}) = templat;
% 
%         elseif strcmp(fdnames{jj},'raw_data')
%             temp.(fdnames{jj}) = [];
%         else
%             temp.(fdnames{jj}) = unique([dca.(fdnames{jj})]);
%         end
%     end
%     
%     one_dataClass = temp;
end
