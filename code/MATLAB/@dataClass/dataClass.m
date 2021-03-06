%% dataContainer

% Other m-files required: 
% Subfunctions: 
% MAT-files required: 
% See also:

% Author: Alastair Garner
% email: alastairgarner@outlook.com
% Website: https://github.com/alastairgarner/
% Sep 2019; Last revision: 


%% TO DO 
%
% - eliminate duplicate entried in salam compile
%
%
%
%
%



%%
classdef dataClass
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        filepath = ""
        pipeline = ""
        behaviourtype = ""
        
        date = ""
        time = ""
        driver = ""
        effector = ""
        protocol1 = ""
        protocol2 = ""
        protocol3 = ""
        protocol4 = ""
        rig = ""
        
        aniID = []
        uniID = []
        timestamp_index = []
        track_start = [];
        track_end = [];

        behaviour = struct();
        timeseries = struct();
%         object = struct();
        
        raw_data = [];
        
        group;
        group_driver;
        group_effector;
        group_effector_protocol;
        group_genotype;
        group_protocol;
        group_experiment;
        
        figure_directory = [];
    end
    
    %% METHODS NORMAL
    methods
        %% dataClass Constructor
        function obj = dataClass(filepaths)
            if nargin ~= 0
                
                timestamps = regexp(filepaths,'\d{8}_\d{6}','match','once');
                filt = ~cellfun(@isempty,timestamps);
                timestamps = string(timestamps(filt));
                filepaths = filepaths(filt);
                
                [split_filepaths,names,exts] = cellfun(@(x) fileparts(x),filepaths,'UniformOutput', false);
                identifiers = strcat(split_filepaths,'@',timestamps);

                [unique_identifiers,~,indicies] = unique(identifiers);
                elems = numel(unique_identifiers);
                obj(elems,1) = dataClass();

                for ii = 1:numel(unique_identifiers)
                    idx = indicies == ii;
                    obj(ii).filepath = string(filepaths(idx))';
                    obj(ii) = obj(ii).update_pipeline();
                    temp = dataClass.parse_filepaths(filepaths(idx));
                    % merge struct to object
                    for fn = fieldnames(temp)'    %enumerat fields
                        try
                          obj(ii).(fn{1}) = temp.(fn{1});   %and copy
                        catch
                          warning('Could not copy field %s', fn{1});
                        end
                    end
                end
                obj = update_groups(obj);
            end
        end
        
        %% Save Data
        function savenames = save_data(obj,params)
            savenames = [];
            for ii = 1:numel(obj)
%                 dC = struct(obj(ii));
                dC = obj(ii);
                savename = fullfile(params.directories.data_compiled,obj(ii).pipeline,strcat(obj(ii).get_full_experiment, '.mat'));
                if ~isdir(fileparts(savename))
                    mkdir(fileparts(savename))
                end
                save(savename,'dC');
                savenames = [savenames, savename];
            end
        end
        
        %% update pipeline
        function obj = update_pipeline(obj)
            for ii = 1:numel(obj.filepath)
                [fold,name,ext] = fileparts(obj.filepath(ii));
                if startsWith(ext,".blo")
                    obj.pipeline = "mwt";
                    return
                elseif startsWith(ext,".txt")
                    if contains(lower(name),"chore")
                        obj.pipeline = "choreography";
                        return
                    elseif contains(lower(name),"animal_stats")
                        obj.pipeline = "salam";
                        return
                    end
                elseif startsWith(ext,".mat")
                    if contains(name,"trx")
                        file_check = dir(fold);
                        if any({file_check.name} == "perframe")
                            obj.pipeline = "jaaba";
                            return
                        else
                            obj.pipeline = "jb";
                            return
                        end
                    elseif contains(fold,"tempdir")
                        [~,typ,~] = fileparts(fold);
                        obj.pipeline = typ;
                        return
                    end
                end
            end
        end
        
        %% Load data
        function obj = load_data(obj)
            if contains(obj.filepath,["tempdir","data_compiled"])
                obj = obj.load_compileddata();
            elseif strcmp(obj.pipeline,"mwt")
                obj = obj.load_blobsdata();
            elseif strcmp(obj.pipeline,"choreography")
                obj = obj.load_choredata();
            elseif strcmp(obj.pipeline,"salam")
                obj = obj.load_salamdata();
            elseif strcmp(obj.pipeline,"jaaba")
                obj = obj.load_jaabadata();
            elseif strcmp(obj.pipeline,"jb")
                obj = obj.load_jbdata();
            end
        end
        
        %% Load blobs data
        function obj = load_compileddata(obj)
            temp = load(obj.filepath(1));
            fields_copy = {'filepath','group','group_driver','group_effector','group_effector_protocol','group_genotype',...
                'group_protocol','group_experiment','figure_directory'};
            for ii = 1:numel(fields_copy)
                temp.dC.(fields_copy{ii}) = obj.(fields_copy{ii});
            end
%             temp.dC.filepath = obj.filepath;
            obj = temp.dC;
        end
            
        %% Load blobs data
        function obj = load_blobsdata(obj)
            filepaths = obj.filepath;
            [fpath,~,~] = fileparts(filepaths(1));
            d = dir(fullfile(fpath,'*.summary'));
            if ~size(d,1)
                return
            end

            delimiter = ' ';
            formatSpec = '%f%f%f%f%f%f%f%f%f%f%s%f%f%f%s%[^\n\r]';
            columns = {'aniID','frame','time',...
                            'x_cent','y_cent',...
                            'n_px','x_vect',...
                            'y_vect','std_orthog',...
                            'len_px','wid_px',...
                            'x_cont','y_cont',...
                            'npts'};

            datfull = []; outlinesfull = [];
            for ii = 1:numel(filepaths)
                fileID = fopen(filepaths(ii));
                datA = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'CollectOutput', 1, 'ReturnOnError', false, 'treatas', '%');
                fclose(fileID);

                datmat = [datA{[1,3]}];
                filt = isnan(datmat(:,1));
                aniID = datmat(filt,2);

                datmat = datmat(~filt,:);
                outlines = datA{4}(~filt,1);
                diffs = diff(find([filt' 1]))-1;
                aniIDs = repelem(aniID,diffs,1);
                datmat = [aniIDs datmat];

                datfull = vertcat(datfull,datmat);
                outlinesfull = [outlinesfull; outlines];
                clear datmat
            end

            [datfull,idx] = sortrows(datfull,[1,2]);
            outlinesfull = outlinesfull(idx);
            datcell = mat2cell(datfull',ones(size(datfull,2),1),size(datfull,1));
            obj.raw_data = cell2struct(datcell,columns);
            obj.raw_data.outline = outlinesfull';
            obj.pipeline = "mwt";
            obj.behaviourtype = "";

            fprintf('\t----%s - blobs loaded \n',obj.get_full_experiment)
        end
        
        %% Load chore data
        function obj = load_choredata(obj)
            if ~check_file(obj)
                return
            end
            filepath = obj.filepath;
            fileID = fopen(filepath(1),'r');
            fline = fgetl(fileID);

            chores = split(fline,' ');
            chname = vertcat('et',chores);
            chores = chores([4:end]);
            len = length(chores);

            delimiter = ' ';
            startRow = 1;

            frewind(fileID)
            formatSpec = ['%f%s' repmat('%f',1,len+1) '%[^\n\r]'];
            dataArray = textscan(fileID, formatSpec,...
                'Delimiter', delimiter,...
                'MultipleDelimsAsOne', true,...
                'HeaderLines' ,startRow,...
                'ReturnOnError', false);
            fclose(fileID);

            dataArray = dataArray(:,[1,3:end-1]);

            obj.raw_data.id = dataArray{1};
            obj.raw_data.et = dataArray{2};
            for jj = 1:length(dataArray)-2
                field_name = chores{jj};
                obj.raw_data.(field_name) = dataArray{jj+2};
            end
            obj.pipeline = "choreography";
            obj.behaviourtype = "";
            
            fprintf('\t----%s - chore loaded \n',obj.get_full_experiment)
        end
        
        %% Load feature extraction (salam) data
        function obj = load_salamdata(obj)
            delimiter = '\t';
            startRow = 1;
            formatSpec = '%s%f%f%f%f%f%f%f%f%f%f%[^\n\r]';

            filenames = obj.filepath;

            expr = '[\_]([A-Za-z0-9]+).txt';
            behaviour = regexp(filenames,expr,'tokens','once');
            if iscell(behaviour)
                behaviour = [behaviour{:}];
            end
            
            for ii = 1:length(filenames)
                fileID = fopen(filenames{ii},'r');
                dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow, 'ReturnOnError', true,'treatAsEmpty', {'NA'});
                fclose(fileID);
                
                if numel(unique(cellfun(@numel, dataArray))) > 1
                    continue
                end

                % Restructure crawl data to match other data
                filler = zeros(length(dataArray{1}),1);
                if ~all(isnan(dataArray{end-1}))
                    dataArray = dataArray([1:5,9,6,10,11,8,7]);
                end

                salam_matrix = [dataArray{[2:11]}];
                obj.raw_data.(behaviour(ii)) = salam_matrix;
            end
            obj.pipeline = "salam";
            obj.behaviourtype = "";
            fprintf('\t----%s - salam loaded \n',obj.get_full_experiment)
        end
        
        %% load jaaba data
        function obj = load_jaabadata(obj)
            for ii = 1:numel(obj.filepath)
                filename = obj.filepath(ii);
                data = load(filename);
                f = fieldnames(data);
                for jj = 1:length(f)
                   obj.raw_data.(f{jj}) = data.(f{jj});
                end
            end
            fprintf('\t----%s - jaaba loaded \n',obj.get_full_experiment)
        end
        
        %% load jb data
        function obj = load_jbdata(obj)
            for ii = 1:numel(obj.filepath)
                data = load(obj.filepath(ii));
                f = fieldnames(data);
                for jj = 1:length(f)
                   obj.raw_data.(f{jj}) = data.(f{jj});
                end
            end
            fprintf('\t----%s - jb loaded \n',obj.get_full_experiment)
        end
        
        %% Compile data
        function obj = compile_data(obj)
            if isempty(obj.raw_data)
                return
            end
            switch obj.pipeline
                case "mwt"
                    obj = obj.compile_mwt();
                case "choreography"
                    obj = obj.compile_choreography();
                case "salam"
                    obj = obj.compile_salam();
                case "jaaba"
                    obj = obj.compile_jaaba();
                case "jb"
                    obj = obj.compile_jb();
                otherwise
                    fprintf("\n unsupported pipeline specified")
            end
        end
        
        
        %% compile salam data
        function obj = compile_salam(obj)
            if isempty(obj.raw_data)
                fprintf('\t\t ...no salam data currently loaded\n')
                data = [];
                return
            end

            fnames = fieldnames(obj.raw_data);
            [rows,cols] = structfun(@size, obj.raw_data);
            behID = repelem([1:length(fnames)]',rows,1);
            dat = struct2cell(obj.raw_data);                
            dat = [vertcat(dat{:}) behID];
            dat = sortrows(dat, [1,2,7]);

            aniID = dat(:,1);
            uniID = unique(aniID);
            tStart = []; tEnd = []; bStart = []; bEnd = []; bType = []; bAmp = []; bFreq = [];

            for jj = 1:size(uniID,1)
                filt = aniID == uniID(jj);
                temp = dat(filt,[2,3,8,11]);
                temp = sortrows(temp,1);
                bS = dat(filt,7);
                bE = dat(filt,7) + dat(filt,5);
                bA = dat(filt,4);
                bF = dat(filt,6);
                minT = min([bS;temp(:,2)]);
                maxT = max([bE;temp(:,1)]);
                
                % detect animals that weren't tracked for a full 15seconds
                if size(temp,1) <= length(fnames)
                    tS = [nan];
                    tE = [nan];
                    % old code - failed when "temp" had fewer than 4 rows
%                     bS = nan(length(fnames),1);
%                     bE = nan(length(fnames),1);
%                     bT = unique(dat(:,end));
                    % assign nan to the 
                    bS = nan(size(temp,1),1);
                    bE = nan(size(temp,1),1);
                    bA = nan(size(temp,1),1);
                    bF = nan(size(temp,1),1);
                    bT = temp(:,end);
                else
                    if minT > (temp(1,2) - temp(1,3))
                        tS = temp(1,2) - temp(1,3);
                    else
                        tS = minT;
                    end

                    if maxT < (temp(end,1) + temp(end,3))
                        tE = temp(end,1) + temp(end,3);
                    else
                        tE = maxT;
                    end
                    bT = temp(:,end);
                end

                tStart = [tStart tS];
                tEnd = [tEnd tE];
                bStart = [bStart; bS];
                bEnd = [bEnd; bE];
                bAmp = [bAmp; bA];
                bFreq = [bFreq; bF];
                bType = [bType; bT];
            end
            % filt out duplicate behaviours
            if all(isnan(aniID))
                fprintf("no usable salam data \n")
                return
            end
            [~,f,~] = unique([aniID,bStart,bEnd,bType],'rows','stable');
            bStart = bStart(f);
            bEnd = bEnd(f);
            bAmp = bAmp(f);
            bFreq = bFreq(f);
            bType = bType(f);
            
            % assign variables to object
            obj.aniID = uniID';
            obj.track_start = tStart;
            obj.track_end = tEnd;
            temp_arr = [];
            for ii = 1:numel(fnames)
                filt = bType == ii;
                temp = struct();
                temp.behaviour = fnames{ii};
                temp.id = aniID(filt)';
                temp.start = bStart(filt)';
                temp.end = bEnd(filt)';
                temp.amplitude = bAmp(filt)';
                if strcmp(fnames{ii},"peran")
                    temp.frequency = bFreq(filt)';
                else
                    temp.frequency = [];
                end
                % count the number of events per animal (occaisionally 0 - have to use "ismember")
                [B,I] = ismember(temp.id,obj.aniID);
                [C,~,ic] = unique(I');
                counts = accumarray(C(ic),1,[numel(obj.aniID),1])';
                % count the number of events per animal (fails for values less than 1)
%                 [~,~,ic] = unique(temp.id');
%                 counts = accumarray(ic,1)';
                
                temp.track_start = repelem(tStart,1,counts);
                temp.track_end = repelem(tEnd,1,counts);
                temp_arr = vertcat(temp_arr,temp);
            end
            
            obj.behaviour = temp_arr;
            obj.behaviourtype = string([fnames{:}]);
            
            obj.raw_data = [];
            fprintf('\t----%s - salam compiled \n',obj.get_full_experiment)
        end
        
        %% compile jaaba
        function obj = compile_jaaba(obj)
            if isempty(obj.raw_data)
                fprintf('\t\t ...no jaaba data currently loaded\n')
                return
            end
            temp = obj.raw_data;
            timestamps = temp.timestamps;

            clear jaaba
            obj.aniID = [temp.trx.id];
%             [~,~,uniID] = unique(obj.aniID,'stable');
%             obj.uniID = uniID';
            if ~any(contains(fieldnames(temp),'allScores'))
                fprintf('\t\t ...no jaaba "allScores" file\n')
                return
            end
            obj.track_start = timestamps(temp.allScores.tStart);
            obj.track_end = timestamps(temp.allScores.tEnd);

            empties = cellfun('isempty',temp.allScores.t0s);
            bStart = cellfun(@(x) timestamps(x), temp.allScores.t0s, 'UniformOutput',false);
            bEnd = cellfun(@(x) timestamps(x), temp.allScores.t1s, 'UniformOutput',false);
            bStart(empties) = {nan};
            bEnd(empties) = {nan};
            
            obj.behaviour.behaviour = temp.behaviorName;
            counts = cellfun(@numel, bEnd);
            obj.behaviour.id = repelem(obj.aniID,counts);
            
            obj.behaviour.start = [bStart{:}];
            obj.behaviour.end = [bEnd{:}];

            obj.behaviour.amplitude = [];
            obj.behaviour.frequency = [];
            
            obj.behaviour.track_start = repelem(obj.track_start,counts);
            obj.behaviour.track_end = repelem(obj.track_end,counts);
            
            obj.timeseries.et = timestamps;
            obj.behaviourtype = string(temp.behaviorName);

            obj.raw_data = [];
            fprintf('\t----%s - jaaba compiled \n',obj.get_full_experiment)
        end
        
        %% Compile JB
        function obj = compile_jb(obj)
            if isempty(obj.raw_data)
                fprintf('\t\t ...no jb data currently loaded\n')
                temp = [];
                return
            end

            trx = obj.raw_data.trx;

            obj.aniID = [trx.numero_larva_num];
            obj.track_start = arrayfun(@(x) min(x.t), trx)';
            obj.track_end = arrayfun(@(x) max(x.t), trx)';

            beh_name = {'run','cast','stop','hunch','back','roll'};
            beh_type = {'t_start_stop',...
                't_start_stop_large_small'};

            beh_type = beh_type(ismember(beh_type,fieldnames(trx)));
            
            if find(contains(beh_type,'t_start_stop_large_small'))
                cells = {trx.t_start_stop_large_small};
                
                t_start_stop_weak = cellfun(@(x) x(1:2:12), cells, 'UniformOutput', false);
                t_start_stop_strong = cellfun(@(x) x(2:2:12), cells, 'UniformOutput', false);
                
                [trx.t_start_stop_weak] = t_start_stop_weak{:};
                [trx.t_start_stop_strong] = t_start_stop_strong{:};
                
                beh_type = beh_type(~strcmp(beh_type,'t_start_stop_large_small'));
                beh_type = [beh_type,'t_start_stop_weak','t_start_stop_strong'];
            end
            
            tokens = regexp(beh_type,'(^\w)|[_](\w)','tokens');
            beh_code = cellfun(@(x) strcat(x{:}), tokens);
            
            %%%%% edit 26/08/20
            extra_behs = {'head_wiggle','wiggling'};
            for ii = 1:numel(extra_behs)
                if ismember(extra_behs{ii},fieldnames(trx))
                    vect_beh = {trx.(extra_behs{ii})};
                    vect_beh = cellfun(@(x) x>0, vect_beh, 'UniformOutput', false);
                    
                    starts = cellfun(@(x) find(diff([0,x']) > 0), vect_beh, 'UniformOutput', false);
                    stops = cellfun(@(x) find(diff([x',0]) < 0), vect_beh, 'UniformOutput', false);
                    
%                     vect_beh{jj}'
%                     for jj = 1:numel(starts)
%                         [starts{jj};stops{jj}]
%                     end

                    beh_cell = cellfun(@(t,x,y) t([x;y])', {trx.t}, starts, stops, 'UniformOutput', false);
                    
                    ff = arrayfun(@(x,y) [x{:},y], {trx.t_start_stop}, beh_cell, 'UniformOutput', false);
                    [trx.t_start_stop] = ff{:};
                end
                beh_name = [beh_name, extra_behs{ii}];
            end
            %%%%%
            
            obj_arr = [];
            for ii = 1:length(beh_type)
                len = length(trx(1).(beh_type{ii}));
                beh_cell = [trx.(beh_type{ii})];
                idx = cellfun(@isempty,beh_cell);
                beh_cell(idx) = {[nan,nan]};
                beh_cell = cellfun(@transpose, beh_cell,'UniformOutput', false);
                beh_starts = cellfun(@(x) x(1,:),beh_cell,'UniformOutput',false, 'ErrorHandler',@(S,varargin) []);
                beh_starts = reshape(beh_starts,len,[]);
                beh_ends = cellfun(@(x) x(2,:),beh_cell,'UniformOutput',false, 'ErrorHandler',@(S,varargin) []);
                beh_ends = reshape(beh_ends,len,[]);
                for jj = 1:len
                    temp = struct();
                    temp.behaviour = [beh_name{jj},'_',beh_code{ii}];
                    counts = cellfun(@numel, beh_starts(jj,:));
                    temp.id = repelem(obj.aniID,1,counts);
                    temp.start = [beh_starts{jj,:}];
                    temp.end = [beh_ends{jj,:}];
                    temp.amplitude = [];
                    temp.frequency = [];
                    temp.track_start = repelem(obj.track_start,1,counts);
                    temp.track_end = repelem(obj.track_end,1,counts);
                    obj_arr = vertcat(obj_arr,temp);
                end
            end
            obj.raw_data = [];  
            obj.behaviour = obj_arr;
            fprintf('\t----%s - jb compiled \n',obj.get_full_experiment)
        end
        
        %% Compile choreography
        function obj = compile_choreography(obj)
            temp = obj.raw_data;
            [id,~,idx] = unique(temp.id);
            counts = accumarray(idx,1);
            fnames = fieldnames(temp);
            fnames = fnames(fnames~="id");

            chore = struct('id',num2cell(id));
            for ii = 1:numel(fnames)
                cells = mat2cell(temp.(fnames{ii})',1,counts);
                [chore.(fnames{ii})] = cells{:};
            end

            obj.aniID = id';
%             obj.uniID = [1:numel(id)];

            [tmin,tmax] = cellfun(@bounds, {chore.et});
            obj.track_start = tmin;
            obj.track_end = tmax;

            obj.behaviourtype = "choreography";
            obj.timeseries = chore;
            obj.raw_data = [];
            fprintf('\t----%s - choreography compiled \n',obj.get_full_experiment)
        end
        
        %% Compile MWT
        function obj = compile_mwt(obj)
            temp = obj.raw_data;
            [id,~,idx] = unique(temp.aniID);
            counts = accumarray(idx,1);
            fnames = {"time","x_cont","y_cont","npts","outline"};

            mwt = struct('id',num2cell(id)');
            for ii = 1:numel(fnames)
                cells = mat2cell(temp.(fnames{ii}),1,counts);
                [mwt.(fnames{ii})] = cells{:};
            end

            obj.aniID = id;
%             obj.uniID = [1:numel(id)];

            [tmin,tmax] = cellfun(@bounds, {mwt.time});
            obj.track_start = tmin;
            obj.track_end = tmax;

            obj.behaviourtype = "blobs";
            obj.timeseries = mwt;
            [obj.raw_data] = [];
            fprintf('\t----%s - mwt compiled \n',obj.get_full_experiment)
        end
        
        %%
        function obj = update_unique_ids(obj)
            % work out the unique ids across timestamps
            counts = arrayfun(@(x) numel(x.aniID),obj);
            ids = [obj.aniID];
            exps = [obj.group_experiment];
            exps = repelem(exps,counts);
            [~,~,timestamp_index] = unique(exps,'stable');
            
            % deal the unique ids back to the class instances
            [lkup,~,unids] = unique([exps;ids]','rows','stable');
            unids_cell = mat2cell(unids',1,counts);
            [obj.uniID] = unids_cell{:};
            index_cell = mat2cell(timestamp_index',1,counts);
            [obj.timestamp_index] = index_cell{:};
            
%             uniID_timestamp = cellfun(@(x) unique(unids(exps == x)),...
%                 num2cell(unique(exps,'stable')),...
%                 'UniformOutput', false);

            % add the unique ids at the behaviour/timeseries level
            for ii = 1:numel(obj)
                if any(contains(fieldnames(obj(ii).timeseries),'id'))
                    counts = arrayfun(@(x) numel(x.id), obj(ii).timeseries);
                    [~,idx] = ismember([obj(ii).timeseries.id],[obj(ii).aniID]);
                    unids = mat2cell(obj(ii).uniID(idx),1,counts);
                    [obj(ii).timeseries.uniID] = unids{:};
                end

                if any(contains(fieldnames(obj(ii).behaviour),'id'))
                    counts = arrayfun(@(x) numel(x.id), obj(ii).behaviour);
                    [~,idx] = ismember([obj(ii).behaviour.id],[obj(ii).aniID]);
                    unids = mat2cell(obj(ii).uniID(idx),1,counts);
                    [obj(ii).behaviour.uniID] = unids{:};
                end
            end
        end
        
        
        %% extract full genotype
        function full_genotype = get_full_genotype(obj)
            full_genotype = strcat([obj.driver],'@',[obj.effector]);
        end

        %% -------------------------------------------------------
        function full_protocol = get_full_protocol(obj)
            full_protocol = strcat([obj.protocol1],...
                '#',[obj.protocol2],...
                '#',[obj.protocol3],...
                '#',[obj.protocol4]);
        end

        %% -------------------------------------------------------
        function full_timestamp = get_full_timestamp(obj)
            full_timestamp = strcat([obj.date],'_',[obj.time]);
        end  
        
        %% -------------------------------------------------------
        function full_experiment = get_full_experiment(obj)
            full_experiment = strcat([obj.get_full_timestamp],'@',...
                [obj.get_full_genotype],'@',...
                [obj.get_full_protocol]);
        end  
        %% Check filesize
        function tf = check_file(obj)
            file = dir(obj.filepath(1));         
            tf = file.bytes ~= 0;
        end
        
        %% generate figure directory
        function obj = update_figure_directories(obj,params)
            [un_prot,~,idx_prot] = unique([obj.group_effector_protocol]);
            
            for ii = 1:length(un_prot)
                idx = idx_prot == ii;

                dates = double([obj(idx).date]);
                n_timestamps = numel(unique(obj(idx).get_full_timestamp));
                effector = [obj(idx).effector];
                protocol = obj(idx).get_full_protocol;

                fig_dir = sprintf("%d_%d@%s@%s@x%d",...
                    min(dates),max(dates),effector(1),...
                    protocol(1),n_timestamps);
                
                fig_dir = fullfile(params.directories.master,'figures',fig_dir);
                
                dirs = num2cell(repelem(fig_dir,numel(idx),1));
                [obj(idx).figure_directory] = dirs{:};
            end
        end
        
        %% merge_dataClasses
        one_dataClass = merge_dataClasses(obj);
                
        dc_genotype = load_bygenotype(obj,genotype_number)
        
        dc_genotype = load_bygroup(obj)
        
        dc_group = load_groups(obj,filter_date_tf)
        
        protocol_struct = parse_protocol(obj)
        
        plot_struct = append_2_plotstruct(obj,plot_struct,params,behaviour,metric,instance,method,frame)
        
        plotStructUpdated = append_choreography_metric(obj,plotStruct,feature,metric,frame,instance)
        
        data_array = get_behaviour_data(obj,behaviour,params)

        obj = filter_by_size(obj,area_threshold)
        
        [ax,figure_path] = plot_ridgeline_timeseries(obj,params,feature,frame_start,frame_end,y_scaler)
        
        [ax,figure_path] = plot_timeseries(obj,params,feature,frame_start,frame_end,y_limits,calculation,colormap)
        
        [ax,fig_path] = plot_timeseries_byday(obj,params,feature,frame_start,frame_end,y_limits,calculation,colormap)
        
        [ax,figure_path] = plot_ethogram(obj,params,behaviours_cell,frame_start,frame_end,offset)
        
        [ax,fig_path] = plot_behaviour_timeseries(obj,params,behaviour,bin_width,x_lims,y_lims,cmap)
                
        plot_stimulation_time(obj,params)
        
        %% update groups
        function obj = update_groups(obj)
            [~,~,grp_prot] = unique(obj.get_full_protocol);
            [~,~,grp_eff] = unique([obj.effector]);
            [~,~,grp_dri] = unique([obj.driver]);
            [~,~,grp_efp] = unique(strcat([obj.effector],'@',obj.get_full_protocol));
            [~,~,grp_gen] = unique(obj.get_full_genotype);
            [~,~,grp_exp] = unique(obj.get_full_experiment');

            grp_prot = num2cell(grp_prot);
            grp_eff = num2cell(grp_eff);
            grp_dri = num2cell(grp_dri);
            grp_efp = num2cell(grp_efp);
            grp_gen = num2cell(grp_gen);
            grp_exp = num2cell(grp_exp);

            [obj.group_genotype] = grp_gen{:};
            [obj.group_driver] = grp_dri{:};
            [obj.group_effector] = grp_eff{:};
            [obj.group_effector_protocol] = grp_efp{:};
            [obj.group_protocol] = grp_prot{:};
            [obj.group_experiment] = grp_exp{:};
        end
        
        %%
        function [n_groups,group_numbers] = get_group_numbers(obj,varargin)
            args = varargin;
            
            if any(contains(args,'protocol'))
                args = [args(~contains(args,'protocol')),...
                    'protocol1','protocol2','protocol3','protocol4'];
            end
            
            if any(contains(args,'timestamp'))
                args = [args(~contains(args,'timestamp')),...
                    'date','time'];
            end

            strngs = cellfun(@(x) [obj.(x)], args, 'UniformOutput', false);
            strngs = [strngs;repelem({'@'},1,numel(args))];

            [C,~,group_numbers] = unique(strcat(strngs{:})','stable');
            n_groups = numel(C);
        end
    
    end
        
    %% METHODS - STATIC
    
    methods (Static)
        
        %% Get Files
        function files = get_files(params)
            files = [];
            fields = struct2cell(params.filetypes);
            fields = [fields{:}];
            n = 1;
            while length(fields) >= n
                f = dir(fullfile(params.directories.data,'**',['*',fields{n}]));
                files = [files; f];
                n = n+1;
            end
            [C,ia,ic] = unique(fullfile({files.folder},{files.name}));
            files = files(ia);
        end
        
        %% get_filetype
        function [types,idx] = get_filetype(filepaths)
            types = [];
            idx = [];
            for ii = 1:numel(filepaths)
                [fold,name,ext] = fileparts(string(filepaths(ii)));
                if startsWith(ext,".blo")
                    types = [types, "mwt"];
                    idx = [idx, ii];
                elseif startsWith(ext,".txt")
                    if contains(lower(name),"chore")
                        types = [types, "choreography"];
                        idx = [idx, ii];
                    elseif contains(lower(name),"animal_stats")
                        types = [types, "salam"];
                        idx = [idx, ii];
                    end
                elseif startsWith(ext,".mat")
                    if contains(name,"trx")
                        file_check = dir(fold);
                        if any({file_check.name} == "perframe")
                            types = [types, "jaaba"];
                            idx = [idx, ii];
                        else
                            types = [types, "jb"];
                            idx = [idx, ii];
                        end
                    end
                end
            end
        end
        
        %% Parse Filepaths
        function [file_details,error_files] = parse_filepaths(filepath_full)
            
            expr = ['(?<date>\d\d\d\d\d\d\d\d)[_]'...
                '(?<time>\d\d\d\d\d\d)[\/\\\@]'...
                '(?<driver>\w+)[@]'...
                '(?<effector>\w+)[@]'...
                '(?<rig>\w+)[@]'...
                '(?<protocol1>[\w\d\_]+)'...
                '[#](?<protocol2>[\w\d\_]+)'...
                '[#](?<protocol3>[\w\d\_]+)'...
                '[#](?<protocol4>[\w\d\_]+)'...
                '|'...
                '[\@\\\/](?<driver>\w+)'...
                '[\@](?<effector>\w+)'...
                '[\@\\\/](?<rig>t\d+)'...
                '[\@\\\/](?<protocol1>[\w\d\_]+)'...
                '[#](?<protocol2>[\w\d\_]+)'...
                '[#](?<protocol3>[\w\d\_]+)'...
                '[#](?<protocol4>[\w\d\_]+)',...
                '[\@\\\/](?<date>\d\d\d\d\d\d\d\d)'...
                '[\_](?<time>\d\d\d\d\d\d)'...
                '|'...
                '[\@\\\/](?<rig>t\d+)'...
                '[\@\\\/](?<driver>\w+)'...
                '[\@](?<effector>\w+)'...
                '[\@\\\/](?<protocol1>[\w\d\_]+)'...
                '[#](?<protocol2>[\w\d\_]+)'...
                '[#](?<protocol3>[\w\d\_]+)'...
                '[#](?<protocol4>[\w\d\_]+)'...
                '[\@]\d+'...
                '[\@\\\/](?<date>\d\d\d\d\d\d\d\d)'...
                '[\_](?<time>\d\d\d\d\d\d)'...
                '|'...
                '(?<date>\d{8})[_]'...
                '(?<time>\d{6})[\/\\\@]'...
                '(?<driver>\w+)[@]'...
                '(?<effector>\w+)[@]'...
                '(?<protocol1>[\w\d\_]+)'...
                '[#](?<protocol2>[\w\d\_]+)'...
                '[#](?<protocol3>[\w\d\_]+)'...
                '[#](?<protocol4>[\w\d\_]+)']; 
            
            file_details = regexp(string(filepath_full),expr,'names');
            
            if length(file_details) > 1
                error_files = cellfun(@isempty,file_details);
                file_details = vertcat(file_details{:});
            else
                error_files = [];
            end
        end
        
        %% Save Data Batch
        function save_data_batch(obj_array,params) 
            savename = fullfile(params.directories.data_compiled,strcat(obj_array(1).get_full_experiment, '.mat'));
            for ii = 1:numel(obj_array)
                content.(obj_array(ii).pipeline).(obj_array(ii).behaviour) = true;
            end

            if isfile(savename)
                old_content = load(savename,"content");
                tf = isequaln(old_content.content,content);
                if ~tf
                    save(savename,"obj_array","content");
                    fprintf('\t----%s - data updated and saved \n',obj_array(1).get_full_experiment)
                else
                    fprintf('\t----%s - data already saved \n',obj_array(1).get_full_experiment)
                end
            else
                save(savename,"obj_array","content");
                fprintf('\t----%s - data saved \n',obj_array(1).get_full_experiment)
            end
        end
        
        %% generate_contents_file
        generate_contents_file(obj,contentfile,pipelines,params);
        
        %% generate_contents_file
        update_contents_file(obj,contentfile,files_to_append,pipelines,params);
        
        %% filterby_contents_file
        filelist_filtered = filterby_contents_file(filelist,contentfile,pipelines,params);
        
        [ax,figure_path] = plot_grouped_beeswarm(obj,plot_struct,sorted,normalise_control_tf,control_genotype,highlight_genotype,y_label)
        
        [ax,figure_path,stats] = plot_violins(obj,plot_struct,sortOrder,normalise_control_tf,control_driver,highlight_driver,y_label)
        
        save_figure_catch(ax,figure_path);
        
    end
end