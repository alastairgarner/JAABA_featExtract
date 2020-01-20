%% dataContainer

% Other m-files required: 
% Subfunctions: 
% MAT-files required: 
% See also:

% Author: Alastair Garner
% email: alastairgarner@outlook.com
% Website: https://github.com/alastairgarner/
% Sep 2019; Last revision: 

classdef dataClass_edit
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
%         uniID = []
        
        tStart = []
        tEnd = []
        
        behaviour = struct();
        timeseries = struct();
        object = struct();
        
        raw_data = [];
    end
    
    %% METHODS NORMAL
    methods
        
        %% Construct object from file paths
        function obj = dataClass(filepaths)
            if nargin ~= 0
                [split_filepaths,names,exts] = cellfun(@(x) fileparts(x),filepaths,'UniformOutput', false);        
                [unique_filepaths,~,indicies] = unique(split_filepaths,'stable');
                issalam = contains(names, "animal_stats");
%                 [~,~,ic] = unique(split_filepaths(isblob),'stable');
                indicies(issalam) = nan;
                [unique_filepaths,~,indicies] = unique(indicies);
                
                elems = numel(unique_filepaths);
                obj(elems,1) = dataClass();
                % 
                for ii = 1:numel(unique_filepaths)
                    idx = indicies == ii;
                    obj(ii).filepath = string(filepaths(idx))';
                    temp = dataContainer.parse_filepaths(filepaths(idx));
                    % merge struct to object
                    for fn = fieldnames(temp)'    %enumerat fields
                        try
                          obj(ii).(fn{1}) = temp.(fn{1});   %and copy
                        catch
                          warning('Could not copy field %s', fn{1});
                        end
                    end
                end
            end
            
        end
        
        %% Save Data
        function save_data(obj)
            obj.create_save_variable()
            savename = fullfile("data_compiled",strcat(obj.get_full_experiment, '.mat'));
            if isfile(savename)
                contents = load(savename,"content");
                content = contents.content;
                content.(obj.pipeline).(obj.behaviour) = true;
                save(savename,obj.behaviour,"content","-append");
            else
                content.(obj.pipeline).(obj.behaviour) = true;
                save(savename,obj.behaviour,"content");
            end
        end
        
        
        %% Create variable for given datacontainer
        function create_save_variable(obj)
            temp = struct(obj);
            assignin("caller",obj.behaviour,temp);
        end
        
        %% Load data
        function obj = load_data(obj)
            for ii = 1:numel(obj.filepath)
                [fold,name,ext] = fileparts(obj.filepath(ii));
                if startsWith(ext,".blo")
                    obj.pipeline = "mwt";
                    obj = obj.load_blobsdata();
                    return
                elseif startsWith(ext,".txt")

                    if contains(lower(name),"chore")
                        obj.pipeline = "choreography";
                        obj = obj.load_choredata();
                        return
                    elseif contains(lower(name),"animal_stats")
                        obj.pipeline = "salam";
                        obj = obj.load_salamdata();
                        return
                    end

                elseif startsWith(ext,".mat")

                    if contains(name,"trx")
                        file_check = dir(fold);
                        if any({file_check.name} == "perframe")
                            obj.pipeline = "jaaba";
                            obj = obj.load_jaabadata();
                            return
                        else
                            obj.pipeline = "jb";
                            obj = obj.load_jbdata();
                            return
                        end
                    end

                end
            end
        end
        
        %% Load blobs data
        function obj = load_blobsdata(obj)
            filepaths = obj.filepath;
            [fpath,~,~] = fileparts(filepaths(1));
            % check if summary file exists - if not then cancel
            d = dir(fullfile(fpath,'*.summary'));
            if ~size(d,1)
                return
            end
            % set formatspec of blobs files and column names
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
                % open and read blobs file
                fileID = fopen(filepaths(ii));
                datA = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'CollectOutput', 1, 'ReturnOnError', false, 'treatas', '%');
                fclose(fileID);
                % get columns containing double-type values
                datmat = [datA{[1,3]}];
                % get animal IDs from rogue rows
                filt = isnan(datmat(:,1));
                aniID = datmat(filt,2);
                % filter out rogue rows
                datmat = datmat(~filt,:);
                % get outlines in string-type
                outlines = datA{4}(~filt,1);
                % make new column with animal IDs
                diffs = diff(find([filt' 1]))-1;
                aniIDs = repelem(aniID,diffs,1);
                datmat = [aniIDs datmat];
                % append new data
                datfull = vertcat(datfull,datmat);
                outlinesfull = [outlinesfull; outlines];
                clear datmat
            end
            % sort data by animal ID
            [datfull,idx] = sortrows(datfull,[1,2]);
            outlinesfull = outlinesfull(idx);
            % reorganise into structure array by animal ID
            [~,~,ic] = unique(datfull(:,1));
            counts = accumarray(ic,1)';
            datastruct = repelem(struct(),numel(counts),1);
            for ii = 1:length(columns)
                datacells = mat2cell(datfull(:,ii)',1,counts);
                [datastruct.(columns{ii})] = datacells{:};
            end
            % add back the outline data
            datastrings = mat2cell(string(outlinesfull)',1,counts);
            [datastruct.outline] = datastrings{:};
            % assign to object
            obj.timeseries = datastruct;
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
            
            [~,~,ic] = unique(dataArray{1});
            counts = accumarray(ic,1)';
            datastruct = repelem(struct(),numel(counts),1);
            columns = vertcat({'aniID';'et'},chores);
            for ii = 1:length(columns)
                datacells = mat2cell(dataArray{ii}',1,counts);
                [datastruct.(columns{ii})] = datacells{:};
            end
            obj.timeseries = datastruct;

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
                dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow, 'ReturnOnError', false,'treatAsEmpty', {'NA'});
                fclose(fileID);

                % Restructure crawl data to match other data
                filler = zeros(length(dataArray{1}),1);
                if ~all(isnan(dataArray{end-1}))
                    dataArray = dataArray([1:5,9,7,10,11,6,8]);
                end

                salam_matrix = [dataArray{[2:11]}];
                obj.raw_data.(behaviour(ii)) = salam_matrix;
            end
            obj.pipeline = "salam";
            obj.behaviour = "";
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
            tStart = []; tEnd = []; bStart = []; bEnd = []; bType = [];

            for jj = 1:size(uniID,1)
                filt = aniID == uniID(jj);
                temp = dat(filt,[2,3,8,11]);
                temp = sortrows(temp,1);
                bS = dat(filt,7);
                bE = dat(filt,7) + dat(filt,5);
                minT = min([bS;temp(:,2)]);
                maxT = max([bE;temp(:,1)]);

                if size(temp,1) <= length(fnames)
                    tS = [nan];
                    tE = [nan];
                    bS = nan(length(fnames),1);
                    bE = nan(length(fnames),1);
                    bT = unique(dat(:,end));
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
                bType = [bType; bT];
            end

            obj.aniID = uniID';
            [~,~,uniID] = unique(obj.aniID,'stable');
            obj.uniID = uniID';
            obj.tStart = tStart;
            obj.tEnd = tEnd;
            obj.bID = aniID';
            obj.bStart = bStart';
            obj.bEnd= bEnd';
            obj.bAmp = [];
            obj.behaviour = string([fnames{:}]);
            
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
            [~,~,uniID] = unique(obj.aniID,'stable');
            obj.uniID = uniID';
            obj.tStart = timestamps(temp.allScores.tStart);
            obj.tEnd = timestamps(temp.allScores.tEnd);

            empties = cellfun('isempty',temp.allScores.t0s);
            bStart = cellfun(@(x) timestamps(x), temp.allScores.t0s, 'UniformOutput',false);
            bEnd = cellfun(@(x) timestamps(x), temp.allScores.t1s, 'UniformOutput',false);
            bStart(empties) = {nan};
            bEnd(empties) = {nan};
            obj.bStart = [bStart{:}];
            obj.bEnd = [bEnd{:}];

            counts = cellfun(@numel, bEnd);
            obj.bID = repelem(obj.aniID,counts);

            obj.timeseries.elapstime = timestamps;
            obj.behaviour = string(temp.behaviorName);

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
            [~,~,uniID] = unique(obj.aniID,'stable');
            obj.uniID = uniID';
            obj.tStart = arrayfun(@(x) min(x.t), trx)';
            obj.tEnd = arrayfun(@(x) max(x.t), trx)';

            beh_name = {'run','cast','stop','hunch','back','roll','beh7','beh8','beh9','beh10','beh11','beh12'};
            beh_type = {'t_start_stop',...
                't_start_stop_large',...
                't_start_stop_large_small'};
            beh_type = {'t_start_stop'};

            tokens = regexp(beh_type,'(^\w)|[_](\w)','tokens');
            tokens = cellfun(@(x) string([x{:}]), tokens, 'UniformOutput',false);
            beh_code = cellfun(@(x) strcat(x{:}),tokens, 'UniformOutput', false);


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
                    temp = obj;
                    temp.raw_data = [];
                    counts = cellfun(@numel, beh_starts(jj,:));
                    temp.bID = repelem(obj.aniID,1,counts);
                    temp.behaviour = [beh_name{jj},'_',beh_code{ii}];
                    temp.bStart = [beh_starts{jj,:}];
                    temp.bEnd = [beh_ends{jj,:}];
                    obj_arr = vertcat(obj_arr,temp);
                end
            end
            obj.raw_data = [];  
            obj = obj_arr;
            fprintf('\t----%s - jb compiled \n',temp.get_full_experiment)
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
            obj.uniID = [1:numel(id)];

            [tmin,tmax] = cellfun(@bounds, {chore.et});
            obj.tStart = tmin;
            obj.tEnd = tmax;

            obj.behaviour = "choreography";
            obj.timeseries = chore;
            obj.raw_data = [];
            fprintf('\t----%s - choreography compiled \n',obj.get_full_experiment)
        end
        
        %%
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
            obj.uniID = [1:numel(id)];

            [tmin,tmax] = cellfun(@bounds, {mwt.time});
            obj.tStart = tmin;
            obj.tEnd = tmax;

            obj.behaviour = "blobs";
            obj.timeseries = mwt;
            [obj.raw_data] = [];
            fprintf('\t----%s - mwt compiled \n',obj.get_full_experiment)
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
        
    end
end

