%% dataStruct

% Other m-files required: 
% Subfunctions: 
% MAT-files required: 
% See also:

% Author: Alastair Garner
% email: alastairgarner@outlook.com
% Website: https://github.com/alastairgarner/
% Sep 2019; Last revision: 

classdef dataStruct
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        date = "";
        time = "";
        driver = "";
        effector = "";
        protocol1 = "";
        protocol2 = "";
        protocol3 = "";
        protocol4 = "";
        rig = "";
        group = [];
        iscontrol = false;
                
        has
        
        data_mwt
        data_chore
        data_salam
        data_jaaba
        data_jb
        
        raw_blobs
        raw_chore
        raw_salam
        raw_jaaba
        raw_jb
        
        files_blobs
        files_chore
        files_salam
        files_jaaba
        files_jb
    end
    
    %% METHODS NORMAL
    methods
        
        %% Construct object from file paths
        function obj = dataStruct(file_paths)
            %  CONSTRUCTOR
            %    Detailed explanation goes here
            
            if nargin == 0
                obj.date = "";
                return
            end
            file_structure = dataStruct.parse_filepaths(file_paths(1));
            
            obj.date = file_structure.date;
            obj.time = file_structure.time;
            obj.driver = file_structure.driver;
            obj.effector = file_structure.effector;
            obj.protocol1 = file_structure.protocol1;
            obj.protocol2 = file_structure.protocol2;
            obj.protocol3 = file_structure.protocol3;
            obj.protocol4 = file_structure.protocol4;
            obj.rig = file_structure.rig;
            
            
            [fpath,fname,fext] = cellfun(@fileparts, file_paths, 'uni', false);
    
            filt_blob = startsWith(fext,'.blo');
            filt_Chore_dat = startsWith(fext,'.dat');
            filt_Chore_txt = endsWith(fname,'compiledChore') & startsWith(fext,'.tx');
            filt_salam = ~endsWith(fname,'compiledChore') & startsWith(fext,'.tx');
            filt_jaaba = contains(fpath,'jaaba') & startsWith(fext,'.mat');
            filt_jb = contains(fpath,'jb') & startsWith(fext,'.mat');
            
            obj.files_blobs = string(file_paths(filt_blob)');
            obj.files_chore = string(file_paths(filt_Chore_txt)');
            obj.files_salam = string(file_paths(filt_salam)');
            obj.files_jaaba = string(file_paths(filt_jaaba)');
            obj.files_jb = string(file_paths(filt_jb)');
            
            obj.has.blobs = ~isempty(obj.files_blobs);
            obj.has.chore = ~isempty(obj.files_chore);
            obj.has.salam = ~isempty(obj.files_salam);
            obj.has.jaaba = ~isempty(obj.files_jaaba);
            obj.has.jb = ~isempty(obj.files_jb);
            
        end
        
        %% Load data, Compile data, Save data
        function load_compile_save(obj,params)
            
            for ii = 1:length(obj)
                temp = obj(ii);
                temp = temp.load_blobsdata();
                temp = temp.load_choredata();    
                temp = temp.load_salamdata();
                temp = temp.load_jaabadata();
                temp = temp.load_jbdata();

                temp = temp.compile_all();
                temp = temp.clear_raw();
                
                fields = fieldnames(temp);
                f = startsWith(fields,'data_');
                datafields = sort(fields(f));
                content = cellfun(@(x) ~isempty(temp.(x)),datafields,'UniformOutput',false);
                content = [datafields,content]';
                content = struct(content{:});
                
                savename = fullfile(params.directories.dataprocessed,strcat(temp.get_full_experiment, '.mat'));
                temp = struct(temp);
                save(savename,'temp','content');
                clear temp
            end
        end

        %% Load blobs data
        function obj = load_blobsdata(obj)
            [fpath,~,~] = fileparts(obj.files_blobs(1));
            d = dir(fullfile(fpath,'*.summary'));
            if isempty(obj.files_blobs) | ~size(d,1)
                return
            end
            
            blobfiles = obj.files_blobs;

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
            for ii = 1:length(blobfiles)
                fileID = fopen(blobfiles(ii));
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
            obj.raw_blobs = cell2struct(datcell,columns);
            obj.raw_blobs.outline = outlinesfull';
            
            fprintf('\t----%s - blobs loaded \n',obj.get_full_experiment)
        end
        
        %% Load chore data
        function obj = load_choredata(obj)
            if isempty(obj.files_chore)
                fprintf('\t\t... no chore data\n')
                return
            end
            
            filepath = obj.files_chore;
            fileID = fopen(filepath,'r');
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

            obj.raw_chore.id = dataArray{1};
            obj.raw_chore.et = dataArray{2};
            for jj = 1:length(dataArray)-2
                field_name = chores{jj};
                obj.raw_chore.(field_name) = dataArray{jj+2};
            end
            fprintf('\t----%s - chore loaded \n',obj.get_full_experiment)
        end
        
        %% Load feature extraction (salam) data
        function obj = load_salamdata(obj)
            if isempty(obj.files_salam)
                fprintf('\t\t... no salam data\n')
                return
            end
            
            delimiter = '\t';
            startRow = 1;
            formatSpec = '%s%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
            
            filenames = obj.files_salam;
            
            expr = '[\_]([A-Za-z0-9]+).txt';
            behaviour = regexp(filenames,expr,'tokens','once');
            behaviour = [behaviour{:}];
            
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
                obj.raw_salam.(behaviour(ii)) = salam_matrix;
            end
            fprintf('\t----%s - salam loaded \n',obj.get_full_experiment)
        end
        
        %% load jaaba data
        function obj = load_jaabadata(obj)
            if isempty(obj.files_jaaba)
                fprintf('\t\t... no jaaba data\n')
                return
            end
            
            for ii = 1:length(obj.files_jaaba)
                data = load(obj.files_jaaba(ii));
                f = fieldnames(data);
                for jj = 1:length(f)
                   obj.raw_jaaba.(f{jj}) = data.(f{jj});
                end
            end
            fprintf('\t----%s - jaaba loaded \n',obj.get_full_experiment)
        end
        
        %% load jb data
        function obj = load_jbdata(obj)
            if isempty(obj.files_jb)
                fprintf('\t\t... no jb data\n')
                return
            end
            for ii = 1:length(obj.files_jb)
                data = load(obj.files_jb(ii));
                f = fieldnames(data);
                for jj = 1:length(f)
                   obj.raw_jb.(f{jj}) = data.(f{jj});
                end
            end
            fprintf('\t----%s - jb loaded \n',obj.get_full_experiment)
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
        
        %% Compile Salam Data
        function obj = compile_salam(obj)
            if isempty(obj.raw_data)
                fprintf('\t\t ...no salam data currently loaded\n')
                data = [];
                return
            end
            
            fnames = fieldnames(obj.raw_data);
            [rows,cols] = structfun(@size, obj.raw_data);
            behID = repelem([1:length(fnames)]',rows);
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
            
            data.pipeline = 'salam';
            data.aniID = uniID';
            [~,~,uniID] = unique(data.aniID,'stable');
            data.uniID = uniID';
            data.tStart = tStart;
            data.tEnd = tEnd;
            
            idcounts = nonzeros(accumarray(aniID,1));
            idcounts(idcounts<length(fnames)) = length(fnames);
            aniID = repelem(uniID,idcounts);
            
            bMerged = [aniID, bStart, bEnd, bType];
            [C,ia,ic] = unique(bMerged,'rows');
            bMerged = bMerged(ia,:);
            for jj = 1:length(fnames)
                [~,~,uni_full] = unique(bMerged(:,1));
                f = bMerged(:,end) == jj;
%                 [C,ia,ic] = unique(uni_full(f,1));
                counts = accumarray(uni_full(f,1),1);
                data.(fnames{jj}).bStart = mat2cell(bMerged(f,2)',1,counts);
                data.(fnames{jj}).bEnd = mat2cell(bMerged(f,3)',1,counts);
            end
            
            obj.data_salam = data;
            fprintf('\t----%s - salam compiled \n',obj.get_full_experiment)
        end
        
        %% Compile Jaaba Data
        function obj = compile_jaaba(obj)
            if isempty(obj.raw_jaaba)
                fprintf('\t\t ...no jaaba data currently loaded\n')
                data = [];
                return
            end
            
            temp = obj.raw_jaaba;
            timestamps = temp.timestamps;

            clear jaaba
            data.aniID = [temp.trx.id];
            [~,~,uniID] = unique(data.aniID,'stable');
            data.uniID = uniID';
            data.tStart = timestamps(temp.allScores.tStart);
            data.tEnd = timestamps(temp.allScores.tEnd);
            data.(temp.behaviorName).bStart = cellfun(@(x) timestamps(x), temp.allScores.t0s, 'UniformOutput',false);
            data.(temp.behaviorName).bEnd = cellfun(@(x) timestamps(x), temp.allScores.t1s, 'UniformOutput',false);
            data.timestamps = timestamps;
            
            obj.data_jaaba = data;
            fprintf('\t----%s - jaaba compiled \n',obj.get_full_experiment)
        end
        
        %% Compile JB Data
        function obj = compile_jb(obj)
            if isempty(obj.raw_jb)
                fprintf('\t\t ...no jb data currently loaded\n')
                data = [];
                return
            end
            
            trx = obj.raw_jb.trx;
            
            data = struct();
            data.aniID = [trx.numero_larva_num];
            [~,~,uniID] = unique(data.aniID,'stable');
            data.uniID = uniID';
            data.tStart = arrayfun(@(x) min(x.t), trx)';
            data.tEnd = arrayfun(@(x) max(x.t), trx)';

            beh_name = {'run','cast','stop','hunch','back','roll','beh7','beh8','beh9','beh10','beh11','beh12'};
            beh_type = {'t_start_stop',...
                't_start_stop_large',...
                't_start_stop_large_small'};
            beh_type = {'t_start_stop'};

            tokens = regexp(beh_type,'(^\w)|[_](\w)','tokens');
            tokens = cellfun(@(x) string([x{:}]), tokens, 'UniformOutput',false);
            beh_code = cellfun(@(x) strcat(x{:}),tokens, 'UniformOutput', false);

            for ii = 1:length(beh_type)
                len = length(trx(1).(beh_type{ii}));
                beh_cell = [trx.(beh_type{ii})];
                beh_cell = cellfun(@transpose, beh_cell,'UniformOutput', false);
                beh_starts = cellfun(@(x) x(1,:),beh_cell,'UniformOutput',false, 'ErrorHandler',@(S,varargin) []);
                beh_starts = reshape(beh_starts,len,[]);
                beh_ends = cellfun(@(x) x(2,:),beh_cell,'UniformOutput',false, 'ErrorHandler',@(S,varargin) []);
                beh_ends = reshape(beh_ends,len,[]);
                for jj = 1:len
                    data.([beh_name{jj},'_',beh_code{ii}]).bStart = beh_starts(jj,:);
                    data.([beh_name{jj},'_',beh_code{ii}]).bEnd = beh_ends(jj,:);
                end
            end
            
            obj.data_jb = data;
            fprintf('\t----%s - jb compiled \n',obj.get_full_experiment)
        end
        %%
        function obj = compile_choreography(obj)
            if isempty(obj.raw_chore)
                fprintf('\t\t ...no choreography data currently loaded\n')
                return
            end

            [C,~,idx] = unique(obj.raw_chore.id);
            unique_ids = num2cell(C);
            parse_struct = struct('aniID',unique_ids);

            counts = accumarray(idx,1);
            fnames = fieldnames(obj.raw_chore);
            fnames = fnames(fnames~="id");
            for ii = 1:length(fnames)
                fields_cells = mat2cell(obj.raw_chore.(fnames{ii})',1,counts);
                [parse_struct.(fnames{ii})] = fields_cells{:};
            end

            obj.data_chore = parse_struct;
            fprintf('\t----%s - choreography compiled \n',obj.get_full_experiment)
        end
        %% Compile All Data
        function obj = compile_all(obj)
            obj = obj.compile_salam();
            obj = obj.compile_jaaba();
            obj = obj.compile_jb();
            obj = obj.compile_choreography();
        end
        
        %% Clear Raw Chore
        function obj = clear_raw(obj)
            obj.raw_blobs = [];
            obj.raw_chore = [];
            obj.raw_jaaba = [];
            obj.raw_jb = [];
            obj.raw_salam = [];
        end
        
        %%
        function behaviour_data = get_behaviour_data(obj,behaviour)
%             if any(size(obj) > 1)
%                 fprintf('\t please specify a single dataStruct to get behaviour data\n')
%                 return
%             end
            temp = obj;
            for ii = 1:length(temp)
                fields = fieldnames(temp(ii));
                subfields = fields( startsWith(fields,'data_') );
                idx = find( cellfun(@(x) isfield(temp(ii).(x),behaviour), subfields),1 );
                if isempty(idx)
                    fprintf('\t could not find specified behaviour: "%s"\n',behaviour)
                    return
                end
                subfields = subfields{idx};

                fields = fieldnames(temp(ii).(subfields));
                filt = ~structfun(@isstruct,temp(ii).(subfields)) | strcmp(fields,behaviour);
                behaviour_data(ii,1) = rmfield(temp(ii).(subfields),fields(~filt));
            end
            
            [behaviour_data(:).behaviour_name] = deal(behaviour);    
            [behaviour_data(:).behaviour] = behaviour_data.(behaviour);
            behaviour_data = rmfield(behaviour_data,behaviour);
        end
        
        %% Sort By Genotype
        function obj = sort_by_genotype(obj,control_genotype)
            genos = strcat(obj.get_full_genotype','@',obj.get_full_protocol');
            
            [obj(:).iscontrol] = deal(false);
            if nargin > 1
                control_idx = ~cellfun(@isempty,regexp(genos,control_genotype));
                [obj(control_idx).iscontrol] = deal(true);
            end
            
            [C,ia,ic] = unique(genos,'stable');
            groups = num2cell(ic,2);

            [obj(:).group] = groups{:};
        end
    end
    
    %% METHODS - STATIC
    
    methods (Static)
        
        %% -------------------------------------------------------
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
        
        %% 
        function structs = make_dataStructs(file_structure)
            fullpaths = fullfile({file_structure.folder},{file_structure.name});

            [file_details, err_files] = dataStruct.parse_filepaths(fullpaths); 
            fullpaths = fullpaths(setdiff(1:length(fullpaths),err_files));
            timestamps = strcat([file_details.date]','@',[file_details.time]');
            [C,~,ic] = unique(timestamps);
            
            for ii = 1:length(C)
                f = ic == ii;
                fullpath = fullpaths(f);
                try
                    structs(ii) = dataStruct(fullpath);
                catch
                    structs(ii) = dataStruct();
                end
            end
        end
        
        %% -------------------------------------------------------
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
                '[\_](?<time>\d\d\d\d\d\d)']; 
            file_details = regexp(string(filepath_full),expr,'names');
            
            if length(file_details) > 1
                error_files = find(cellfun(@isempty,file_details));
                file_details = vertcat(file_details{:});
            else
                error_files = [];
            end
        end
        
        %% 
        function unique_id = get_unique_ids(animal_ids)
            unique_id = cumsum(diff([animal_ids])~=0);
        end
        
    end
end

