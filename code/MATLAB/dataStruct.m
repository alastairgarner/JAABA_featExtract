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
        exp_date
        exp_time
        exp_driver
        exp_effector
        exp_protocol1
        exp_protocol2
        exp_protocol3
        exp_protocol4
        rig
        
        files_blobs
        files_choreography
        files_featextract
        files_jaaba
        files_jb
        
        data_blobs
        data_choreography
        data_featextract
        data_jaaba
        data_jb
        
    end
    
    %% METHODS NORMAL
    methods
        
        %% Construct object from file paths
        function obj = dataStruct(file_paths)
            %  CONSTRUCTOR
            %    Detailed explanation goes here
            file_structure = dataStruct.parse_filepaths(file_paths(1));
            
            obj.exp_date = file_structure.date;
            obj.exp_time = file_structure.time;
            obj.exp_driver = file_structure.driver;
            obj.exp_effector = file_structure.effector;
            obj.exp_protocol1 = file_structure.protocol1;
            obj.exp_protocol2 = file_structure.protocol2;
            obj.exp_protocol3 = file_structure.protocol3;
            obj.exp_protocol4 = file_structure.protocol4;
            obj.rig = file_structure.rig;
            
            
            [fpath,fname,fext] = cellfun(@fileparts, file_paths, 'uni', false);
    
            filt_blob = startsWith(fext,'.blo');
            filt_Chore_dat = startsWith(fext,'.dat');
            filt_Chore_txt = endsWith(fname,'compiledChore') & startsWith(fext,'.tx');
            filt_salam = ~endsWith(fname,'compiledChore') & startsWith(fext,'.tx');
            filt_jaaba = contains(fpath,'jaaba') & startsWith(fext,'.mat');
            filt_jb = contains(fpath,'jb') & startsWith(fext,'.mat');
%             if any(filt_jb)
%                 fpath'
%                 filt_jb
%             end
            
            obj.files_blobs = string(file_paths(filt_blob)');
            obj.files_choreography = string(file_paths(filt_Chore_txt)');
            obj.files_featextract = string(file_paths(filt_salam)');
            obj.files_jaaba = string(file_paths(filt_jaaba)');
            obj.files_jb = string(file_paths(filt_jb)');
        end

        %% Load blobs data
        function obj = load_blobsdata(obj)
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

            datfull = sortrows(datfull,[1,2]);
            datcell = mat2cell(datfull',ones(size(datfull,2),1),size(datfull,1));
            obj.data_blobs = cell2struct(datcell,columns);
            obj.data_blobs.outline = outlinesfull';
            
            fprintf('----%s - blobs loaded \n',obj.get_full_experiment)
        end
        
        %% Load choreography data
        function obj = load_choreographydata(obj)
            filepath = obj.files_choreography;
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

            obj.data_choreography.id = dataArray{1};
            obj.data_choreography.et = dataArray{2};
            for jj = 1:length(dataArray)-2
                field_name = chores{jj};
                obj.data_choreography.(field_name) = dataArray{jj+2};
            end
            fprintf('----%s - choreography loaded \n',obj.get_full_experiment)
        end
        
        %% Load feature extraction (salam) data
        function obj = load_salamdata(obj)
            delimiter = '\t';
            startRow = 1;
            formatSpec = '%s%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
            
            filenames = obj.files_featextract;
            
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
                obj.data_featextract.(behaviour(ii)) = salam_matrix;
            end
            fprintf('----%s - salam loaded \n',obj.get_full_experiment)
        end
        
        %% load jaaba data
        function obj = load_jaabadata(obj)
            for ii = 1:length(obj.files_jaaba)
                data = load(obj.files_jaaba(ii));
                f = fieldnames(data);
                for jj = 1:length(f)
                   obj.data_jaaba.(f{jj}) = data.(f{jj});
                end
            end
            fprintf('----%s - jaaba loaded \n',obj.get_full_experiment)
        end
        
        %% load jaaba data
        function obj = load_jbdata(obj)
            for ii = 1:length(obj.files_jb)
                data = load(obj.files_jb(ii));
                f = fieldnames(data);
                for jj = 1:length(f)
                   obj.data_jb.(f{jj}) = data.(f{jj});
                end
            end
            fprintf('----%s - jb loaded \n',obj.get_full_experiment)
        end
        
        %% extract full genotype
        function full_genotype = get_full_genotype(obj)
            full_genotype = strcat(obj.exp_driver,'@',obj.exp_effector);
        end

        %% -------------------------------------------------------
        function full_protocol = get_full_protocol(obj)
            full_protocol = strcat('#',obj.exp_protocol1,...
                '#',obj.exp_protocol2,...
                '#',obj.exp_protocol3,...
                '#',obj.exp_protocol4);
        end

        %% -------------------------------------------------------
        function full_timestamp = get_full_timestamp(obj)
            full_timestamp = strcat(obj.exp_date,'_',obj.exp_time);
        end  
        
        %% -------------------------------------------------------
        function full_experiment = get_full_experiment(obj)
            full_experiment = strcat(obj.get_full_timestamp,'@',...
                obj.get_full_genotype,'@',...
                obj.get_full_protocol,'@');
        end  
        
        %% -------------------------------------------------------
        function salam_data = compile_salam(obj)
            fnames = fieldnames(obj.data_featextract);
            [rows,cols] = structfun(@size, obj.data_featextract);
            behID = repelem([1:length(fnames)]',rows);
            dat = struct2cell(obj.data_featextract);
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

                if size(temp,1) == 1
                    tS = [nan];
                    tE = [nan];
                    bS = [nan];
                    bE = [nan];
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
                end

                tStart = [tStart tS];
                tEnd = [tEnd tE];
                bStart = [bStart; bS];
                bEnd = [bEnd; bE];
                bType = [bType; temp(:,end)];
            end
            
            salam_data.aniID = uniID';
            salam_data.tStart = tStart;
            salam_data.tEnd = tEnd;

            bMerged = [aniID, bStart, bEnd, bType];
            [C,ia,ic] = unique(bMerged,'rows');
            bMerged = bMerged(ia,:);
            for jj = 1:length(fnames)
                f = bMerged(:,end) == jj;
                [C,ia,ic] = unique(bMerged(f,1));
                counts = accumarray(ic,1);
                salam_data.(fnames{jj}).bStart = mat2cell(bMerged(f,2),counts,1)';
                salam_data.(fnames{jj}).bEnd = mat2cell(bMerged(f,3),counts,1)';
            end
        end
        
        %% -------------------------------------------------------
        function jaaba = compile_jaabadata(obj,fields_wanted)
            jaaba = struct();
            for ii = 1:length(fields_wanted)
                args = cellstr(fields_wanted{ii});
                len = length(args);
                try
                    if len == 1
                        jaaba.(args{end}) = [obj.data_jaaba.(args{1})];
                    elseif len == 2
                        jaaba.(args{end}) = [obj.data_jaaba.(args{1}).(args{2})];
                    elseif len == 3
                        jaaba.(args{end}) = [obj.data_jaaba.(args{1}).(args{2}).(args{3})];
                    end
                catch
                    fprintf('... No field called jaaba_data.%s \n',strjoin(strcat(args,'.')));
                end
            end
            
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

            file_details = dataStruct.parse_filepaths(fullpaths);           
            timestamps = strcat([file_details.date]','@',[file_details.time]');
            [C,~,ic] = unique(timestamps);
            
            for ii = 1:length(C)
                f = ic == ii;
                fullpath = fullpaths(f);                
                structs(ii) = dataStruct(fullpath);
            end
        end
        
        %% -------------------------------------------------------
        function file_details = parse_filepaths(filepath_full)
            
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
                file_details = vertcat(file_details{:});
            end
        end
        
    end
end

