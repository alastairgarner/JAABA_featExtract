%% larvaPlot

% Other m-files required: 
% Subfunctions: 
% MAT-files required: 
% See also:

% Author: Alastair Garner
% email: alastairgarner@outlook.com
% Website: https://github.com/alastairgarner/
% Sep 2019; Last revision: 

classdef larvaPlot
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        date
        time
        driver
        effector
        protocol1
        protocol2
        protocol3
        protocol4
        rig
        
        data
        
        mwt
        chore
        salam
        jaaba
        jb
        
    end
    
    methods        
        %%
        function obj = larvaPlot(varargin)
            p = inputParser;
            
            addOptional(p,'dataStruct',[]);
            
            parse(p,varargin{:})
            
            if isa(p.Results.dataStruct,'dataStruct')
                fprintf('\n Good Struct \n');
                dS = p.Results.dataStruct;
                obj.date = dS.exp_date;
                obj.time = dS.exp_time;
                obj.driver = dS.exp_driver;
                obj.effector = dS.exp_effector;
                obj.protocol1 = dS.exp_protocol1;
                obj.protocol2 = dS.exp_protocol2;
                obj.protocol3 = dS.exp_protocol3;
                obj.protocol4 = dS.exp_protocol4;
                obj.rig = dS.rig;
                
%                 obj.salam = larvaPlot.compile_salam(dS);
%                 obj.jaaba = larvaPlot.compile_jaaba(dS,fields);
                obj = compile_alldata(obj,dS);
                
            elseif ~isa(p.Results.dataStruct,'dataStruct')
                fprintf('\n Bad Struct \n');
            end
        end
        
        %%
        function obj = compile_alldata(obj,dStruct)
            obj.salam = larvaPlot.compile_salam(dStruct);
            obj.jaaba = larvaPlot.compile_jaaba(dStruct);
            obj.jb = larvaPlot.compile_jb(dStruct);
        end
        
        %%
        function [ts_tracked,uniIDs] = get_tracked(obj,pipeline,behaviour)
            temp = obj.(pipeline{:});
            bins = 0:.2:max(ceil(temp.tEnd));
            flt = ~isnan(temp.tStart);

            aniIDs = temp.aniID(flt);
            uniIDs = 1:length(aniIDs);
            [binStart,~] = discretize(temp.tStart(flt),bins);
            [binEnd,~] = discretize(temp.tEnd(flt),bins);

            spStart = sparse(uniIDs',binStart',[1],max(uniIDs),length(bins));
            spEnd = sparse(uniIDs',binEnd',[-1],max(uniIDs),length(bins));
            spBoth = cumsum((spStart+spEnd),2);
            ts_tracked = full(sum(spBoth,1));
        end
        
        %%
        function [ts_behaved,uniIDs] = get_behaved(obj,pipeline,behaviour)
            temp = obj.(pipeline{:});
            bins = 0:.2:max(ceil(temp.tEnd));
            flt = ~isnan(temp.tStart);

            aniIDs = temp.aniID(flt);
            uniIDs = 1:length(aniIDs);

            temp = obj.(pipeline{:}).(behaviour{:});
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
            ts_behaved = full(sum(spBoth,1));
        end
        
    end
    
    %%
    methods (Static)
        
        
        %%
        function salam_data = compile_salam(dataStruct)
            if isempty(dataStruct.data_featextract)
                fprintf('\t\t ...no salam data \n')
                salam_data = [];
                return
            end
            
            fnames = fieldnames(dataStruct.data_featextract);
            [rows,cols] = structfun(@size, dataStruct.data_featextract);
            behID = repelem([1:length(fnames)]',rows);
            dat = struct2cell(dataStruct.data_featextract);
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
            
            salam_data.pipeline = 'salam';
            salam_data.aniID = uniID';
            salam_data.tStart = tStart;
            salam_data.tEnd = tEnd;
            
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
                salam_data.(fnames{jj}).bStart = mat2cell(bMerged(f,2)',1,counts);
                salam_data.(fnames{jj}).bEnd = mat2cell(bMerged(f,3)',1,counts);
            end
        end
        
        %% 
        function jaaba = compile_jaaba(dataStruct)
%             get_fields = {["behaviorName"];
%                     ["timestamps"];
%                     ["trx","id"];
%                     ["allScores","tStart"];
%                     ["allScores","tEnd"];
%                     ["allScores","t0s"];
%                     ["allScores","t1s"];
%                     ["allScores","postprocessed"]};
%             new_fields = {[""]};
%                 
%             jaaba = struct();
%             for ii = 1:length(get_fields)
%                 args = cellstr(get_fields{ii});
%                 len = length(args);
%                 try
%                     if len == 1
%                         jaaba.(args{end}) = [dataStruct.data_jaaba.(args{1})];
%                     elseif len == 2
%                         jaaba.(args{end}) = [dataStruct.data_jaaba.(args{1}).(args{2})];
%                     elseif len == 3
%                         jaaba.(args{end}) = [dataStruct.data_jaaba.(args{1}).(args{2}).(args{3})];
%                     end
%                 catch
%                     fprintf('... No field called jaaba_data.%s \n',strjoin(strcat(args,'.')));
%                 end
%             end
            if isempty(dataStruct.data_jaaba)
                fprintf('\t\t ...no jaaba data \n')
                jaaba = [];
                return
            end
            
            temp = dataStruct.data_jaaba;
            timestamps = temp.timestamps;

            clear jaaba
            jaaba.aniID = [temp.trx.id];
            jaaba.tStart = timestamps(temp.allScores.tStart);
            jaaba.tEnd = timestamps(temp.allScores.tEnd);
            jaaba.(temp.behaviorName).bStart = cellfun(@(x) timestamps(x), temp.allScores.t0s, 'UniformOutput',false);
            jaaba.(temp.behaviorName).bEnd = cellfun(@(x) timestamps(x), temp.allScores.t1s, 'UniformOutput',false);
            jaaba.timestamps = timestamps;
        end
        
        %% compile jb
        function jb = compile_jb(dataStruct)
            if isempty(dataStruct.data_jb)
                fprintf('\t\t ...no jb data \n')
                jb = [];
                return
            end
            
            trx = dataStruct.data_jb.trx;
            
            jb = struct();
            jb.aniID = [trx.numero_larva_num];
            jb.tStart = arrayfun(@(x) min(x.t), trx)';
            jb.tEnd = arrayfun(@(x) max(x.t), trx)';

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
                    jb.([beh_name{jj},'_',beh_code{ii}]).bStart = beh_starts(jj,:);
                    jb.([beh_name{jj},'_',beh_code{ii}]).bEnd = beh_ends(jj,:);
                end
            end            
        end
        
    end
    
    
end