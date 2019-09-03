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
        
    end
    
    %%
    methods (Static)
        
        
        %%
        function salam_data = compile_salam(dataStruct)
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
        
        %% 
        function jaaba = compile_jaaba(dataStruct)
            get_fields = {["behaviorName"];
                    ["timestamps"];
                    ["trx","id"];
                    ["allScores","tStart"];
                    ["allScores","tEnd"];
                    ["allScores","t0s"];
                    ["allScores","t1s"];
                    ["allScores","postprocessed"]};
            new_fields = {[""]};
                
            jaaba = struct();
            for ii = 1:length(get_fields)
                args = cellstr(get_fields{ii});
                len = length(args);
                try
                    if len == 1
                        jaaba.(args{end}) = [dataStruct.data_jaaba.(args{1})];
                    elseif len == 2
                        jaaba.(args{end}) = [dataStruct.data_jaaba.(args{1}).(args{2})];
                    elseif len == 3
                        jaaba.(args{end}) = [dataStruct.data_jaaba.(args{1}).(args{2}).(args{3})];
                    end
                catch
                    fprintf('... No field called jaaba_data.%s \n',strjoin(strcat(args,'.')));
                end
            end
        end
        
        %% compile jb
        function jb = compile_jb(dataStruct)
            trx = dataStruct.data_jb.trx;
            
            jb = struct();
            jb.aniID = [trx.numero_larva_num];
            jb.tStart = arrayfun(@(x) min(x.t), trx)';
            jb.tEnd = arrayfun(@(x) max(x.t), trx)';

            beh_name = {'run','cast','stop','hunch','back','roll','beh7','beh8','beh9','beh10','beh11','beh12'};
            beh_type = {'t_start_stop',...
                't_start_stop_large',...
                't_start_stop_large_small'};

            tokens = regexp(beh_type,'(^\w)|[_](\w)','tokens');
            tokens = cellfun(@(x) string([x{:}]), tokens, 'UniformOutput',false);
            beh_code = cellfun(@(x) strcat(x{:}),tokens, 'UniformOutput', false);

            for ii = 1:length(beh_type)
                len = length(trx(1).(beh_type{ii}));
                beh_cell = [trx.(beh_type{ii})];
                beh_starts = cellfun(@(x) x(:,1),beh_cell,'UniformOutput',false, 'ErrorHandler',@(S,varargin) []);
                beh_starts = reshape(beh_starts,len,[]);
                beh_ends = cellfun(@(x) x(:,2),beh_cell,'UniformOutput',false, 'ErrorHandler',@(S,varargin) []);
                beh_ends = reshape(beh_ends,len,[]);
                for jj = 1:len
                    jb.([beh_name{jj},'_',beh_code{ii}]).bStart = beh_starts(jj,:);
                    jb.([beh_name{jj},'_',beh_code{ii}]).bEnd = beh_ends(jj,:);
                end
            end            
        end
        
    end
    
    
end