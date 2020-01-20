%% blobClass

% Other m-files required: 
% Subfunctions: 
% MAT-files required: 
% See also:

% Author: Alastair Garner
% email: alastairgarner@outlook.com
% Website: https://github.com/alastairgarner/
% Sep 2019; Last revision: 

classdef blobClass
    %%
    properties
        aniID
        frames
        mmperpx = 0.088
        time
        
        fMin
        fMax
        
        x_pt
        y_pt
        npts
        outline
    end
    
    %%
    methods
        %%
        function obj = blobClass(blobData)
            temp = blobData;
            
            obj.fMin = min(temp.frame);
            obj.fMax = max(temp.frame);
            obj.aniID = temp.aniID;
            obj.frames = temp.frame;
            obj.time = temp.time;
            obj.x_pt = temp.x_cont;
            obj.y_pt = temp.y_cont;
            obj.npts = temp.npts;
            obj.outline = temp.outline;

        end  
        
        %%
        function blobSelect = get_animals(obj,animal_ids)
            ids = [obj.aniID];
            ids = unique(ids);
            filt = any(animal_ids' == ids,1);
            blobSelect = obj(filt);
        end
    end
    
    %%
    methods (Static)
        function obj = blobClass_batch(blobData)
            temp = blobData;

            aniIDs = temp.aniID;
            counts = accumarray(aniIDs',1);
            counts = counts(find(counts));

            fnames = fieldnames(temp);
            fcells = struct2cell(temp);
            fcells = cellfun(@(x) mat2cell(x,1,counts), fcells, 'UniformOutput', false);
            fcells = vertcat(fcells{:});
            temp = cell2struct(fcells,fnames,1);
            
            for ii = 1:length(temp)
                obj(ii) = blobClass(temp(ii));
            end
        end
    end
    
    
end