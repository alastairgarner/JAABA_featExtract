%% master.m

% Other m-files required: 
% Subfunctions: 
% MAT-files required: 
% See also:

% Author: Alastair Garner
% email: alastairgarner@outlook.com
% Website: https://github.com/alastairgarner/
% September 2019; Last revision: 

function merged = merge_structarray(structs)

%structs = temp;

merged = struct();
for ii = 1:size(structs,1)
    temp = structs(ii);
    if ii==1
        merged = temp;
        continue
    end
    
    lev1 = fieldnames(temp);
%     goodfields = setdiff(lev1,ignore_fields);
    for jj = 1:length(lev1)
        if ~isstruct(temp.(lev1{jj}))
            merged.(lev1{jj}) = [merged.(lev1{jj}) temp.(lev1{jj})];
            
        else
            temp2 = temp.(lev1{jj});
            lev2 = fieldnames(temp2);
            for kk = 1:length(lev2)
                if ~isstruct(temp2.(lev2{kk}))
                    merged.(lev1{jj}).(lev2{kk}) = [temp2.(lev2{kk}) temp2.(lev2{kk})];
                else
                    temp3 = temp2.(lev2{kk});
                    lev3 = fieldnames(temp3);
                    for pp = 1:length(lev3)
                        if ~isstruct(temp3.(lev3{pp}))
                            merged.(lev1{jj}).(lev2{kk}).(lev3{pp}) = [temp3.(lev3{pp}) temp3.(lev3{pp})];
                        else

                        end
                    end
                    
                end
            end
        end
    end
end
        
        


