%% merge structure fields

% Other m-files required: 
% Subfunctions: 
% MAT-files required: 
% See also:

% Author: Alastair Garner
% email: alastairgarner@outlook.com
% Website: https://github.com/alastairgarner/
% October 2019; Last revision:

function merged_struct = merge_structure_fields(old_struct)

temp = old_struct;

fields = fieldnames(temp);
for ii = 1:length(fields)
    field_data = {temp.(fields{ii})};
    if size(field_data{1},1) > 1
        field_data = cellfun(@transpose, field_data,'UniformOutput',false);
    end
    if isstr(field_data{1})
        if ~length(unique(field_data))
            fprintf('\t "%s" field has inconsistent values - please check! \n',fields{ii})
            return
        end
        merged_struct.(fields{ii}) = string(field_data);
    elseif isnumeric(field_data{1})
        merged_struct.(fields{ii}) = [field_data{:}];
    elseif iscell(field_data{1})
        merged_struct.(fields{ii}) = [field_data{:}];
    elseif isstruct(field_data{1})
        nested_struct = vertcat(temp.(fields{ii}));
        nested_merged = merge_structure_fields(nested_struct);
        merged_struct.(fields{ii}) = nested_merged;
    elseif islogical(field_data{1})
        merged_struct.(fields{ii}) = [field_data{:}];
    else
        fprintf('\t "%s" field could not be merged! \n',fields{ii})
        merged_struct.(fields{ii}) = [];
    end
end




