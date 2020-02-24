%% update_contents_file

% Other m-files required: 
% Subfunctions: 
% MAT-files required: 
% See also:

% Author: Alastair Garner
% email: alastairgarner@outlook.com
% Website: https://github.com/alastairgarner/
% Sep 2019; Last revision: 

%%
function update_contents_file(contentfile,newfilenames,pipelines,params)
    timestamps = regexp(newfilenames,'(\d{8}_\d{6})', 'tokens', 'once');
    if isempty(timestamps)
        return
    end

    if ~isfile(fullfile(params.directories.data_compiled,contentfile))
        dataClass.generate_contents_file(contentfile,pipelines,params)
    else
        timestamps = regexp(newfilenames,'(\d{8}_\d{6})', 'tokens', 'once');
        timestamps = string([timestamps{:}]);
        [timestamps_sorted,I] = sort(timestamps);
        newfilenames = newfilenames(I);

        filebits = split(newfilenames,filesep,1);
        file_pipelines = squeeze(filebits(end-1,:,:));
        content_matrix = file_pipelines == pipelines;

        [C,ia,ic] = unique(timestamps_sorted);
        filenames_sorted = newfilenames(ia)';
        counts = accumarray(ic,1);
        content_cells = mat2cell(content_matrix,counts,size(content_matrix,2));
        content_out = cellfun(@(x) any(x,1), content_cells,'UniformOutput', false);
        content_out = vertcat(content_out{:});

        [file_details,error_files] = dataClass.parse_filepaths(filenames_sorted);
        details_cell = struct2cell(file_details)';
        details_cell = details_cell(:,[1,4,3,5,6,7,8,9]);
        details_cell = cellfun(@char, details_cell, 'UniformOutput', false);

        formatSpec = '%d%d%d%d%d%s%s%s%s%s%s%s%s%[^\n\r]';
        delimiter = ",";
        fid = fopen(fullfile(params.directories.data_compiled,contentfile));
        csv = textscan(fid, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'ReturnOnError', false);
        fclose(fid);

        cont_mat = [csv{1:numel(pipelines)}];
        match = all(ismember([csv{numel(pipelines)+1:end-1}],details_cell),2);
        if any(match)
            cont_mat(match,:) = cont_mat(match,:) | content_out;
            out_cell = [num2cell(cont_mat),csv{6:end-1}]';
        else
            csv(1:numel(pipelines)) = cellfun(@num2cell, csv(1:numel(pipelines)), 'UniformOutput', false);
            out_cell = [csv{1:end-1}];
            out_cell = vertcat(out_cell,[num2cell(content_out),details_cell]);
            out_cell = sortrows(out_cell,[6,7])';
        end

        fid = fopen(fullfile(params.directories.data_compiled,contentfile),"w");
        fprintf(fid,"%d,%d,%d,%d,%d,%s,%s,%s,%s,%s,%s,%s,%s\n",out_cell{:});
        fclose(fid);
    end

end
