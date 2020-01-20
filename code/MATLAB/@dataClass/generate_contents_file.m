%% generate_contents_file

% Other m-files required: 
% Subfunctions: 
% MAT-files required: 
% See also:

% Author: Alastair Garner
% email: alastairgarner@outlook.com
% Website: https://github.com/alastairgarner/
% Sep 2019; Last revision: 

%%
function generate_contents_file(contentfile,pipelines,params)
    % contentfile = "contents.csv";
    searchterm = ["mwt","choreography","salam","jaaba","jb"];
    searchterm = pipelines;

    d = dir(fullfile(params.directories.data_compiled,"**","*.mat"));
    if numel(d) == 0
        return
    end
    [~,pipelines] = cellfun(@fileparts, {d.folder}, 'UniformOutput', false);
    pipelines = string(pipelines)';
    content_matrix = pipelines == searchterm;

    fullpaths = fullfile({d.folder},{d.name});
    timestamps = regexp(fullpaths,'(\d{8}_\d{6})', 'tokens', 'once');
    timestamps = string([timestamps{:}]);
    [timestamps_sorted,I] = sort(timestamps);
    fullpaths = fullpaths(I);

    content_matrix = content_matrix(I,:);

    [C,ia,ic] = unique(timestamps_sorted);
    fullpaths_sorted = fullpaths(ia)';
    counts = accumarray(ic,1);
    content_cells = mat2cell(content_matrix,counts,size(content_matrix,2));
    content_out = cellfun(@(x) any(x,1), content_cells,'UniformOutput', false);
    content_out = vertcat(content_out{:});

    [file_details,error_files] = dataClass.parse_filepaths(fullpaths_sorted);
    details_cell = struct2cell(file_details)';
    details_cell = details_cell(:,[1,4,3,5,6,7,8,9]);

    out_cell = [num2cell(content_out),details_cell]';

    fid = fopen(fullfile(params.directories.data_compiled,contentfile),"w");
    fprintf(fid,"%d,%d,%d,%d,%d,%s,%s,%s,%s,%s,%s,%s,%s\n",out_cell{:});
    fclose(fid);

end
