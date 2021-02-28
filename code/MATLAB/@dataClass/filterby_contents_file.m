%% filterby_contents_file

% Other m-files required: 
% Subfunctions: 
% MAT-files required: 
% See also:

% Author: Alastair Garner
% email: alastairgarner@outlook.com
% Website: https://github.com/alastairgarner/
% Sep 2019; Last revision: 

%%
function filelist_filtered = filterby_contents_file(filelist,contentfile,pipelines, params)
    % get filelist for data files
    if ~isfile(fullfile(params.directories.data_compiled,contentfile))
        dataClass.generate_contents_file(contentfile,pipelines,params);
        
        if ~isfile(fullfile(params.directories.data_compiled,contentfile))
            filelist_filtered = filelist;
            return
        end
    end
    
    [types,idx] = dataClass.get_filetype(filelist);
    filelist = filelist(idx);
    [filedeets,err] = dataClass.parse_filepaths(filelist);
    filelist = filelist(~err);
    filedeets = rmfield(filedeets,["rig"]);
    filedeets = orderfields(filedeets,{'date','time','driver','effector',...
        'protocol1','protocol2','protocol3','protocol4'});
    filesprocess = [cellstr(types(~err))',cellstr(struct2cell(filedeets)')];

    % get previously processed files
    formatSpec = '%d%d%d%d%d%s%s%s%s%s%s%s%s%[^\n\r]';
    delimiter = ",";
    fid = fopen(fullfile(params.directories.data_compiled,contentfile));
    csv = textscan(fid, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'ReturnOnError', false);
    fclose(fid);

    mat = [csv{1:5}];
    counts = sum(mat,2);
    [r,c] = find(mat');
    donefiles = [cellstr(pipelines(r))',repelem([csv{6:end-1}],counts,1)];

    % compare lists to generate logical filter
%     filter = ~all(ismember(filesprocess,donefiles),2);
    filter = ~ismember(...
        string(filesprocess(:,[1,2,3])),...
        string(donefiles(:,[1,2,4])),...
        'rows');
    
    filelist_filtered = filelist(filter);

end
