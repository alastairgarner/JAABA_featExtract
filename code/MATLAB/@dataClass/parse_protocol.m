%% plot_timeseries

% Other m-files required: 
% Subfunctions: 
% MAT-files required: 
% See also:

% Author: Alastair Garner
% email: alastairgarner@outlook.com
% Website: https://github.com/alastairgarner/
% Sep 2019; Last revision: 

%% 
function protocol_struct = parse_protocol(obj)

    % get stim protocol
    fnames = fieldnames(obj);   
    fnames = fnames(startsWith(fnames,"protocol"));
    protocols = cellfun(@(x) obj.(x), fnames);
    exp = '[_](?<start>\d+)s(?<reps>\d+)x(?<length>\d+)s(?<interval>\d+)';
    prot = regexp(protocols,exp,'names');
    prot = vertcat(prot{:});
    protocol_struct = arrayfun(@(x) structfun(@double,x,'UniformOutput',false), prot);
     
end