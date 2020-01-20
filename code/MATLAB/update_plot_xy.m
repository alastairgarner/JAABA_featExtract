%% update_plot_xy

% Other m-files required: 
% Subfunctions: 
% MAT-files required: 
% See also:

% Author: Alastair Garner
% email: alastairgarner@outlook.com
% Website: https://github.com/alastairgarner/
% October 2019; Last revision:

function [x,y] = update_plot_xy(linehandle,blobsStruct,unique_ids,frame)

[x,y] = unpackKerr_xy_batch(blobsStruct(unique_ids).outline(frame),...
            blobsStruct(unique_ids).npts(frame),...
            blobsStruct(unique_ids).x_cent(frame),...
            blobsStruct(unique_ids).y_cent(frame));
        
set(linehandle,'XData',x{:},'YData',y{:});