%% unpackKerr_xy.m

% Other m-files required: 
% Subfunctions: 
% MAT-files required: 
% See also:

% Author: Alastair Garner
% email: alastairgarner@outlook.com
% Website: https://github.com/alastairgarner/
% August 2019; Last revision: 

function [x,y] = unpackKerr_xy(kerrstring,npts,x_point,y_point)

if isempty(kerrstring)
    x = nan;
    y = nan;
    return
end

ccstr = uint8(kerrstring) - '0';
npt = npts;

bc = reshape([bitand(bitshift(ccstr, -4), 3); bitand(bitshift(ccstr, -2), 3); bitand(ccstr, 3)], 1, []);
bc = bc(1:(npt));

dxx = sum([-1; 1].*(bc == [0;1]), 1);
dyy = sum([-1; 1].*(bc == [2;3]), 1);

x = [0 cumsum(dxx) 0] + x_point;
y = [0 cumsum(dyy) 0] + y_point;