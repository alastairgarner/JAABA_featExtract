%% flexi_legend

% Other m-files required: 
% Subfunctions: 
% MAT-files required: 
% See also:

% Author: Alastair Garner
% email: alastairgarner@outlook.com
% Website: https://github.com/alastairgarner/
% Sep 2019; Last revision: 

%%
function [leg,icns,plts,txt] = flexi_legend(varargin)
    
    f = strcmpi(varargin,'RelativePosition');
    relPos = varargin{find(f)+1};
    f(find(f)+1) = true;
    args = varargin(~f);

%     pbp = plotboxpos;
%     sz = pbp(3:end);
%     topright = sum(reshape(pbp,2,[])',1);
%     
% %     annotation('rectangle',pbp,'Color','red');
%     prop = [sz,sz].*relPos;
%     legpos = prop+[pbp(1:2),0,0];
%     args = [args,{'Position'},{legpos}];

    [leg,icns,plts,txt] = legend(args{:});
    
    lines = findobj(icns,'Type','Line');
    odds = [1:2:numel(lines)];
    xdata = get(lines(1),'XData');
    xdata = max(xdata)+[-.1 0];
    set(lines,'XData',xdata);

end