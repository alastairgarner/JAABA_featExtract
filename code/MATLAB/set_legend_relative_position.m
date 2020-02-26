%% set_legend_relative_position

% Other m-files required: 
% Subfunctions: 
% MAT-files required: 
% See also:

% Author: Alastair Garner
% email: alastairgarner@outlook.com
% Website: https://github.com/alastairgarner/
% Sep 2019; Last revision: 

%%
function set_legend_relative_position(RelativePosition)

iconscaler = [3;1];

legs = findobj(gcf,'Type','Legend');
for xx = 1:numel(legs)
    ax = legs(xx).Axes;
    
    % Change line length
    IconLength = get(legs(xx),'ItemTokenSize');
    set(legs(xx),'ItemTokenSize',ceil(IconLength./iconscaler));
    
    % Change position
    pbp = plotboxpos(ax);
    % annotation('rectangle',pos,'Color','red');
    sz = pbp(3:end);
    topright = sum(reshape(pbp,2,[])',1);
    prop = [sz,sz].*RelativePosition;
    legpos = prop+[pbp(1:2),0,0];
    legs(xx).Position = legpos;
end

% ax = findobj(gcf,'Type','Axes');
% for xx = 1:numel(ax)
%     leg = legend(ax(xx));
%     pbp = plotboxpos(ax(xx));
%     % annotation('rectangle',pos,'Color','red');
%     sz = pbp(3:end);
%     topright = sum(reshape(pbp,2,[])',1);
%     prop = [sz,sz].*RelativePosition;
%     legpos = prop+[pbp(1:2),0,0];
%     leg.Position = legpos;
% end

end