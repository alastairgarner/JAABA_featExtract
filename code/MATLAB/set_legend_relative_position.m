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
function set_legend_relative_position(anchor,RelativePosition)

iconscaler = [3;1];

% anchor = 'bottomleft';

% legs = findobj(gcf,'Type','Legend');
% for xx = 1:numel(legs)
%     ax = legs(xx).Axes;
%     
%     % Change line length
%     IconLength = get(legs(xx),'ItemTokenSize');
%     set(legs(xx),'ItemTokenSize',ceil(IconLength./iconscaler));
%     
%     % Change position
%     pbp = plotboxpos(ax);
%     % annotation('rectangle',pos,'Color','red');
%     sz = pbp(3:end);
%     topright = sum(reshape(pbp,2,[])',1);
%     prop = [sz,sz].*RelativePosition;
%     legpos = prop+[pbp(1:2),0,0];
%     legpos(3:4) = max([legpos(3:4); legs(xx).Position(3:4)],[],1);
%     legs(xx).Position = legpos;
% end


iconscaler = [3;1];

legs = findobj(gcf,'Type','Legend');
for xx = 1:numel(legs)
    ax = legs(xx).Axes;
    
    % Change line length
    IconLength = get(legs(xx),'ItemTokenSize');
    set(legs(xx),'ItemTokenSize',ceil(IconLength./iconscaler));
    
    % Change position
    pbp = plotboxpos(ax);
    sz = pbp(3:end);
    prop = sz.*RelativePosition(1:2);
    anchorPos = prop+pbp(1:2);
    
    startPos = legs(xx).Position;
    startCorners = startPos + [0,0,startPos(1:2)];
    switch anchor
        case 'bottomleft'
            posDiff = anchorPos - startCorners([1,2]);
            newLegPos = startPos + [posDiff,0,0];
        case 'bottomright'
            posDiff = anchorPos - startCorners([3,2]);
            newLegPos = startPos + [posDiff,0,0];
        case 'topleft'
            posDiff = anchorPos - startCorners([1,4]);
            newLegPos = startPos + [posDiff,0,0];
        case 'topright'
            posDiff = anchorPos - startCorners([3,4]);
            newLegPos = startPos + [posDiff,0,0];
    end
    legs(xx).Position = newLegPos;
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