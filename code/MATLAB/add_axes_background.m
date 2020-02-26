%% add_axes_background

% Other m-files required: 
% Subfunctions: 
% MAT-files required: 
% See also:

% Author: Alastair Garner
% email: alastairgarner@outlook.com
% Website: https://github.com/alastairgarner/
% Sep 2019; Last revision: 

%%
function add_axes_background(axHandle,BackgroundColor,LineColor,LineWidth,FaceAlpha)

ax = axHandle;

xlm = get(ax,'XLim');
ylm = get(ax,'YLim');
xtk = get(ax,'XTick');
ytk = get(ax,'YTick');

xp = [xlm fliplr(xlm)];
yp = repelem(ylm,2);

xtk = repmat(setdiff(xtk,xlm),2,1);
ytk = repmat(setdiff(ytk,ylm),2,1);

xl = [xtk,repmat(xlm',1,size(ytk,2))];
yl = [repmat(ylm',1,size(xtk,2)),ytk];

hold(ax,'on')
patch(ax,xp,yp,BackgroundColor,'FaceAlpha',FaceAlpha,'EdgeColor','none','DisplayName','background_patch');
plot(ax,xl,yl,'Color',LineColor,'LineWidth',LineWidth,'DisplayName','background_line');
hold(ax,'off')

background_elements = findobj(ax.Children,'-regexp','DisplayName','background*');
other_elements = findobj(ax.Children,'-not',{'-regexp','DisplayName','background*'});
set(ax,'Children',[other_elements;background_elements],...
    'XLim',xlm,...
    'YLim',ylm);

end