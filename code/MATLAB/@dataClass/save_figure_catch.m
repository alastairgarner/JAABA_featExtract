%% save_figure_catch

% Other m-files required: 
% Subfunctions: 
% MAT-files required: 
% See also:

% Author: Alastair Garner
% email: alastairgarner@outlook.com
% Website: https://github.com/alastairgarner/
% Sep 2019; Last revision: 

%%
function save_figure_catch(ax,figure_path)

if ~isempty(ax) & ~isempty(figure_path)
    if ~isdir(fileparts(figure_path))
        mkdir(fileparts(figure_path))
    end
    
    set(gcf, 'InvertHardCopy', 'off',...
        'PaperPositionMode', 'auto',...
        'Color', [1 1 1]);
    
    shrink_fig2axes(gcf)
    
    print(gcf,figure_path,'-dsvg','-painters');
%     savefig(figure_path);
%     print(gcf,figure_path,'-dpdf','-painters','-fillpage');
else
    fprintf("no figure to save \n")
end











