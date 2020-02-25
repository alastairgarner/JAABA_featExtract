%% fix_axes

% Other m-files required: 
% Subfunctions: 
% MAT-files required: 
% See also:

% Author: Alastair Garner
% email: alastairgarner@outlook.com
% Website: https://github.com/alastairgarner/
% Sep 2019; Last revision: 

%%
function fix_axes(handle,font_size,grid_tf)
%     elem = gcf;
    elem = handle;
%     fprintf("rep \n")
    if strcmpi(elem.Type,'figure')
%         fprintf('no colour \n')
%         set(elem, 'Color', 'none')
    end
    
    if strcmpi(elem.Type,'axes')
        try
            [leg,icons,~,~] = legend(elem);
            txts = findobj(icons,'Type','Text');
            set(leg,'FontSize',font_size);
            set(txts,'FontSize',font_size);
        end
        switch grid_tf
            case true
%                 fprintf("axis found \n")
%                 grid(elem,'on')
%                 set(elem, 'Color', [.95 .95 .95],...
%                     'GridColor',[1 1 1],...
%                     'GridAlpha',1);
                add_axes_background(elem,[.85 .85 .85],[1 1 1], .5, .35)

            case false
                set(elem, 'Color', [.95 .95 .95])
        end
    end
    
    if ~isempty(elem.Children)
%         children = findobj(elem,'-not','Type','Line');
        for ii = 1:numel(elem.Children)
            if ~any(strcmpi(elem.Children(ii).Type,{'line','patch'}))
                fix_axes(elem.Children(ii),font_size,grid_tf)
            end
        end
    end
    
    if ~strcmpi(elem.Type,'text')
        try
    %         fprintf("rep \n")
            set(elem,'FontSize',font_size); 
        end
    end
end