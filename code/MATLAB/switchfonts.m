function switchfonts(state)
% LARGEFONTS ON   Turn on large fonts
% LARGEFONTS OFF  Turn off large fonts
% LARGEFONTS      Toggle large fonts
% LARGEFONTS FS   Turn on large fonts, with font size FS.  

% Get largefonts preferences
FS = getpref('LargeFonts');
if isempty(FS)  % Create preferences if necessary
    setpref('LargeFonts','State','Off'); 
end

% Parse input arguments
if nargin==0            % Toggle state
    state = lower(getpref('LargeFonts','State'));
    switch state
        case 'on'
            state = 'off';
        case 'off'
            state = 'on';
    end
elseif ~isstr(state)  % Specified font size.  Turn on.
    FontSize = state;
    state = 'on';       % Turn on
else                    % Go to specified state
    state= lower(state);
end

if ~exist('FontSize','var')
    FontSize = 30;      % Default size for large fonts
end

%%Toggle font sizes
switch state
    case 'on'
        % Big-
        com.mathworks.services.FontPrefs.setCodeFont(java.awt.Font('Consolas',java.awt.Font.PLAIN,FontSize))
        com.mathworks.services.FontPrefs.setTextFont(java.awt.Font('Calibri',java.awt.Font.PLAIN,FontSize))
        setpref('LargeFonts','State','On');
        set(0, 'DefaultUIControlFontSize',20);
    case 'off'
        % Small-
        com.mathworks.services.FontPrefs.setCodeFont(java.awt.Font('Consolas',java.awt.Font.PLAIN,16))
        com.mathworks.services.FontPrefs.setTextFont(java.awt.Font('Calibri',java.awt.Font.PLAIN,16))
        setpref('LargeFonts','State','Off');
        set(0, 'DefaultUIControlFontSize',10);
end;