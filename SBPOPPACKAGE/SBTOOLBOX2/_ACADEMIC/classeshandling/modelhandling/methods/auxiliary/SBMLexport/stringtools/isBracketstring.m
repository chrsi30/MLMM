function [isBracketstring] = isBracketstring(string)
%isBracketstring
% a simple test wether the given string contains squared brackets, curly
% brackets or parenthesis
%
% USAGE:
% ======
%
% [boolean] = isBracketString(string)
%
% string: test string
%
% boolean: true if string contains (, ), {, }, [ or ]
%          false otherwise

% Information:
% ============
% Author: Gunnar Drews, gunnar.drews@uni-rostock.de

    if strfind(string, '('),
        isBracketstring = 1;
        return;
    end
    if strfind(string, ')'),
        isBracketstring = 1;
        return;
    end
    if strfind(string, '['),
        isBracketstring = 1;
        return;
    end
    if strfind(string, ']'),
        isBracketstring = 1;
        return;
    end
    if strfind(string, '{'),
        isBracketstring = 1;
        return;
    end
    if strfind(string, '}'),
        isBracketstring = 1;
        return;
    end
    isBracketstring = 0;
return