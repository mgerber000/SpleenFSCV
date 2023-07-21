function tf = contains_(str1, str2, varargin)

if verLessThan('matlab', '9.1')
    p = inputParser;
    p.CaseSensitive = false;
    addRequired(p, 'str1', @isstr);
    addRequired(p, 'str2', @isstr);
    addParameter(p, 'IgnoreCase', false, @islogical);    
    parse(p, str1, str2, varargin{:});
    ignoreCase = p.Results.IgnoreCase;
    
    if ignoreCase
        tf = ~isempty(strfind(lower(p.Results.str1), lower(p.Results.str2)));  %#ok
    else
        tf = ~isempty(strfind(p.Results.str1, p.Results.str2));  
    end
else
    tf = contains(lower(str1),lower(str2));
    
%     tf = builtin('contains', str1, str2, varargin{:});
end
