function str = capitalize(str, option)

if nargin < 2
    option = 'first';
end

str = char(str);
switch option
    case 'all'
        idx = regexp([' ' str],'(?<=\s+)\S','start')-1;
    otherwise
        idx = 1;
end
if ~isempty(str)
    str(idx) = upper(str(idx));
end

end