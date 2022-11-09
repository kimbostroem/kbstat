function str = sigprint(p)
str = cell(size(p));
for iP = 1:length(p)
    if p(iP) <= 0.001
        str{iP} = '***';
    elseif p(iP) <= 0.01
        str{iP} = '**';
    elseif p(iP) <= 0.05
        str{iP} = '*';
    else
        str{iP} = '';
    end
end

if length(p)==1
    str = str{1};
end

end