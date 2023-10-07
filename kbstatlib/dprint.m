function out = dprint(ds)

p = [.2 .5 .8];
q = [mean([0,p(1)]), mean([p(1),p(2)]), mean([p(2),p(3)]), mean([p(3),1])];
r = [q(1), mean([p(1),q(2)]), mean([q(2),p(2)]), mean([p(2),q(3)]), mean([q(3),p(3)]), mean([p(3),1])];

nDs = length(ds);
out = cell(size(ds));
for iD = 1:nDs
    d = ds(iD);
    if 0 <= d && d < r(1)
        out{iD} = 'very small';
    elseif r(1) <= d && d < r(2)
        out{iD} = 'small';
    elseif r(2) <= d && d < r(3)
        out{iD} = 'small to medium';
    elseif r(3) <= d && d < r(4)
        out{iD} = 'medium';
    elseif r(4) <= d && d < r(5)
        out{iD} = 'medium to large';
    elseif r(5) <= d && d < r(6)
        out{iD} = 'large';
    elseif r(6) <= d && d <= 1
        out{iD} = 'very large';
    else
        if d < 0
            warning('Cohen''s d must be between 0 and 1, and is %f', d);
            out{iD} = 'too small';
        else
            warning('Cohen''s d must be between 0 and 1, and is %f', d);
            out{iD} = 'too large';
        end
    end
end

if nDs == 1
    out = out{1};
end

end