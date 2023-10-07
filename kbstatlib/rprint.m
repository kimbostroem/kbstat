function out = rprint(rs)

p = [.2 .5 .8];
q = [mean([0,p(1)]), mean([p(1),p(2)]), mean([p(2),p(3)]), mean([p(3),1])];
ranges = [q(1), mean([p(1),q(2)]), mean([q(2),p(2)]), mean([p(2),q(3)]), mean([q(3),p(3)]), mean([p(3),1])];

nRs = length(rs);
out = cell(size(rs));
for iR = 1:nRs
    r = abs(rs(iR));
    if 0 <= r && r < ranges(1)
        out{iR} = 'very small';
    elseif ranges(1) <= r && r < ranges(2)
        out{iR} = 'small';
    elseif ranges(2) <= r && r < ranges(3)
        out{iR} = 'small to medium';
    elseif ranges(3) <= r && r < ranges(4)
        out{iR} = 'medium';
    elseif ranges(4) <= r && r < ranges(5)
        out{iR} = 'medium to large';
    elseif ranges(5) <= r && r < ranges(6)
        out{iR} = 'large';
    elseif ranges(6) <= r && r <= 1
        out{iR} = 'very large';
    else
        if r < 0
            warning('Cohen''s d must be between 0 and 1, and is %f', r);
            out{iR} = 'too small';
        else
            warning('Cohen''s d must be between 0 and 1, and is %f', r);
            out{iR} = 'too large';
        end
    end
end

if nRs == 1
    out = out{1};
end

end