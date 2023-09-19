function out = etaprint(etas)

p = [.01 .06 .14];
q = [mean([0,p(1)]), mean([p(1),p(2)]), mean([p(2),p(3)]), mean([p(3),1])];
r = [q(1), mean([p(1),q(2)]), mean([q(2),p(2)]), mean([p(2),q(3)]), mean([q(3),p(3)]), mean([p(3),1])];

nEtas = length(etas);
out = cell(size(etas));
for iEta = 1:nEtas
    eta = etas(iEta);
    if 0 <= eta && eta < r(1)
        out{iEta} = 'very small';
    elseif r(1) <= eta && eta < r(2)
        out{iEta} = 'small';
    elseif r(2) <= eta && eta < r(3)
        out{iEta} = 'small to medium';
    elseif r(3) <= eta && eta < r(4)
        out{iEta} = 'medium';
    elseif r(4) <= eta && eta < r(5)
        out{iEta} = 'medium to large';
    elseif r(5) <= eta && eta < r(6)
        out{iEta} = 'large';
    elseif r(6) <= eta && eta <= 1
        out{iEta} = 'very large';
    else
        if eta < 0
            warning('Eta must be between 0 and 1, and is %f',eta);
            out{iEta} = 'too small';
        else
            warning('Eta must be between 0 and 1, and is %f',eta);
            out{iEta} = 'too large';
        end
    end
end

if nEtas == 1
    out = out{1};
end

end