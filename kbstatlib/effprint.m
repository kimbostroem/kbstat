function out = effprint(effects, effectType)
%% Interpret effect sizes of a certain types in terms of 'small', 'medium' 
%% and 'large' according to Cohen's rule of thumb.
% Source: Cohen, J. 1988. Statistical Power Analysis for the Behavioral
% Sciences. Vol. 2nd. Hillsdale, New Jersey: Lawrence Erlbaum Associates.
%
% SYNTAX
%   out = effprint(effects, effectType)
%
% INPUT
%   effects     (Nx1 or 1xN double) Effect sizes
%   effectType  (char)  Type of effect size. Possible values:
%                   'eta2' eta squared (for one-way ANOVA) or partial 
%                           eta squared (for multi-way ANOVA)
%                   'd'     Cohen's d
%                   'r'     Pearson's correlation coefficient r
%                   'rho'   Spearman's correlation coefficient rho
%                   'R2'    Coefficient of determination, the square of 
%                           the Pearson correlation r.
%                   'f2'    Amount of bias in an F-test (ANOVA or multiple 
%                           regression)
%                   'delta' Cliff's delta, effect size for ordinal data.
%                           Measure of how often the values in one
%                           distribution are larger than the values in a
%                           second distribution
%                   'omega' Effect size used for chi-squared test
%                   'g'     Hedges' g, based on a standardized difference.

switch effectType
    case 'eta2'
        p = [.01 .06 .14];
    case 'd'
        p = [.2 .5 .8];
    case 'r'
        p = [.1 .3 .5];
    case 'rho'
        p = [.1 .3 .5];
    case 'R2'
        p = [.1 .3 .5].^2;
    case 'f2'
        p = [.02 .15 .35];
    case 'delta'
        p = [.25 .75 1.25];
    case 'omega'
        p = [.1 .3 .5];
    case 'g'
        p = [.05 .15 .25];
end

q = [mean([0,p(1)]), mean([p(1),p(2)]), mean([p(2),p(3)]), mean([p(3),1])];
ranges = [q(1), mean([p(1),q(2)]), mean([q(2),p(2)]), mean([p(2),q(3)]), mean([q(3),p(3)]), mean([p(3),1])];

nEffects = length(effects);
out = cell(size(effects));
for iEffect = 1:nEffects
    effect = abs(effects(iEffect));
    if 0 <= effect && effect < ranges(1)
        out{iEffect} = 'very small';
    elseif ranges(1) <= effect && effect < ranges(2)
        out{iEffect} = 'small';
    elseif ranges(2) <= effect && effect < ranges(3)
        out{iEffect} = 'small to medium';
    elseif ranges(3) <= effect && effect < ranges(4)
        out{iEffect} = 'medium';
    elseif ranges(4) <= effect && effect < ranges(5)
        out{iEffect} = 'medium to large';
    elseif ranges(5) <= effect && effect < ranges(6)
        out{iEffect} = 'large';
    elseif ranges(6) <= effect && effect <= 1
        out{iEffect} = 'very large';
    else
        if effect < 0
            warning('Eta must be between 0 and 1, and is %f',effect);
            out{iEffect} = 'too small';
        else
            warning('Eta must be between 0 and 1, and is %f',effect);
            out{iEffect} = 'too large';
        end
    end
end

if nEffects == 1
    out = out{1};
end

end