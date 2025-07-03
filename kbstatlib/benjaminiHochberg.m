function [adjusted_pvals, rejected] = benjaminiHochberg(pvals, fdr)
%% Adjust p-values according to the Benjamini-Hochberg method
% This method aims to control the false detection rate (FDR) of a set of hypothesis tests to a given
% level. The FDR is the expected proportion of false discoveries among all discoveries, where a
% "discovery" is a rejection of the null hypothesis. The difference to other methods like Bonferroni
% or Holm-Bonferroni is that the latter aim at controlling the Family-Wise Error Rate (FWER), which
% is the probability of making at least one false discovery.
%
% SYNTAX
%   [adjusted_pvals, rejected] = benjaminiHochberg(pvals, fdr)
%
% INPUT
%   pvals   - Array of p-values (vector)
%   fdr     - Desired FDR level (e.g., 0.05)
%
% OUTPUT
%   adjusted_pvals  - Benjamini-Hochberg adjusted p-values
%   rejected        - Logical array indicating which hypotheses are rejected

if nargin < 2
    fdr = 0.05;
end

origSize = size(pvals);

pvals = pvals(:); % Ensure column vector
m = length(pvals); % Number of tests
[p_sorted, sort_idx] = sort(pvals);
[~, unsort_idx] = sort(sort_idx); % To restore original order

% Compute BH thresholds
bh_thresh = (1:m)' / m * fdr;

% Find largest k such that p(k) <= (k/m)*q
below_thresh = p_sorted <= bh_thresh;
if any(below_thresh)
    k_max = find(below_thresh, 1, 'last');
    rejected_sorted = false(m,1);
    rejected_sorted(1:k_max) = true;
else
    rejected_sorted = false(m,1);
end

% Adjusted p-values (step-up procedure)
p_adj_sorted = min( cummin((m ./ (1:m)') .* p_sorted), 1 );

% Reorder and reshape to match original input
rejected = reshape(rejected_sorted(unsort_idx), origSize);
adjusted_pvals = reshape(p_adj_sorted(unsort_idx), origSize);

end
