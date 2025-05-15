function [hasDuplicates, duplicatePairs] = checkDuplicatesTol(arr, tol)
%CHECKDUPLICATESTOL Checks for approximate duplicates within a tolerance
%
%   [hasDuplicates, duplicatePairs] = checkDuplicatesTol(arr, tol)
%
%   Inputs:
%       arr - A numerical array (vector or matrix)
%       tol - Tolerance for comparing floating-point numbers
%
%   Outputs:
%       hasDuplicates - Logical true if near-duplicates exist
%       duplicatePairs - Nx2 array of index pairs with near-duplicate values

    if nargin < 2
        tol = 1e-12;  % Default tolerance
    end

    arrVec = arr(:);
    n = length(arrVec);
    duplicatePairs = [];

    % Brute-force comparison with tolerance
    for i = 1:n-1
        for j = i+1:n
            if abs(arrVec(i) - arrVec(j)) < tol
                duplicatePairs = [duplicatePairs; i, j]; %#ok<AGROW>
            end
        end
    end

    hasDuplicates = ~isempty(duplicatePairs);
end
