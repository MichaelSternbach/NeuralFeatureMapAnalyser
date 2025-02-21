function A = removeNanColumns(A)
    A(:, any(isnan(A), 1)) = [];
end