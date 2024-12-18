function result = flattenCellArray(cellArray)
    % FLATTENCELLARRAY Converts a 2D cell array of vectors into a single vector
    %
    % Input:
    %   cellArray - A 2D cell array where each cell contains a vector of numbers
    %
    % Output:
    %   result - A single row vector containing all the values from the cell array

    % Ensure the input is a cell array
    if ~iscell(cellArray)
        error('Input must be a cell array.');
    end

    % Flatten the cell array to a single column
    flattened = cellArray(:);

    % Ensure all elements are row vectors
    flattened = cellfun(@(x) x(:)', flattened, 'UniformOutput', false);

    % Concatenate all values into a single row vector
    result = [flattened{:}];
end