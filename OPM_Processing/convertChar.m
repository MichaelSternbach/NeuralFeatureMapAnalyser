function result = convertChar(inputStr)
    seperator = '|';
    % Try to convert the input to a numeric value
    numericValue = str2double(inputStr);
    if ~isnan(numericValue)
        % If the input is numeric, return it as a number
        result = numericValue;
    elseif isstruct(inputStr)
        result = inputStr;
    elseif isnumeric(inputStr)
        result = inputStr;
    elseif contains(inputStr, ':') && contains(inputStr, seperator)
        % If the input contains ':' and ',' assume it's a struct format
        result = charToStruct(inputStr,seperator);
    else
        disp(inputStr)
        error('Input format not recognized. It should be either numeric or a struct-like string.');
    end
end

function result = charToStruct(inputStr,seperator)
    % Helper function to convert struct-like strings to MATLAB structures
    pairs = strsplit(inputStr, seperator);
    result = struct(); % Initialize an empty structure

    for i = 1:length(pairs)
        % Split each pair by the colon
        keyValue = strsplit(pairs{i}, ':');
        if length(keyValue) == 2
            key = strtrim(keyValue{1}); % Remove any whitespace around the key
            value = str2double(keyValue{2}); % Convert the value to a number if possible
            if isnan(value) % If not a number, keep as a string
                value = strtrim(keyValue{2});
            end
            % Assign the value to the corresponding field in the structure
            result.(key) = value;
        else
            error('Invalid format for key-value pair: %s', pairs{i});
        end
    end
end
