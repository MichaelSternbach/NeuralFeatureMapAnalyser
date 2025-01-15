function numList = checkFormatNumList(inputStr)
    % checkFormatNumList Converts a delimited string into a numerical array.
    %
    %   numList = checkFormatNumList(inputStr)
    %   Converts a string formatted as '1|2.3|4|5' into a numerical array [1 2.3 4 5].
    %
    %   Input:
    %       inputStr - A string with numbers separated by '|' characters.
    %
    %   Output:
    %       numList - A numerical array corresponding to the input string.
    %
    %   Example:
    %       checkFormatNumList('1|2.3|4|5')
    %       returns [1 2.3 4 5]

    % Check if inputStr is a valid string
    if ~ischar(inputStr) && ~isstring(inputStr)
        error('Input must be a string.');
    end

    % Replace '|' with spaces to prepare for conversion
    formattedStr = strrep(inputStr, '|', ' ');

    % Convert the formatted string to a numerical array
    numList = str2num(formattedStr); %#ok<ST2NM> Used for concise conversion

    % Ensure the conversion succeeded
    if isempty(numList) && ~isempty(inputStr)
        error('Input string could not be converted to a numerical array.');
    end
end
