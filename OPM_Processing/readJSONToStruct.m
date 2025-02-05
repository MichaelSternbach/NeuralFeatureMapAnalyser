function structData = readJSONToStruct(filename)
    % Reads a JSON file and converts it to a MATLAB structure
    % Inputs:
    %   filename: Name of the JSON file to read
    % Output:
    %   structData: MATLAB structure containing the parsed JSON data
    
    % Read JSON string from the file
    fid = fopen(filename, 'r');
    if fid == -1
        error('Cannot open file %s for reading.', filename);
    end
    raw = fread(fid, inf, 'char');
    fclose(fid);
    jsonString = char(raw');
    
    % Convert JSON string to MATLAB structure
    structData = jsondecode(jsonString);
end