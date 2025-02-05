function writeStructToJSON(structData, filename)
    % Converts a MATLAB structure to a JSON file
    % Inputs:
    %   structData: MATLAB structure to be converted
    %   filename: Name of the JSON file to save the data
    
    % Convert structure to JSON string
    jsonString = jsonencode(structData);
    
    % Write JSON string to the file
    fid = fopen(filename, 'w');
    if fid == -1
        error('Cannot open file %s for writing.', filename);
    end
    fprintf(fid, '%s', jsonString);
    fclose(fid);
end