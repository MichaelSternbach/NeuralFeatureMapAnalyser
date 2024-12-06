function fieldTable = countFieldOccurrencesWithTypes(structArray)
    % Function to count occurrences of each field name, store data types, and compute relative frequency.
    % INPUT:
    %   structArray - Cell array of structures to be analyzed.
    % OUTPUT:
    %   fieldTable - Table summarizing counts, data types, and relative frequency of fields.

    if ~iscell(structArray) || isempty(structArray) || ~all(cellfun(@isstruct, structArray))
        error('Input must be a non-empty cell array of structures.');
    end
    
    % Initialize a map to track field counts and data types.
    fieldDataSummary = struct();
    totalStructures = numel(structArray);
    
    % Loop over all structures in the input array
    for i = 1:totalStructures
        % Get field names for the current structure
        fields = fieldnames(structArray{i});
        
        % Update counts and data types for each field
        for j = 1:numel(fields)
            field = fields{j};
            % Determine data type of the current field value
            dataType = class(structArray{i}.(field));
            
            if isfield(fieldDataSummary, field)
                % Increment occurrence count
                fieldDataSummary.(field).count = fieldDataSummary.(field).count + 1;
                % Add data type if not already in the list
                if ~ismember(dataType, fieldDataSummary.(field).dataTypes)
                    fieldDataSummary.(field).dataTypes{end+1} = dataType; %#ok<AGROW>
                end
            else
                % Initialize count and data type list for new fields
                fieldDataSummary.(field).count = 1;
                fieldDataSummary.(field).dataTypes = {dataType};
            end
        end
    end
    
    % Prepare data for table output
    fieldNames = fieldnames(fieldDataSummary);
    counts = zeros(numel(fieldNames), 1);
    dataTypes = cell(numel(fieldNames), 1);
    relativeOccurrences = zeros(numel(fieldNames), 1);
    
    for i = 1:numel(fieldNames)
        counts(i) = fieldDataSummary.(fieldNames{i}).count;
        dataTypes{i} = strjoin(fieldDataSummary.(fieldNames{i}).dataTypes, ', ');
        relativeOccurrences(i) = (counts(i) / totalStructures) * 100;
    end
    
    % Create output table
    fieldTable = table(fieldNames, counts, relativeOccurrences, dataTypes, ...
        'VariableNames', {'FieldName', 'Count', 'RelativeOccurrencePercent', 'DataTypes'});
    
    % Sort the table by Count in descending order
    fieldTable = sortrows(fieldTable, 'Count', 'descend');
end
