function structUI()
    % Create a UI figure
    fig = uifigure('Name', 'Struct Variable Editor', 'Position', [100, 100, 500, 500]);

    % Create a blank struct
    myStruct = struct();

    % Define the initial variable names and default values
    variableNames = {'Variable1', 'Variable2', 'Variable3'};
    defaultValues = {0, 'Text', 3.14};

    % Store the struct and dynamic fields count in the figure's app data
    setappdata(fig, 'myStruct', myStruct);
    setappdata(fig, 'dynamicFieldCount', 0);

    % Create a panel for the existing fields
    existingFieldsPanel = uipanel(fig, ...
        'Title', 'Edit Existing Variables', ...
        'Position', [20, 200, 460, 270]);

    % Create UI components for each variable
    for i = 1:length(variableNames)
        uilabel(existingFieldsPanel, ...
            'Text', variableNames{i}, ...
            'Position', [10, 230 - 40 * i, 100, 22]);

        field = uieditfield(existingFieldsPanel, 'text', ...
            'Value', num2str(defaultValues{i}), ...
            'Position', [120, 230 - 40 * i, 200, 22], ...
            'Tag', variableNames{i});
        
        % Set the value change callback to update the struct
        field.ValueChangedFcn = @(src, event) updateStructValue(src, fig);
    end

    % Create a panel for dynamic fields
    dynamicFieldsPanel = uipanel(fig, ...
        'Title', 'Add and Edit Dynamic Variables', ...
        'Position', [20, 60, 460, 130]);

    % Add UI for the first dynamic field
    addDynamicField(dynamicFieldsPanel, fig);
    
    % Create a save button
    uibutton(fig, 'Text', 'Save', ...
        'Position', [200, 20, 100, 30], ...
        'ButtonPushedFcn', @(btn, event) saveStruct(fig, variableNames));
end

function addDynamicField(parentPanel, fig)
    % Retrieve the current dynamic field count
    fieldCount = getappdata(fig, 'dynamicFieldCount');
    fieldCount = fieldCount + 1;

    % Update the dynamic field count
    setappdata(fig, 'dynamicFieldCount', fieldCount);

    % Y-position for new dynamic fields
    yPos = 90 - 40 * (fieldCount - 1);

    % Create input for new key
    keyField = uieditfield(parentPanel, 'text', ...
        'Position', [10, yPos, 150, 22], ...
        'Placeholder', 'Enter key', ...
        'Tag', sprintf('DynamicKey%d', fieldCount));

    % Create input for new value
    valueField = uieditfield(parentPanel, 'text', ...
        'Position', [170, yPos, 150, 22], ...
        'Placeholder', 'Enter value', ...
        'Tag', sprintf('DynamicValue%d', fieldCount));

    % Set the value change callback for key and value
    keyField.ValueChangedFcn = @(src, event) updateStructValue(src, fig);
    valueField.ValueChangedFcn = @(src, event) updateStructValue(src, fig);
end

function updateStructValue(field, fig)
    % Retrieve the struct from app data
    myStruct = getappdata(fig, 'myStruct');

    % Get the key and value from the field
    key = field.Tag;
    newValue = field.Value;

    % Update the struct with the new value
    myStruct.(key) = newValue;

    % Store the updated struct in app data
    setappdata(fig, 'myStruct', myStruct);
end

function saveStruct(fig, variableNames)
    % Get the struct stored in the app data
    myStruct = getappdata(fig, 'myStruct');

    % Get updated values from the UI for existing variables
    for i = 1:length(variableNames)
        field = findobj(fig, 'Tag', variableNames{i});
        value = field.Value;

        % Convert numeric fields back to numbers if possible
        numValue = str2double(value);
        if ~isnan(numValue)
            value = numValue;
        end

        % Assign the value to the struct
        myStruct.(variableNames{i}) = value;
    end

    % Display the updated struct
    disp('Updated Struct:');
    disp(myStruct);

    % Save updated struct back to app data
    setappdata(fig, 'myStruct', myStruct);

    % Notify the user of successful saving
    uialert(fig, 'Struct updated successfully!', 'Success');
end
