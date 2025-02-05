function result = getBooleanInput(question)
    while true
        response = input([question ' (y/n): '], 's'); % Get input as a string
        if strcmpi(response, 'y') % Accept 'Y' (case-insensitive)
            result = true;
            break;
        elseif strcmpi(response, 'n') % Accept 'N' (case-insensitive)
            result = false;
            break;
        else
            disp('Invalid input. Please enter "y" or "n".'); % Ask again for valid input
        end
    end
end