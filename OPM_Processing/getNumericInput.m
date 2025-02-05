function result = getNumericInput(question)
    while true
        response = input(question); % Get input
        if isnumeric(response) % Check if the input is numeric
            result = response; % Assign the numeric input to the result
            break; % Exit the loop if valid input
        else
            disp('Invalid input. Please enter numeric values (e.g., a number or an array).');
        end
    end
end
