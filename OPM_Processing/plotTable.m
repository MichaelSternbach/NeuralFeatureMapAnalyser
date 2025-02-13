function fig = plotTable(tbl)
    % Function to plot a MATLAB table based on its content and return the figure handle
    fig = figure;
    hold on;
    
    % Check number of columns
    [numRows, numCols] = size(tbl);
    
    % Convert categorical variables to strings if needed
    for i = 1:numCols
        if iscategorical(tbl{:, i})
            tbl{:, i} = string(tbl{:, i});
        end
    end
    
    % Identify numeric and categorical columns
    numericCols = varfun(@isnumeric, tbl, 'OutputFormat', 'uniform');
    catCols = varfun(@iscategorical, tbl, 'OutputFormat', 'uniform');
    
    % If all columns are numeric, plot them
    if all(numericCols)
        plot(tbl.Variables);
        legend(tbl.Properties.VariableNames, 'Interpreter', 'none');
        xlabel('Row Index');
        ylabel('Value');
        title('Line Plot of Numeric Data');
    
    % If one categorical and one numeric column, use boxplot or bar
    elseif sum(catCols) == 1 && sum(numericCols) == 1
        catVar = tbl{:, catCols};
        numVar = tbl{:, numericCols};
        categories = unique(catVar);
        means = zeros(size(categories));
        
        for i = 1:length(categories)
            means(i) = mean(numVar(catVar == categories(i)));
        end
        
        bar(categorical(categories), means);
        xlabel(tbl.Properties.VariableNames{catCols});
        ylabel(tbl.Properties.VariableNames{numericCols});
        title('Bar Plot of Categorical vs Numeric Data');
    
    % If categorical and numerical mix, use scatter
    elseif sum(catCols) > 0 && sum(numericCols) > 0
        numVars = tbl{:, numericCols};
        catVars = tbl{:, catCols};
        
        scatter(categorical(catVars), numVars, 'filled');
        xlabel(tbl.Properties.VariableNames{catCols});
        ylabel(tbl.Properties.VariableNames{numericCols});
        title('Scatter Plot of Categorical vs Numeric Data');
    
    else
        error('Unsupported table format for plotting.');
    end
    
    hold off;
end
