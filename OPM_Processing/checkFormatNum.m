function var = checkFormatNum(var)
    if isstring(var)
        var = str2num(var);
    elseif ischar(var)
        var = str2num(convertCharsToStrings(var));
    end

end