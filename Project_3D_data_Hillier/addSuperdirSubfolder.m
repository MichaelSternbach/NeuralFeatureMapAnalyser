function addSuperdirSubfolder(subfolderName)
% ADDSUPERDIRSUBFOLDER Adds a subfolder (inside the superdirectory of this script) to the MATLAB path.
%
% Usage:
%   addSuperdirSubfolder('utils')
%
% This adds the folder <superdir>/utils to the MATLAB path, where <superdir>
% is one level above the folder containing this function or script.

    % Get the full path of the calling script or function
    callerPath = mfilename('fullpath');

    % Get the directory containing the current file
    currentDir = fileparts(callerPath);

    % Get the superdirectory (parent of the current folder)
    superDir = fileparts(currentDir);

    % Construct full path to the subfolder
    targetSubfolder = fullfile(superDir, subfolderName);

    % Check if the folder exists
    if isfolder(targetSubfolder)
        addpath(targetSubfolder);
        disp(['[addSuperdirSubfolder] Added to path: ', targetSubfolder]);
    else
        warning('[addSuperdirSubfolder] Folder does not exist: %s', targetSubfolder);
    end
end
