function varargout = blkImportDlg(varargin)
% BLKIMPORTDLG M-file for blkImportDlg.fig
%      BLKIMPORTDLG, by itself, creates a new BLKIMPORTDLG or raises the existing
%      singleton*.
%
%      H = BLKIMPORTDLG returns the handle to a new BLKIMPORTDLG or the handle to
%      the existing singleton*.
%
%      BLKIMPORTDLG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BLKIMPORTDLG.M with the given input arguments.
%
%      BLKIMPORTDLG('Property','Value',...) creates a new BLKIMPORTDLG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before blkImportDlg_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to blkImportDlg_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help blkImportDlg

% Last Modified by GUIDE v2.5 11-Jan-2010 11:13:56

% $Id: blkImportDlg.m,v 1.1 2008/04/30 04:56:07 shaunc Exp $

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @blkImportDlg_OpeningFcn, ...
                   'gui_OutputFcn',  @blkImportDlg_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before blkImportDlg is made visible.
function blkImportDlg_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to blkImportDlg (see VARARGIN)

% set path to search
blkFilePath = pwd;
if nargin >= 4,
  blkFilePath = varargin{2};
end;
handles.blkFilePath = blkFilePath;

condIdFilePath = blkFilePath;
handles.condIdFilePath = condIdFilePath;

% populate the pathMenu
tmp = get(handles.pathMenu, 'String');
tmp = addToList(tmp, blkFilePath);
set(handles.pathMenu, 'String', tmp);
set(handles.pathMenu, 'Value', 1);

% default to single condition BLK files
handles.singleCondBlkFiles = 1;

% populate the condIdFilePathMenu
tmp = get(handles.condIdFilePathMenu, 'String');
tmp = addToList(tmp, condIdFilePath);
set(handles.condIdFilePathMenu, 'String', tmp);
set(handles.condIdFilePathMenu, 'Value', 1);

% Choose default command line output for blkImportDlg
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% initialise uicontrols
updateUIControls(handles.blkImportDlg);

% UIWAIT makes blkImportDlg wait for user response (see UIRESUME)
% uiwait(handles.blkImportDlg);


% --- Outputs from this function are returned to the command line.
function varargout = blkImportDlg_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in quitButton.
function quitButton_Callback(hObject, eventdata, handles)
% hObject    handle to quitButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% close the figure
close(handles.blkImportDlg);


% --- Executes on button press in importButton.
function importButton_Callback(hObject, eventdata, handles)
% hObject    handle to importButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% 1. get the selected file names from the listbox
blkFileIdx = get(handles.blkFileListbox, 'Value');

% 2. create the full paths to the files
blkFilePath = handles.blkFilePath;
%blkFileInfo = handles.blkFileInfo;

names = cell(size(blkFileIdx));
%[names{:}] = deal(blkFileInfo(idx).name);
names = get(handles.blkFileListbox,'String');
names = cellfun(@(x) fullfile(blkFilePath, x), names(blkFileIdx), 'UniformOutput', false);

% 3. open the files (i.e., get the fid's)
blkFileIds = cellfun(@(x) fopen(x,'r'), names);

fprintf(1, 'Importing %i .BLK files... ', length(blkFileIds));

% 4 import the imaging data (i.e., call blkImport( ))

singleCondBlkFiles = handles.singleCondBlkFiles;
if ~singleCondBlkFiles,
  data = blkImport(blkFileIds);
else,
  % load condIds, group by condition, load data by calling blkImport()
  condIdFileIdx = get(handles.condIdFileListbox, 'Value');
  
  condIdFilePath = handles.condIdFilePath;

  names = cell(size(condIdFileIdx));
%   [names{:}] = deal(condIdFileInfo(idx).name);
  names = get(handles.condIdFileListbox,'String');
  names = cellfun(@(x) fullfile(condIdFilePath, x), names(condIdFileIdx), 'UniformOutput', false);
   
  condIdFileIds = cellfun(@(x) fopen(x,'r'), names);
  [conds, idx] = expoStimImport(condIdFileIds(1)); % FIXME: only one!
%   [conds, idx] = alexStimImport(condIdFileIds(1)); % Alex's stimbox

  jnk = arrayfun(@(x) fclose(x), condIdFileIds);
  
  % uncomment for 'offline' testing only
%  conds = [0:30:330, NaN]; % [contrast; direction]
%  idx = randperm(13)';
  
  % repeat presentation order if we have more files than condition indicies
  idx = repmat(idx, ceil(length(blkFileIds)/length(idx)), 1);
  idx = idx(1:length(blkFileIds));
    
  for i = 1:length(conds),
    j = find(idx == i);
    data(i,:) = blkImport(blkFileIds(j));
%    data{i,:} = j(:)';
  end
end

fprintf(1, 'Ok! (%i conditions, %i trials).\n', size(data));

% 5. copy the imaging data to Matlab's 'base' workspace
assignin('base','images',data);

% 6. close the files
jnk = arrayfun(@(x) fclose(x), blkFileIds);


% --- Executes on button press in reloadButton.
function reloadButton_Callback(hObject, eventdata, handles)
% hObject    handle to reloadButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

updateUIControls(handles.blkImportDlg);


% --- Executes on selection change in blkFileListbox.
function blkFileListbox_Callback(hObject, eventdata, handles)
% hObject    handle to blkFileListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns blkFileListbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from blkFileListbox


% --- Executes during object creation, after setting all properties.
function blkFileListbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to blkFileListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pathMenu.
function pathMenu_Callback(hObject, eventdata, handles)
% hObject    handle to pathMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns pathMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pathMenu

tmp = get(handles.pathMenu, 'String');
idx = get(handles.pathMenu, 'Value');

blkFilePath = tmp{idx};
handles.blkFilePath = blkFilePath;
guidata(hObject, handles);

updateUIControls(handles.blkImportDlg);


% --- Executes during object creation, after setting all properties.
function pathMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pathMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pathBrowseButton.
function pathBrowseButton_Callback(hObject, eventdata, handles)
% hObject    handle to pathBrowseButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

blkFilePath = handles.blkFilePath;

tmp = 'Select folder to search for .BLK files';
blkFilePath = uigetdir(blkFilePath, tmp);

if blkFilePath == 0,
  return
end

tmp = get(handles.pathMenu, 'String');
tmp = addToList(tmp, blkFilePath);

set(handles.pathMenu, 'String', tmp);
set(handles.pathMenu, 'Value', length(tmp));

handles.blkFilePath = blkFilePath;
guidata(hObject,handles);

updateUIControls(handles.blkImportDlg);


% --- Executes on button press in singleCondBlkFileCheckbox.
function singleCondBlkFileCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to singleCondBlkFileCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of singleCondBlkFileCheckbox

% get application data - hObject is the application window object
handles = guidata(hObject);

singleCondBlkFiles = get(hObject,'Value');

handles.singleCondBlkFiles = singleCondBlkFiles;
guidata(hObject,handles);

updateUIControls(handles.blkImportDlg);


% --- Executes on selection change in condIdFilePathMenu.
function condIdFilePathMenu_Callback(hObject, eventdata, handles)
% hObject    handle to condIdFilePathMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns condIdFilePathMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from condIdFilePathMenu

tmp = get(handles.condIdFilePathMenu, 'String');
idx = get(handles.condIdFilePathMenu, 'Value');

condIdFilePath = tmp{idx};
handles.condIdFilePath = condIdFilePath;
guidata(hObject, handles);

updateUIControls(handles.blkImportDlg);


% --- Executes during object creation, after setting all properties.
function condIdFilePathMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to condIdFilePathMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in condIdFilePathBrowseButton.
function condIdFilePathBrowseButton_Callback(hObject, eventdata, handles)
% hObject    handle to condIdFilePathBrowseButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

condIdFilePath = handles.condIdFilePath;

tmp = 'Select folder to search for .TXT files';
condIdFilePath = uigetdir(condIdFilePath, tmp);

if condIdFilePath == 0,
  return
end

tmp = get(handles.condIdFilePathMenu, 'String');
tmp = addToList(tmp, condIdFilePath);

set(handles.condIdFilePathMenu, 'String', tmp);
set(handles.condIdFilePathMenu, 'Value', length(tmp));

handles.condIdFilePath = condIdFilePath;
guidata(hObject,handles);

updateUIControls(handles.blkImportDlg);


% --- Executes on selection change in condIdFileListbox.
function condIdFileListbox_Callback(hObject, eventdata, handles)
% hObject    handle to condIdFileListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns condIdFileListbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from condIdFileListbox


% --- Executes during object creation, after setting all properties.
function condIdFileListbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to condIdFileListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%
%%% private functions...
%%%

function updateUIControls(hObject)

% get application data - hObject is the application window object
handles = guidata(hObject);

blkFilePath = handles.blkFilePath;

singleCondBlkFiles = handles.singleCondBlkFiles;
condIdFilePath = handles.condIdFilePath;

blkFileInfo = dir(fullfile(blkFilePath, '*.BLK'));

names = cell(size(blkFileInfo));
[names{:}] = deal(blkFileInfo.name);

%
pat = 'E(?<exp>\d+)B(?<blk>\d+)';
t = regexp(names, pat, 'names');

expIds = cellfun(@(x) str2num(x.exp), t);
blkIds = cellfun(@(x) str2num(x.blk), t);

[y,i] = sortrows([expIds, blkIds]);

names = names(i);
%

set(handles.blkFileListbox, 'String', names);

%handles.blkFileInfo = blkFileInfo; % FIXME: resort this!?

% enable/disable single condition block file controls
set(handles.singleCondBlkFileCheckbox,'Value',singleCondBlkFiles);
if singleCondBlkFiles,
  set(handles.condIdFilePathEntryBoxLabel,'Enable','on');
  set(handles.condIdFilePathMenu,'Enable','on');
  set(handles.condIdFilePathBrowseButton,'Enable','on');
  set(handles.condIdFileListbox,'Enable','on');
else,
  set(handles.condIdFilePathEntryBoxLabel,'Enable','off');
  set(handles.condIdFilePathMenu,'Enable','off');
  set(handles.condIdFilePathBrowseButton,'Enable','off');
  set(handles.condIdFileListbox,'Enable','off');
end

condIdFileInfo = dir(fullfile(condIdFilePath, '*.txt'));

names = cell(size(condIdFileInfo));
[names{:}] = deal(condIdFileInfo.name);

set(handles.condIdFileListbox, 'String', names);

%handles.condIdFileInfo = condIdFileInfo;
guidata(hObject,handles);


function tmp = addToList(list, item, maxLen),

if nargin < 3,
  maxLen = 5;
end

if isempty(list),
  list = {item};
else,
  list{end+1} = item;
end

if length(list) > maxLen,
  list = list(2:maxLen+1);
end

tmp = list;





