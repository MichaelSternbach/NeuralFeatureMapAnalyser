%function read_Wallaby_data()
    %{

    %}
    %% Parameters
    addpath '/home/michael/Cloud//git/vone/MatlabCode/PatchAnalysis'
    folder ='~/Cloud/PhD/data/data share/Wallaby data/Maps/WallabyH/';
    DataFile = 'WallabyH_Image.mat';

    %load([folder,DataFile])

    % how many frames are per stimulus and trial
    frames_per_stim = 50;
    ImageData = dimg(2:9,:);

    % data frames to use
    frames_to_use = {11:16; ... % 1 sec before stim onset (blank)
        26:30; 31:35;... % 2 sec before stim offset %%%after stim onset
        36:40; 41:45}; % 2 sec after stim offset

    % %% Get stimulation protocol
    % [filename, pathname] = uigetfile('*.mat','Pick GratingIndex');
    % protocol = load([pathname,filename]);
    % 
    % % get list of presented conditions
    % stim_conditions = mod(protocol.Z + 90,360);
    % stim_conditions_unique = unique(stim_conditions);
    % frames_per_block = frames_per_stim * length(stim_conditions_unique);
    % 
    % %% Get files names
    % [filename, pathname] = uigetfile('*.qcamraw',['Pick qcam files for ',filename],'MultiSelect', 'on');
    % if ~iscell(filename)
    %     filename = {filename};
    % end

    %% Extract data from files
    stim_conditions_unique = (0:getStimulusNumber(ImageData)-1)/getStimulusNumber(ImageData)*180;
    num_blocks = getTrialNumber(ImageData);
    [pix_x, pix_y] = getMapSize(ImageData);

    % fill data array
    data = zeros([pix_x, pix_y,length(stim_conditions_unique),num_blocks,length(frames_to_use)]);
    %zeros(pix_y, pix_x, length(stim_conditions_unique), num_blocks, length(frames_to_use));

    for stim_ii=1:length(stim_conditions_unique)
        % find stim index where current condition is shown
        %stim_idx = find(stim_conditions == stim_conditions_unique(stim_ii));
        for block_ii=1:num_blocks
            % start index for frames to read
            %frame_idx = (stim_idx(block_ii)-1) * frames_per_stim;
            for frame_ii = 1:length(frames_to_use)
                % list of frames, remove the ones that go above num_frames
%                 frame_list = frame_idx + frames_to_use{frame_ii};
%                 frame_list(frame_list>num_frames) = [];
%                 assert(numel(frame_list)>4, 'Not enough images to make an averaged frame.')
                % read frames
                trial=ReturnTrial(ImageData,block_ii,stim_ii);
                % add average
                data(:,:,stim_ii,block_ii,frame_ii) = mean(trial(:,:,frames_to_use{frame_ii}),3);
            end
        end
    end



    %% Pre-process data array

    % extract blank (first frame)
    blanks = data(:,:,:,:,1);
    data(:,:,:,:,1) = [];

    %{
        % blank normalization  -(frame - first)./first
        data = - bsxfun(@rdivide,data,blanks) + 1;
    %}

    % block averaged blank normalization
    data = - bsxfun(@rdivide,data,mean(blanks,3)) + 1;

    % combine frames and blocks to have more blocks
    blanks = reshape(blanks,[pix_y, pix_x, 1, length(stim_conditions_unique)*num_blocks]);
    data = reshape(data,[pix_y, pix_x, length(stim_conditions_unique), num_blocks * size(data,5)]);

    % add shuffled blanks    
    [~,idx] = sort(rand(size(blanks,4),1));
    data = cat(3,data,blanks(:,:,1,idx(1:size(data,4))));
    stim_list = [stim_conditions_unique NaN];  

% %     % permute x and y coordinates
% %     data = permute(data,[2 1 3 4]);    

%     %% Save results
% 
%     save([folder,'WallabyH.mat'],'data','stim_list');

%%Plot Preprocessed data

plot_map(makeMap(data,stim_list))

%end

%%  HELPER FUNCTIONS TO READ HEADER AND FRAMES

function map =makeMap(data,stimuli_cond)
    map = zeros(size(data,1:2));
    for ii_condition = 1:size(stimuli_cond)
        if isnan(stimuli_cond(ii_condition))
            continue
        end
        tmp_img = mean(data(:,:,ii_condition,:),4);
        map = map + exp(1i*2*pi/180*real(stimuli_cond(ii_condition)))*tmp_img;
    end
end

function trial=ReturnTrial(ImageData,TrialIndex,StimulusIndex)
    trial = ImageData{StimulusIndex,TrialIndex};
end
function [px, py] = getMapSize(ImageData)
    MapSize = size(ImageData{1,1},1:2);
    px=MapSize(1);
    py=MapSize(2);
end

function TrialNumber = getTrialNumber(ImageData)
    TrialNumber = size(ImageData,2);
end

function FrameNumber = getFrameNumber(ImageData)
    FrameNumber = size(ImageData{1,1},3);
end

function StimulusNumber = getStimulusNumber(ImageData)
    StimulusNumber = size(ImageData,1);
end

function h = read_qcam_header(fileName)
%{
    Function to read the header of qcam raw data type
    
    Input:
        fileName = name of file to open with FULL path
    Output:
        h = struct with header details
%}

% open file
fid = fopen(fileName, 'r');
if fid<0
    error('Error opening file')
end

% add general info to header
h = [];
h.fileName = fileName;
fseek(fid, 0, 'eof');
h.fileSize = ftell(fid);

% get rest of header infomation from the file
frewind(fid);
gotThreeData = 0;
while gotThreeData < 3
    
    % read line and find info type
    tline = fgets(fid);
    [left, rem] = strtok( tline, ':');
    
    if strcmp(left, 'Fixed-Header-Size')
        right = strtok(rem, ':');
        h.headerSize = str2double(strtok(right));
        gotThreeData = gotThreeData +1;
        
    elseif strcmp(left, 'Frame-Size')
        right = strtok(rem, ':');
        h.frameSize = str2double(strtok(right));
        gotThreeData = gotThreeData +1;
        
    elseif strcmp(left, 'ROI')
        right = strtok(rem, ':');
        [~, right]=strtok(right, ',');
        [~, right]=strtok(right, ',');
        [width, right]=strtok(right, ',');
        height = strtok(right, ',');
        h.width = str2double(width);
        h.height = str2double(height);
        gotThreeData = gotThreeData +1;
        
    else % empty spaces
        continue;
    end
end
h.numberOfFrames = (h.fileSize - h.headerSize)/h.frameSize;

% close file
fclose(fid);

end

function stack = read_qcam_frames(h,frameNumber)
%{
    Function to read the frames of a qcam raw data
    
    Input: 
        h = struct with header details
        frameNumber = vector with frame IDs to extract
    Output:
        stack = 3D stack of extracted frames
%}

% check that number of frames matches requested frames
if max(frameNumber) > h.numberOfFrames
    error(['Frame number (' num2str(max(frameNumber)) ') exceeded maximum (' num2str(h.numberOfFrames) ').']);
end

% pre-allocate array
stack = zeros(h.width, h.height, length(frameNumber),'uint16');

% fill
for ii = 1:length(frameNumber)
    fseek(h.fid, h.headerSize + h.frameSize*(frameNumber(ii)-1), 'bof');
    frameData = fread(h.fid, [h.width, h.height], 'uint16');
    stack(:,:,ii) = frameData;
end

end