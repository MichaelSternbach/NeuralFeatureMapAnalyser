function h = blkFileHeader(fid),
%BLKFILEHEADER Read the file header from a Vdaq .blk file.
%
%  H = BLKFILEHEADER(FID) reads and returns the file header
%  from the specified file handle.

% $Id: blkFileHeader.m,v 1.1 2008-04-30 04:56:07 shaunc Exp $

% read the file header

try,
  frewind(fid);
catch,
  
  warning(['blkFileHeader:' ferror(fid) 'Invalid file handle?' ]);
  h = [];
  return;
end;

% data integrity
h.fileSize = fread(fid, 1, 'uint32');
h.headerCheckSum = fread(fid, 1, 'uint32');
h.dataCheckSum = fread(fid, 1, 'uint32');

% common to all data files
h.headerLen = fread(fid, 1, 'uint32'); % bytes
h.versionId = fread(fid, 1, 'uint32');
h.fileType = fread(fid, 1, 'uint32');
h.fileSubtype = fread(fid, 1, 'uint32');
h.dataType = fread(fid, 1, 'uint32');
h.sizeOf = fread(fid, 1, 'uint32');
h.frameWidth = fread(fid, 1, 'uint32'); % pixels
h.frameHeight = fread(fid, 1, 'uint32'); % pixels
h.framesPerStim = fread(fid, 1, 'uint32');
h.numStim = fread(fid, 1, 'uint32');
h.initBinFactor = fread(fid, 2, 'uint32'); % x, y
h.binFactor = fread(fid, 2, 'uint32'); % x, y

%h.userName = fread(fid, 32, 'uchar');
%h.date = fread(fid, 16, 'uchar'); % 'mm/dd/?yy       '
h.userName = fscanf(fid, '%c', 32);
h.userName = h.userName(h.userName ~= 0);
h.date = fscanf(fid, '%c', 16); % 'mm/dd/?yy       '
h.date = h.date(h.date ~= 0);

h.roi = fread(fid, 4, 'uint32'); % x1, y1, x2, y2

% data and reference frames
h.stimOffset = fread(fid, 1, 'uint32');
h.stimSize = fread(fid, 1, 'uint32');
h.frameOffset = fread(fid, 1, 'uint32');
h.frameSize = fread(fid, 1, 'uint32'); % bytes

% note: the Imager 3001 has no reference frame, these fields will be zero
h.refOffset = fread(fid, 1, 'uint32');
h.refSize = fread(fid, 1, 'uint32');
h.refWidth = fread(fid, 1, 'uint32');
h.refHeight = fread(fid, 1, 'uint32');

% compression or summing
h.blocks = fread(fid, 16, 'uint16'); % max 256 blocks per experiment
h.frames = fread(fid, 16, 'uint16'); % max 256 frames per condition

% data analysis
h.lowClip = fread(fid, 1, 'int32');
h.highClip = fread(fid, 1, 'int32');
h.lowPass = fread(fid, 1, 'uint32');
h.highPass = fread(fid, 1, 'uint32');
%h.ops = fread(fid, 64, 'uchar');
h.ops = fscanf(fid, '%c', 64);
h.ops = h.ops(h.ops ~= 0);

% ora specific (not used by Vdaq)
h.ora.magnification = fread(fid, 1, 'int32');
h.ora.gain = fread(fid, 1, 'uint16');
h.ora.wavelength = fread(fid, 1, 'uint16');
h.ora.exposureTime = fread(fid, 1, 'uint32');
h.ora.numReps = fread(fid, 1, 'uint32');
h.ora.acquisitionDelay = fread(fid, 1, 'uint32');
h.ora.interStimInterval = fread(fid, 1, 'uint32');
%h.ora.creationDate = fread(fid, 16, 'uchar');
%h.ora.filename = fread(fid, 64, 'uchar');
h.ora.creationDate = fscanf(fid, '%c', 16);
h.ora.creationDate = h.ora.creationDate(h.ora.creationDate ~= 0);
h.ora.filename = strtrim(fscanf(fid, '%c', 64));
h.ora.filename = h.ora.filename(h.ora.filename ~= 0);
h.ora.reserved = fread(fid, 256, 'uchar');

% Vdaq specific
h.refFramePresent = fread(fid, 1, 'uint32'); % 0 or 1
%h.stimIds = fread(fid, 256, 'uchar');
h.stimIds = sscanf(fscanf(fid, '%c', 256), '%i')';
h.videoFramesPerDataFrame = fread(fid, 1, 'uint32');
h.numTrials = fread(fid, 1, 'uint32');
h.scaleFactor = fread(fid, 1, 'uint32');
h.meanAmpGain = fread(fid, 1, 'float32');
h.meanAmpDC = fread(fid, 1, 'float32');

h.baseline = fread(fid, 2, 'uint8'); % [first, last] frame number
h.activity = fread(fid, 2, 'uint8'); % [first, last] frame number

h.digitizerBits = fread(fid, 1, 'uint8');
h.systemId = fread(fid, 1, 'uint8');

h.dummy2 = fread(fid, 1, 'uint8'); % wah!?
h.dummy3 = fread(fid, 1, 'uint8');

h.superPixel = fread(fid, 4, 'uint32'); % [x1, y1, x2, y2]

h.frameDuration = fread(fid, 1, 'float32'); % milliseconds
h.validFrames = fread(fid, 1, 'uint32');

h.reserved = fread(fid, 224, 'uchar'); % reserved for Vdaq

% timeing information
h.startTime = fread(fid, 16, 'uint8');
h.endTime = fread(fid, 16, 'uint8');

% reformat time into 6 element date vectors
h.startTime = bitshift(h.startTime(2:2:16), 8) + h.startTime(1:2:15);
h.startTime(7) = h.startTime(7) + h.startTime(8)/1e3;
h.startTime([3,8]) = []; % delete week day and milliseconds fields

h.endTime = bitshift(h.endTime(2:2:16), 8) + h.endTime(1:2:15);
h.endTime(7) = h.endTime(7) + h.endTime(8)/1e3;
h.endTime([3,8]) = []; % delete week day and milliseconds fields

% user defined
h.user = fread(fid, 224, 'uchar');

% comment
%h.comment = fread(fid, 256, 'uchar');
h.comment = fscanf(fid, '%c', 256);
h.comment = h.comment(h.comment ~= 0);

% find block file type
if ((h.fileType == 11) & (h.fileSubtype == 11) & (h.sizeOf == 2) & ...
    (h.dataType == 12) & (h.refSize ~= 0)),
  % type A, 16-bit differential (with reference) (Vdaq w/ Imager 2001)
  h.dataPrecision = 'uint16';
end

if ((h.fileType == 12) & (h.fileSubtype == 11) & (h.sizeOf == 4) & ...
    (h.dataType == 14) & (h.refSize == 0)),
  % type B, 32-bit DC float (no reference) (Block Convert)
  h.dataPrecision = 'float32';
end

if ((h.fileType == 12) & (h.fileSubtype == 11) & (h.sizeOf == 4) & ...
    (h.dataType == 13) & (h.refSize == 0)),
  % type C, 32-bit integer DC without reference (Vdaq w/ Imager 3001)
  h.dataPrecision = 'uint32';
end

if ((h.fileType == 12) & (h.fileSubtype == 11) & (h.sizeOf == 2) & ...
    (h.dataType == 12) & (h.refSize == 0)),
  % type D, 16-bit integer DC without reference (Vdaq w/ Imager 3001)
  h.dataPrecision = 'uint16';
end
