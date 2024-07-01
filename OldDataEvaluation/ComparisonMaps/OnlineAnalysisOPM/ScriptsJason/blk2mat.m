% batch import/conversion of optical imaging (.BLK) data files

% 24/6/2012 - Shaun L Cloherty <s.cloherty@ieee.org>

fprintf(1,'WARNING: This script converts .BLK files to .mat files, consolidating blocks\n');
fprintf(1,'         (.BLK files) within each experiment. It reads and converts all .BLK\n');
fprintf(1,'         files in, and writes .mat files to, the current working directory.\n\n');

r = input('Any existing .mat files will be overwritten. Do you wish to proceed? [Y/N/n]: ','s');

if ~strcmp(r,'Y'),
  return
end

% srcPath = '/Volumes/NVRI/IN VIVO LAB/Data/Cat Data/CatAB/OI';
srcPath = 'U:\WallabyA\Incoming';
% dstPath = '~/Documents/Data/Optical Imaging/Optical Imaging @ NVRI/M1378/';
dstPath = 'U:\WallabyA\Optical Imaging';

d = dir(fullfile(srcPath, '*.BLK'));

% extract expPrefix, expId and blkId from each filename 
finfo = arrayfun(@(x) regexp(x.name,'(?<expPrefix>.*)_E(?<expId>\d+)B(?<blkId>\d+)\.BLK','tokens'), d);

for expPrefix = unique(cellfun(@(x) x(1), finfo))',
  i = find(cellfun(@(x) strcmp(x(1),expPrefix),finfo));

  expIds = cellfun(@(x) str2num(x{2}), finfo(i));
  blkIds = cellfun(@(x) str2num(x{3}), finfo(i));

%   for expId = intersect(unique(expIds)', [18:26]),
%   for expId = setdiff(unique(expIds)', [18:26]),
  for expId = unique(expIds)',
    j = find(expIds == expId);

    [~,k] = sort(blkIds(j)); % try and preserve order of acquisition...
  
    tmp = expPrefix{1}; tmp(1) = lower(tmp(1));
    matFile = sprintf('%sexp%i.mat',tmp,expId);

    fprintf(1,'Exp. %i, %i blks --> %s... ', expId, length(j), matFile);

    % note: we do the import in n blocks of m files to avoid any limit placed
    %       on the number of files which may be open simultaneously
    m = min(length(j),50);
  
    n = floor(length(j)/m);
    for p = 0:n,
      ii = p*m + 1;
      jj = min((p+1)*m, length(j));
      fids = arrayfun(@(x) fopen(fullfile(srcPath,x.name),'r'), d(j(k(ii:jj))));
    
      try,
        images_ = blkImport(fids);
        cnt_ = sum(cellfun(@(x) ~isempty(x), images_),2);
      catch,
        error('Call to blkImport() failed!');
      end
    
      fclose('all');
    
      if p == 0,
        images = images_;
        cnt = cnt_;
      else,
        if (length(cnt) < length(cnt_)),
          cnt(length(cnt_)) = 0; % pad with zeros
        end
        
        for q = 1:length(cnt_),
          images(q,cnt(q)+[1:cnt_(q)]) = images_(q,1:cnt_(q));
        end
      
      
        cnt(1:length(cnt_)) = cnt(1:length(cnt_)) + cnt_;
      end
      clear ii jj fids images_ cnt_ q
    end

    s = whos('images');
  
    saveVer = '';
    if (log(s.bytes)/log(2) >= 31), % 2^31 = 2Gb
      saveVer = '-v7.3';
    end
  
    try,
      save(fullfile(dstPath,matFile),'images',saveVer);
    catch,
      error(sprintf('Failed to save ''images'' to %s.', matFile));
    end
    fprintf(1,'Ok\n');

    clear j k tmp matFile m n images cnt p s saveVer
  end
  
  clear i expIds blkIds expId tmp
end

clear r srcPath dstPath d finfo expPrefix
