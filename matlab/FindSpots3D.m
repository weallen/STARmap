function centroids = FindSpots3D( imgstack, varargin)
%FINDSPOTS3D Find spots in single plane
useGPU = false;
p = inputParser();
p.addParameter('thresholdRel', 0.2, @isnumeric);
p.addParameter('ksize', 3, @isnumeric);
p.addParameter('kdims', [5 5], @isnumeric);
p.parse(varargin{:});

thresholdRel = p.Results.thresholdRel;
ksize = p.Results.ksize;
kdims = p.Results.kdims;

%colorMax = max(imgstack,[],4);
imgstack = single(imgstack);
%imgstack = (imgstack - mean(imgstack(:)))./std(imgstack(:)); % z score intensities
%colorMax = colorMax(500:1000, 500:1000, :);
thresholdRel = thresholdRel * max(imgstack(:));
fprintf('Using threshold value of %d\n', thresholdRel);
[X,Y,Z,nCh] = size(imgstack);

h = fspecial('log',kdims, ksize);
if useGPU
    filteredImg = gpuArray(single(zeros(size(imgstack))));
else
     filteredImg = single(zeros(size(imgstack)));
end
for i=1:Z
    msg = sprintf('\t--> Filtering image [Z=%02d / %02d]\n', i, Z);
    if i>1
        fprintf(repmat('\b',1,numel(msg)));
    end
    fprintf(msg);
    % apply LoG filtering to each plane
    if useGPU
        filteredImg(:,:,i) = imfilter(gpuArray(single(imgstack(:,:,i))), h);
    else
        filteredImg(:,:,i) = imfilter(single(imgstack(:,:,i)), h);
    end
end


rmax = boolean(zeros(size(imgstack)));

for i=1:Z
    msg = sprintf('\t--> Computing regional min [Z=%02d / %02d]\n', i, Z);
    if i>1
        fprintf(repmat('\b',1,numel(msg)));
    end
    fprintf(msg);
    rmax(:,:,i) = boolean(gather(imregionalmin(filteredImg(:,:,i))));
end

%rmax = boolean(imregionalmin(gather(filteredImg)));


for i=1:Z
    rmax(:,:,i) = rmax(:,:,i) & (imgstack(:,:,i) > thresholdRel);
end
fprintf(1, '\t--> Finding centroids\n');
CC = bwconncomp(rmax);
centroidProps = regionprops(CC, 'centroid');

centroids = zeros(length(centroidProps), 3,'int16');
if numel(centroidProps(1).Centroid) == 2
    parfor i=1:length(centroidProps)
        centroids(i,:) = [int16(centroidProps(i).Centroid) 0];
    end
else    
    for i=1:length(centroidProps)    
        centroids(i,:) = int16(centroidProps(i).Centroid);
    end
end
temp = centroids(centroids(:,3)==10,:);
disp('done');

end

