function [allImgs, images] = LoadMaxProjImages( input_dir, prefix_name, downscale, Nround, Zrange )
images = cell(Nround,1);
for i=1:Nround
    
    fprintf('Loading round %d...\n', i);
    
    dirname = fullfile(input_dir, ['round' num2str(i)]);
    fnames = cell(5,1);
    for j=1:4 % 4 colors without DAPI
        fnames{j} = fullfile(dirname, sprintf('%sseq%d_decon_ch0%d.tif', prefix_name, i, j-1));
    end
    currImg = LoadMultipageTiff(fnames{1});
    if ~isempty(Zrange)
        disp(['Using Z range ' num2str(Zrange{i})]);
        currImg = currImg(:,:,Zrange{i});
    end
    firstFrame = max(currImg,[],3);
    if downscale < 1
        firstFrame = imresize(firstFrame, downscale);
    end
    currImgs = zeros(size(firstFrame,1), size(firstFrame,2), 4);
    currImgs(:,:,1) = firstFrame;
    for j=2:4        
        currFrame = LoadMultipageTiff(fnames{j});
        if ~isempty(Zrange)
            currFrame = currFrame(:,:,Zrange{i});
        end
        temp = max(currFrame,[],3);
        if downscale < 1
            temp = imresize(temp, downscale);
        end
        currImgs(:,:,j) = temp;
    end
    images{i} = currImgs;
end


%

% collapse to common sized array
maxX = 1E10; maxY = 1E10;
for i=1:Nround
   [currX,currY,~] = size(images{i}); 
   if currX < maxX
       maxX = currX;
   end
   if currY < maxY
       maxY = currY;
   end
end

fprintf('Collapsed to size %d by %d\n', maxX, maxY);
allImgs = uint16(zeros(maxX, maxY, 5, Nround)); % 5th channel is max projection
for i=1:Nround
    temp = images{i};
    allImgs(:,:,1:4,i) = temp(1:maxX, 1:maxY, :);
    allImgs(:,:,5,i) = squeeze(max(allImgs(:,:,1:4,i),[],3));
end


