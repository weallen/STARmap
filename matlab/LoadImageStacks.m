function images = LoadImageStacks( input_dir, prefix_name, Nround, zrange)
%LOADIMAGESTACKS Load image stacks for each round
warning('off','all');
if nargin < 4
    zrange = [];
end

images = cell(Nround,1);
for i=1:Nround
    
    fprintf('Loading round %d...\n', i);
    
    dirname = fullfile(input_dir, ['round' num2str(i)]);
    fnames = cell(5,1);
    for j=1:4 % 4 colors without DAPI
        fnames{j} = fullfile(dirname, sprintf('%sseq%d_decon_ch0%d.tif', prefix_name, i, j-1));
    end
    currImg = uint8(LoadMultipageTiff(fnames{1}));
    if ~isempty(zrange)
        currImg = currImg(:,:,zrange);
    end
    currImgs = zeros(size(currImg,1), size(currImg,2), size(currImg,3), 4,'uint8');
    currImgs(:,:,:,1) = currImg;
    for j=2:4        
        currFrame = uint8(LoadMultipageTiff(fnames{j}));
        if isempty(zrange)
            currImgs(:,:,:,j) = currFrame;
        else
            currImgs(:,:,:,j) = currFrame(:,:,zrange);
        end
    end
    images{i} = currImgs;
end

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



fprintf('Collapsed to size %d by %d by %d\n', maxX, maxY, size(images{1},3));
for i=1:Nround
    fprintf('Collapsing round %d\n', i);
    temp = images{i};
    images{i} = temp(1:maxX, 1:maxY, :,:);
end


end

