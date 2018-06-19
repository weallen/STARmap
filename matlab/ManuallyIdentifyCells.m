function cellLocs = ManuallyIdentifyCells( basePath, useReg )

if nargin < 2
    useReg = true;
end


if useReg
    avgImg = imread(fullfile(basePath, 'output', 'registered_seg_dapi.tif'));
else
    avgImg = imread(fullfile(basePath, 'output', 'seg_dapi.tif'));
end
imagesc(avgImg); colormap gray;
cellLocs = [];
while true
    
    [x,y] = ginput(1);
    cellLocs = [cellLocs; x y]
    if x < size(avgImg,2) && y < size(avgImg,1)
        imagesc(avgImg); colormap gray; hold on;        
        plot(cellLocs(:,1),cellLocs(:,2),'r.');        
    else
        break;
    end    
end

end

