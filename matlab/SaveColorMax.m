function SaveColorMax(dirName, imgs )
%SAVECOLORMAX Summary of this function goes here
%   Detailed explanation goes here
for i=1:numel(imgs)
    i
    temp = uint8(max(imgs{i},[],4));
    fname = fullfile(dirName, sprintf('reg_round%02d.tif',i));
    if exist(fname, 'file') == 2
        delete(fname);
    end
    for j=1:size(temp,3)        
        imwrite(temp(:,:,j), fname, 'writemode', 'append');        
    end
end

