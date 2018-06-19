function  SaveRegistered(dirName, imgs )

%SAVECOLORMAX Summary of this function goes here
%   Detailed explanation goes here
for i=1:numel(imgs)
    i
    temp = uint8(imgs{i});
    for ch=1:4
        fname = fullfile(dirName, sprintf('reg_round%02d_ch%02d.tif',i,ch));
        if exist(fname, 'file') == 2
            delete(fname);
        end
        for j=1:size(temp,3)        
            imwrite(temp(:,:,j,ch), fname, 'writemode', 'append');        
        end
    end
end

