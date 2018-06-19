function  imgs = SaveRegistered(dirName, Nround)

%SAVECOLORMAX Summary of this function goes here
%   Detailed explanation goes here
info = imfinfo(fullfile(dirName, sprintf('reg_round%02d_ch%02d.tif',1,1)));
Y = info(1).Width; X = info(1).Height; Z = numel(info);
nCh = 4;
imgs = cell(Nround, 1);
for i=1:Nround
    i
    temp = zeros(X,Y,Z,nCh,'uint8');
    for ch=1:4
        for k=1:Z        
            temp(:,:,k,ch) = uint8(imread( fullfile(dirName, sprintf('reg_round%02d_ch%02d.tif',i,ch)), k));
        end
    end
    imgs{i} = temp;
end
