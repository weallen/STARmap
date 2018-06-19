function [colorSeq] = ExtractColorSeq( imgs, points, ksize )
%EXTRACTCOLORSEQ Extract color sequence for each point across sequencing
%rounds and color channels
% ksize is [x,y,z] vector
if nargin < 3
    ksize = [4 4 4];
end
[X Y Z nCh] = size(imgs{1});
Npoints = size(points,1);
Nrounds = numel(imgs);
colorSeq = zeros(Npoints, Nrounds, 4);

if size(imgs{1},3) == 1
    for r=1:Nrounds
        fprintf('Processing round %d\n', r);
        currRound = imgs{r};
        for i=1:Npoints            
            p = points(i,:);
            extentsX = GetExtents(p(2), ksize(1), X);
            extentsY = GetExtents(p(1), ksize(2), Y);
            currVol = squeeze(currRound(extentsX,...
                extentsY, 1,:));
            for j=1:4
                if numel(extentsX) > 1 && numel(extentsY) > 1
                    temp = currVol(:,:,j);
                    colorSeq(i,r,j) = squeeze(mean(temp(:)));
                end
            end
        end        
    end
else        
    upd = textprogressbar(Npoints);
    for i=1:Npoints
        if mod(i,100) == 0            
            upd(i);
        end

        p = points(i,:);
        extentsX = GetExtents(p(2), ksize(1), X);
        extentsY = GetExtents(p(1), ksize(2), Y);                    
        extentsZ = GetExtents(p(3), ksize(3), Z);    
        for r=1:Nrounds    
            
            currRound = imgs{r};        
            
            currVol = currRound(extentsX, extentsY, extentsZ, :);
            if size(currVol,3) == 1
                colorSeq(i,r,:) = single(squeeze(sum(sum(currVol))));
            else
                colorSeq(i,r,:) = single(squeeze(sum(sum(sum(currVol)))));
            end
        end
        for r=1:Nrounds
            colorSeq(i,r,:) = colorSeq(i,r,:) ./ (sqrt(sum(squeeze(colorSeq(i,r,:)).^2)) + 1E-6);
        end
    end
end
fprintf('\n');

function e = GetExtents(pos, ksize, lim)
if pos-ksize < 1 
    e1 = 1;
else
    e1 = pos-ksize;
end
if pos+ksize > lim
    e2 = lim;
else
    e2 = pos+ksize;
end
e = e1:e2;

