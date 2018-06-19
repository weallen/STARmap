function [params, regImg] = DFTRegister2D(fixedImage, movingImage, preFFT)
if nargin < 3
    preFFT = false;
end
[nr,nc] = size(movingImage);
Nr = ifftshift(-fix(nr/2):ceil(nr/2)-1);
Nc = ifftshift(-fix(nc/2):ceil(nc/2)-1);

%
if preFFT
    CC = ifft2(fixedImage .* conj(movingImage));
else
    CC = ifft2(fftn(fixedImage) .* conj(fftn(movingImage)));
end
CCabs = abs(CC);

[~,ix] = max(CCabs(:));
[i,j,k] = ind2sub(size(CCabs), ix);%ftr = fftFixed .* fftMoving;
CCmax = CC(i,j);
diffphase = angle(CCmax);

rowShift = Nr(i); colShift = Nc(j); 

params = struct();
params.shifts = [rowShift, colShift];
params.diffphase = diffphase;

if nargout > 1   
    regImg = DFTApply2D(movingImage, params,preFFT);            
end

