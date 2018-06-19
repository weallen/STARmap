function [params, regImg] = DFTRegister3D(fixedVolume, movingVolume, preFFT)
if nargin < 3
    preFFT = false;
end



[nr,nc,nz] = size(movingVolume);
Nr = ifftshift(-fix(nr/2):ceil(nr/2)-1);
Nc = ifftshift(-fix(nc/2):ceil(nc/2)-1);
Nz = ifftshift(-fix(nz/2):ceil(nz/2)-1);

%
if preFFT
    CC = ifftn(fixedVolume .* conj(movingVolume));
else
    CC = ifftn(fftn(fixedVolume) .* conj(fftn(movingVolume)));
end
CCabs = abs(CC);

[~,ix] = max(CCabs(:));
[i,j,k] = ind2sub(size(CCabs), ix);%ftr = fftFixed .* fftMoving;
CCmax = CC(i,j,k);
diffphase = angle(CCmax);

rowShift = Nr(i); colShift = Nc(j); zShift = Nz(k);

params = struct();
params.shifts = [rowShift, colShift, zShift];
params.diffphase = diffphase;

if nargout > 1   
    regImg = DFTApply3D(movingVolume, params,preFFT);            
end