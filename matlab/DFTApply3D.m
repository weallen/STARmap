function reg = DFTApply3D(movingVolume, params,preFFT)
if nargin < 3
    preFFT = false;
end

[nr,nc,nz] = size(movingVolume);
Nr = ifftshift(-fix(nr/2):ceil(nr/2)-1);
Nc = ifftshift(-fix(nc/2):ceil(nc/2)-1);
Nz = ifftshift(-fix(nz/2):ceil(nz/2)-1);

[Nc,Nr,Nz] = meshgrid(Nc,Nr,Nz);
if preFFT
    reg = movingVolume.*exp(1i.*2.*pi.*(-params.shifts(1).*Nr/nr-params.shifts(2).*Nc/nc-params.shifts(3).*Nz/nz));
else
    reg = fftn(movingVolume).*exp(1i.*2.*pi.*(-params.shifts(1).*Nr/nr-params.shifts(2).*Nc/nc-params.shifts(3).*Nz/nz));
end
reg = abs(ifftn(reg*exp(1i*params.diffphase)));

% set wrapped to 0
if params.shifts(1) > 0
    reg(1:ceil(params.shifts(1)),:,:) = 0;
else
    reg(nr+floor(params.shifts(1)):nr,:,:) = 0;
end
if params.shifts(2) > 0
    reg(:,1:ceil(params.shifts(2)),:) = 0;
else
    reg(:,nc+floor(params.shifts(2)):nc,:) = 0;
end
if params.shifts(3) > 0
    reg(:,:,1:ceil(params.shifts(3))) = 0;
else
    reg(:,:,nz+floor(params.shifts(3)):nz) = 0;
end
end

