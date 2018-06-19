function reg = DFTApply2D(movingImage, params,preFFT)
if nargin < 3
    preFFT = false;
end

[nr,nc] = size(movingImage);
Nr = ifftshift(-fix(nr/2):ceil(nr/2)-1);
Nc = ifftshift(-fix(nc/2):ceil(nc/2)-1);

[Nc,Nr] = meshgrid(Nc,Nr);
if preFFT
    reg = movingImage.*exp(1i.*2.*pi.*(-params.shifts(1).*Nr/nr-params.shifts(2).*Nc/nc));
else
    reg = fft2(movingImage).*exp(1i.*2.*pi.*(-params.shifts(1).*Nr/nr-params.shifts(2).*Nc/nc));
end
reg = abs(ifft2(reg*exp(1i*params.diffphase)));

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
end