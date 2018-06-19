function RegisterNissl( basePath, prefix, Nround, dapiRound )
%locs = ManuallyIdentifyCells(outPath,Nround);
allImgs = LoadMaxProjImages(basePath, prefix, 1, Nround, []);

% register dapi round to image
% register raw sequencing data -> registered
%
[X Y Nch Nround] = size(allImgs);
outPath = fullfile(basePath,'output');
disp('Loading images');
regRound = imread(fullfile(outPath, 'dftreg', sprintf('f3_T_%02d_C_04.tif', dapiRound-1)));
seqRound = squeeze(allImgs(:,:,5,dapiRound));

% register nissl dapi to seq dapi

% register nissl -> raw sequencing data
dapiSeq = imread(fullfile(outPath, sprintf('dapi_round%d.tif', dapiRound)));
[dX, dY] = size(dapiSeq);

dapiSeq = dapiSeq(1:X,1:Y);

dapiSeg = imread(fullfile(outPath, 'seg_dapi.tif'));
dapiSeg = dapiSeg(1:X,1:Y);

nisslSeg = imread(fullfile(outPath, 'seg_nissl.tif'));
nisslSeg = nisslSeg(1:X, 1:Y);
% apply to dapi

disp('Computing transforms');
regSegToReg = dftregistration(fft2(regRound), fft2(seqRound));
regSegToSeq = dftregistration(fft2(dapiSeq), fft2(dapiSeg));

% apply transform nissl -> raw -> registered
disp('applying transform');
temp = dftapply(fft2(dapiSeg), regSegToSeq);
registeredDapi = dftapply(fft2(temp), regSegToReg);
temp = dftapply(fft2(nisslSeg), regSegToSeq);
registeredNissl = dftapply(fft2(temp), regSegToReg);

%
disp('saving files');
imwrite(uint16(registeredDapi), fullfile(outPath, 'registered_seg_dapi.tif'))
imwrite(uint16(registeredNissl), fullfile(outPath, 'registered_seg_nissl.tif'));

end

