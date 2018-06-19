function [imgs,avgReg,params] = JointRegister3D( imgs,nblocks,useAverage)
% given cell array of images, register multiple images to reference (last
% image in array)

if nargin < 2    
    nblocks = [1 1];    
end

nblocksX = nblocks(1);
nblocksY = nblocks(2);

if nargin < 3
    useAverage = false;
end

%{
for i=1:numel(imgs)
    temp = imgs{i};
    imgs{i} = temp(1:1000,1:1000,:,:);
end
%}
tic;
% part 1: register to last round
ol = 0.1; % 10% overlap
colorMax = cell(numel(imgs),1);
regImgs = cell(numel(imgs),1);
[m,n,z,nCh] = size(imgs{1});
B = numel(imgs);
stepsize_m = floor(m./nblocksX);
mod_m = mod(m,nblocksX);
stepsize_n = floor(n./nblocksY);
fprintf('Using block size of %f x %f \n', stepsize_m, stepsize_n);
mod_n = mod(n,nblocksY);
% Compute overlap number of pixels
overlap_m = floor(stepsize_m.*ol);
overlap_n = floor(stepsize_n.*ol);


%for i=1:B
%    colorMax{i} = max(imgs{i},[],4);
%end

% compute average image
%if useAverage
%    avgReg = colorMax{i};
%    for i=2:numel(colorMax)
%        avgReg = avgReg + colorMax{i};
%    end
%    avgReg = avgReg/numel(colorMax);
%    avgReg = fftn(avgReg);
%end

regImgs{1} = uint8(max(imgs{1},[],4));

% register to first round
%for i=1:numel(colorMax)
%    msg = sprintf('Precomputing FFT %d [time=%02f]...\n', i, toc);
%    if i>1
%        fprintf(repmat('\b',1,numel(msg)));
%    end
%    fprintf(msg);
    %colorMax{i} = fft2(max(max(imgs{i},[],4),[],3));
%    colorMax{i} = fftn(colorMax{i});        
%end
fprintf('%s\n', repmat('-',1,20));
if useAverage
    params = cell(numel(imgs),1);

    newAvgReg = zeros(size(colorMax{1}));
    for i=1:B
       msg = sprintf('Registering round [%d / %d] [time=%02f]\n', i, numel(imgs)-1, toc);
        fprintf(msg);      
    
        [params{i}, currReg] = DFTRegister3D(single(avgReg),colorMax{i},true);
        %params{i-1}.shifts = -params{i-1}.shifts;;
        disp(['Shifting by ' num2str(params{i}.shifts)]);
        %p.shifts
        regImgs{i} = uint8(currReg);
        newAvgReg = newAvgReg + currReg;
    end
    avgReg = newAvgReg / B;
else
    params = cell(numel(imgs)-1,1);
    %round1 = uint8(max(max(imgs{1},[],4),[],3));
    round1 = uint8(max(imgs{1},[],4));
    for i=2:B
        roundN = uint8(max(imgs{i},[],4));    
        msg = sprintf('Registering round [%d / %d] [time=%02f]\n', i, numel(imgs)-1, toc);
        fprintf(msg);          
        [params{i-1}, currReg] = DFTRegister3D(round1,roundN,false);
        %params{i-1}.shifts = -params{i-1}.shifts;;
        disp(['Shifting by ' num2str(params{i-1}.shifts)]);
        %p.shifts
        regImgs{i} = uint8(currReg);
    end
end


if nblocksX > 1 || nblocksY > 1
    if useAverage
        blockParams = cell(B, nblocksX, nblocksY);    
        for t = 1:B
            fprintf('Processing timepoint %d\n', t);
            % Loop through each block generating the correct subset and registration
            for i = 1:nblocksX
                for j = 1:nblocksY
                    % Store m block range
                    [bin_m, bin_n] = GetBlocks(i,j);
                    % Compute fourier transforms
                    block1 = avgReg(bin_m, bin_n,:);
                    block2 = regImgs{t}(bin_m, bin_n,:);
                    % Compute offset for block
                    [blockParams{t,i,j},~] = DFTRegister3D(block1, block2, false);
                    %blockParams{t,i,j}.shifts = [ blockParams{t,i,j}.shifts 0];
                    % Store output
                    disp(['Shifting (' num2str(i) ' ' num2str(j) ') by ' num2str(blockParams{t,i,j}.shifts)]);
                end
            end
        end        
    else
        blockParams = cell(B-1, nblocksX, nblocksY);    

        for t = 2:B
            fprintf('Processing timepoint %d\n', t);
            % Loop through each block generating the correct subset and registration
            for i = 1:nblocksX
                for j = 1:nblocksY
                    % Store m block range
                    [bin_m, bin_n] = GetBlocks(i,j);
                   % disp(size(regImgs{1}));
                    %disp(size(regImgs{t}));
                    % Compute fourier transforms
                    block1 = regImgs{1}(bin_m, bin_n,:);
                    block2 = regImgs{t}(bin_m, bin_n,:);
                    % Compute offset for block
                    [blockParams{t-1,i,j},~] = DFTRegister3D(block1, block2, false);
                    %blockParams{t-1,i,j}.shifts = [ blockParams{t-1,i,j}.shifts 0];
                    % Store output
                    disp(['Shifting (' num2str(i) ' ' num2str(j) ') by ' num2str(blockParams{t-1,i,j}.shifts)]);
                end
            end
        end
    end
end

if useAverage
    for i=1:numel(imgs)
        temp = imgs{i};
        for j=1:4        
            msg = sprintf('Applying registration to round [%d / %d], channel %d [time=%02f]\n', i, numel(imgs),j,toc);
            fprintf(msg);        
            temp(:,:,:,j) = uint8(DFTApply3D(temp(:,:,:,j), params{i}, false));
        end        
        imgs{i} = temp;
    end
else
    for i=2:numel(imgs)
        temp = imgs{i};
        for j=1:4        
            msg = sprintf('Applying registration to round [%d / %d], channel %d [time=%02f]\n', i, numel(imgs),j,toc);
            fprintf(msg);        
            temp(:,:,:,j) = uint8(DFTApply3D(temp(:,:,:,j), params{i-1}, false));
        end        
        imgs{i} = temp;
    end
end

if nblocksX > 1 || nblocksY > 1
    if useAverage
        for t = 1:B
            fprintf('Processing timepoint %d [time=%02f]\n', t, toc);
            % Loop through each block generating the correct subset and registration
            temp = imgs{t};
            for kk=1:4
                for i = 1:nblocksX
                    for j = 1:nblocksY
                        [bin_m, bin_n] = GetBlocks(i,j);
                        % Compute fourier transforms
                        % Compute offset for block
                        temp(bin_m, bin_n,:,kk) = uint8(DFTApply3D(temp(bin_m, bin_n,:,kk), blockParams{t,i,j}));                  
                    end
                end
            end
            imgs{t} = temp;
        end            
    else
        for t = 2:B
            fprintf('Processing timepoint %d [time=%02f]\n', t, toc);
            % Loop through each block generating the correct subset and registration
            temp = imgs{t};
            for kk=1:4
                for i = 1:nblocksX
                    for j = 1:nblocksY
                        [bin_m, bin_n] = GetBlocks(i,j);
                        % Compute fourier transforms
                        % Compute offset for block
                        temp(bin_m, bin_n,:,kk) = uint8(DFTApply3D(temp(bin_m, bin_n,:,kk), blockParams{t-1,i,j}));                  
                    end
                end
            end
            imgs{t} = temp;
        end    
    end
end

% recompute average
avgReg = zeros([m,n,z]);

for i=1:B
    avgReg = avgReg + single(max(imgs{i},[],4));
end
avgReg = avgReg / B;


function [bin_m, bin_n] = GetBlocks(i,j)
    % Store m block range
    if nblocksX > 1
        if i == 1
            bin_m = 1:(stepsize_m + overlap_m);
        elseif i == nblocksX
            bin_m = (stepsize_m.*(i-1) - overlap_m):(stepsize_m.*i + mod_m);
        else
            bin_m = (stepsize_m.*(i-1) - overlap_m):(stepsize_m.*i + overlap_m);
        end
    elseif nblocksX == 1
        bin_m = 1:stepsize_m;
    end
    if nblocksY > 1
        % Store n block range
        if j == 1
            bin_n = 1:(stepsize_n + overlap_n);
        elseif j == nblocksY
            bin_n = (stepsize_n.*(j-1) - overlap_n):(stepsize_n.*j + mod_n);
        else
            bin_n = (stepsize_n.*(j-1) - overlap_n):(stepsize_n.*j + overlap_n);
        end
    elseif nblocksY == 1
        bin_n = 1:stepsize_n;
    end
end
end


%%
%}

% part 1: 2D registration of max projections
%{
regImgs{1} = max(max(imgs{1},[],4),[],3);
avgReg = regImgs{1};
for i=1:numel(imgs)
    regImgs{i}  = max(max(imgs{i},[],4),[],3);
    avgReg = avgReg + regImgs{i};    
end
avgReg = avgReg / numel(imgs);

for i=1:niter
    for kk=1:numel(regImgs)
        msg = sprintf('Precomputing FFT %d [time=%02f]...\n', kk, toc);
        if kk>1
            fprintf(repmat('\b',1,numel(msg)));
        end
        fprintf(msg);
        %regImgs{kk} = fft2(regImgs{kk});
    end
    
    fprintf('%s\n', repmat('-',1,20));
    fprintf('Iteratively registering to average round %d\n', i);
    %avgReg = fft2(avgReg);
    newAvgReg = zeros(size(avgReg));
    totalShift = 0;
    for j=1:numel(regImgs)
        msg = sprintf('\t--> Registering round [%d / %d] [time=%02f]\n', j, numel(imgs), toc);
        fprintf(msg);      
        [params{j}, currReg] = DFTRegister2D(avgReg, regImgs{j}, false);
        params{j}.shifts = [params{j}.shifts 0];
        disp(['Current shift ' num2str(params{j}.shifts)]);
        totalShift = totalShift + sum(abs(params{j}.shifts));
        newAvgReg = newAvgReg + currReg;
        regImgs{j} = currReg;
    end
    avgReg = newAvgReg / numel(imgs);
    fprintf('TOTAL SHIFT %d\n', totalShift);
    if totalShift == 0
        break;
    end    
end
%}

%}
%}
% part 2: Full 3D registration
% part 2: register to average multiple times
%{
regImgs{1} = max(imgs{1},[],4);
avgReg = regImgs{1};
for i=2:numel(imgs)
    regImgs{i}  = max(imgs{i},[],4);
    avgReg = avgReg + regImgs{i};    
end
avgReg = avgReg / numel(imgs);
%}
%{
for i=1:niter 
    for kk=1:numel(regImgs)
        msg = sprintf('Precomputing FFT %d [time=%02f]...\n', kk, toc);
        if kk>1
            fprintf(repmat('\b',1,numel(msg)));
        end
        fprintf(msg);
        regImgs{kk} = fftn(regImgs{kk});
    end
    
    fprintf('%s\n', repmat('-',1,20));
    fprintf('Iteratively registering to average round %d\n', i);
    avgReg = fftn(avgReg);
    newAvgReg = zeros(size(avgReg));
    totalShift = 0;
    for j=1:numel(regImgs)
        msg = sprintf('\t--> Registering round [%d / %d] [time=%02f]\n', j, numel(imgs), toc);
        fprintf(msg);      
        [params{j}, currReg] = DFTRegister3D(avgReg, regImgs{j}, true);
        params{j}.shifts
        totalShift = totalShift + sum(abs(params{j}.shifts));
        newAvgReg = newAvgReg + currReg;
        regImgs{j} = currReg;
    end
    avgReg = newAvgReg / numel(imgs);
    if totalShift == 0
        break;
    end    
end

%}
%{
for i=2:numel(imgs)
    temp = imgs{i};
    for j=1:4        
        msg = sprintf('Applying registration to round [%d / %d], channel %d [time=%02f]\n', i, numel(imgs),j,toc);
        fprintf(msg);
        temp(:,:,:,j) = DFTApply3D(temp(:,:,:,j), params{i-1}, false);
    end        
    imgs{i} = temp;
end
%}


