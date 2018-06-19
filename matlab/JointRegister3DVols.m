function [regVols,params] = JointRegister3D( resizedVols,nblocks)
% imgs = Nrounds cell array per time point
% each element in Nx x Ny x Nz x 5 (4 color channels + DAPI)

if nargin < 2    
    nblocks = [1 1];    
end

nblocksX = nblocks(1);
nblocksY = nblocks(2);

tic;
% part 1: register to last round
ol = 0.05; % 10% overlap
[m,n,z,nCh] = size(resizedVols{1}{1});
B = numel(resizedVols);
stepsize_m = floor(m./nblocksX);
mod_m = mod(m,nblocksX);
stepsize_n = floor(n./nblocksY);
fprintf('Using block size of %f x %f \n', stepsize_m, stepsize_n);
mod_n = mod(n,nblocksY);
% Compute overlap number of pixels
overlap_m = floor(stepsize_m.*ol);
overlap_n = floor(stepsize_n.*ol);
tic;
params = {};
firstRound = resizedVols{1}{5};
for i=2:B    
    msg = sprintf('Registering round [%d / %d] [time=%02f]\n', i, numel(resizedVols)-1, toc);    
    currRound = resizedVols{i}{5};
    size(currRound)
    [params{i}, ~] = DFTRegister3D(single(firstRound),single(currRound),false);
    params{i}
end


regImgs{1} = resizedVols{1};

%% apply registration to sequencing rounds
regVols = cell(7,1);
regVols{1} = resizedVols{1};
for i=2:numel(resizedVols)

    currRound = cell(5,1);
    for j=1:5
        msg = sprintf('Applying registration to round [%d / %d], channel %d [time=%02f]\n', i, numel(resizedVols),j,toc);        
        img = resizedVols{i}{j};
        currRound{j} = DFTApply3D(img, params{i}, false);
    end
    regVols{i} = currRound;
end
%%


if nblocksX > 1 || nblocksY > 1
    blockParams = cell(B-1, nblocksX, nblocksY);    

    for t = 2:B
        fprintf('Processing timepoint %d\n', t);
        % Loop through each block generating the correct subset and registration
        for i = 1:nblocksX
            for j = 1:nblocksY
                % Store m block range
                [bin_m, bin_n] = GetBlocks(i,j);
                % Compute fourier transforms
                block1 = regVols{1}{5}(bin_m, bin_n,:);
                block2 = regVols{t}{5}(bin_m, bin_n,:);
                % Compute offset for block
                [blockParams{t-1,i,j},~] = DFTRegister3D(block1, block2, false);
                %blockParams{t-1,i,j}.shifts = [ blockParams{t-1,i,j}.shifts 0];
                % Store output
                disp(['Shifting (' num2str(i) ' ' num2str(j) ') by ' num2str(blockParams{t-1,i,j}.shifts)]);
            end
        end
    end
end


if nblocksX > 1 || nblocksY > 1
    for t = 2:B
        fprintf('Processing timepoint %d [time=%02f]\n', t, toc);
        % Loop through each block generating the correct subset and registration        
        for kk=1:5
            temp = regVols{t}{kk};
            for i = 1:nblocksX
                for j = 1:nblocksY
                    [bin_m, bin_n] = GetBlocks(i,j);
                    % Compute fourier transforms
                    % Compute offset for block
                    temp(bin_m, bin_n,:) = DFTApply3D(temp(bin_m, bin_n,:), blockParams{t-1,i,j});                  
                end
            end
            regVols{t}{kk} = temp;
        end
    end    
end


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
