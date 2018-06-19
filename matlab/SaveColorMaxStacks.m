function SaveColorMaxStacks( input_dir, prefix_name, output_dir,Nround)
images = cell(Nround,1);
for i=1:Nround
    fprintf('Loading round %d...\n', i);
    dirname = fullfile(input_dir, ['round' num2str(i)]);
    fnames = cell(5,1);
    for j=1:4 % 4 colors without DAPI
        fnames{j} = fullfile(dirname, sprintf('%sseq%d_decon_ch0%d.tif', prefix_name, i, j-1));
    end
    firstFrame = LoadMultipageTiff(fnames{1});
    [X,Y,Z] = size(firstFrame);
    currImgs = zeros(size(firstFrame,1), size(firstFrame,2), Z, 4);
    currImgs(:,:,:,1) = firstFrame;
    for j=2:4        
        currFrame = LoadMultipageTiff(fnames{j});        
        currImgs(:,:,:,j) = currFrame;
    end
    finalImg = max(currImgs,[],4);
    %t = Tiff(fullfile(output_dir, ['stack_round' num2str(i) '.tif']),'w');
    disp('Saving ');
    for z=1:Z        
        imwrite(uint16(finalImg(:,:,z)), fullfile(output_dir, ['stack_round' num2str(i) '.tif']), 'tif', 'WriteMode', 'append', 'compress','none');
    end
    %t.close();
end


%
