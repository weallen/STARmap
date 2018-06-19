function [ imgStack ] = EqualizeHist3D( imgStack )


Nround = numel(imgStack);
for t=1:Nround 
    fprintf('Equalizing round %d\n', t);
    currStack = uint8(imgStack{t});        
    for i=1:4                
        i
        currStack(:,:,:,i) = uint8(imhistmatchn(uint8(currStack(:,:,:,i)), uint8(currStack(:,:,:,1))));
        %currStack(:,:,:,i) = gather(histeq(gpuArray(currStack(:,:,:,i))));
    end    
    imgStack{t} = currStack;
end

end

