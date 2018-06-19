function [colorCodes, colorMap] = GenerateBarcodeColorSeq( Nround )
Nbarcodes = 4^Nround;
bitSeqs = dec2bin(0:(Nbarcodes-1));
colorCodes = zeros(Nbarcodes, Nround);
colorMap = containers.Map;
keys = cell(Nbarcodes,1);
vals = zeros(Nbarcodes,1);
for i=1:Nbarcodes
    k = 1;
    for j=1:2:size(bitSeqs,2)
        colorCodes(i,k) = Bin2Base(bitSeqs(i,j:(j+1)));
        k = k + 1;
    end
    keys{i} = num2str(colorCodes(i,:));
    vals(i) = i;
end
colorMap = containers.Map(keys, vals);

function base = Bin2Base(val)
if strcmp(val, '00')
    base = 0;
elseif strcmp(val, '01')
    base = 1;
elseif strcmp(val, '10')
    base = 2;
elseif strcmp(val, '11')
    base = 3;
end
    


