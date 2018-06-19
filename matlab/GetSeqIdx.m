function [ whichBase ] = GetSeqIdx( bases, target )
whichBase = zeros(numel(bases),1);
for i=1:numel(bases)
    if strcmp(bases{i}, target)
        whichBase(i) = 1;
    end
end


end

