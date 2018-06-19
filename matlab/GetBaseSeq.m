function [bases, maxColors, baseScores] = GetBaseSeq( colorSeq, startBase)
%GETBASESEQ Summary of this function goes here
%   Detailed explanation goes here
[Npoint,Nround,Nch] = size(colorSeq);
maxColors = zeros(Npoint, Nround);
baseScores = zeros(Npoint, Nround);
% determine maximal color for each round
prb = textprogressbar(Npoint);
for i=1:Npoint
    if mod(i,500) == 0
        prb(i);
    end
    for j=1:Nround
        currMax = max(colorSeq(i,j,:),[],3);
        if ~isnan(currMax)
            m = find(colorSeq(i,j,:) == currMax); 
            % choose first in case of ties
            maxColors(i,j) = m(1);
            baseScores(i,j) = -log(currMax);
        else
            maxColors(i,j) = -1;
            baseScores(i,j) = Inf;
        end
    end    
end
fprintf('\n');
bases = {};
%prb = textprogressbar(Npoint);

% do alexa decoding
parfor i=1:Npoint
   % if mod(i,500) == 0
   %     prb(i);
   % end
    
    if ~any(isinf(baseScores(i,:)))
        %bases{i} = cell2mat(DecodeBases(maxColors(i,:),startBase));            
        bases{i} = Colorseq2Str(maxColors(i,:));
    end
end
fprintf('\n');

