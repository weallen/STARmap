function [bases] = DecodeBases( maxColors, startBase)
%GETBASESEQ Summary of this function goes here
%   Detailed explanation goes here
[Npoint,Nround] = size(maxColors);
% do alexa decoding
bases = cell(Npoint,1);
%prb = textprogressbar(Npoint);
k = {1,2,3,4};
v = {{'AA','CC','GG','TT'},...
    {'AC','CA','GT','TG'},....
    {'AG','CT','GA','TC'},....
    {'AT','CG','GC','TA'}};

reverseMap = containers.Map(k,v);

parfor i=1:Npoint
    currSeq = maxColors(i,:);
    currStr = '';
    possible = reverseMap(currSeq(1));
    for j=1:4
        p = possible{j};
        if strcmp(p(1), startBase)
            currStr = p;
        end
    end
    for kk=2:numel(currSeq)
        possible = reverseMap(currSeq(kk));
        for j=1:4
            p = possible{j};
            if strcmp(p(1), currStr(kk))
                currStr = [currStr p(2)];
            end
        end
    end
    bases{i} = currStr;
    
end
%{
function p = PossibleBases(c)
switch c
    case 1
        p = {'AA', 'CC', 'GG', 'TT'};
    case 2
        p = {'AC', 'CA', 'GT', 'TG'};
    case 3
        p = {'AG', 'CT', 'GA', 'TC'};
    case 4
        p = {'AT', 'CG', 'GC', 'TA'};
end
%}