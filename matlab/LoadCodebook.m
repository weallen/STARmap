function [ nameToSeq, seqToName ] = LoadCodebook( dirname, doFlip)
if nargin < 2
    doFlip = true;
end
genes = {};
seqs = {};
fname = fullfile(dirname, 'genes.csv');
f = fopen(fname,'r');
tline = fgets(f);
fields = strsplit(tline,',');

genes{end+1} = fields{1}; 
if doFlip
    seqs{end+1} = fliplr(strtrim(fields{2}));
else
    seqs{end+1} = strtrim(fields{2});
end
while ischar(tline)
    tline = fgets(f);    
    if tline ~= -1
        fields = strsplit(tline,',');
        if isempty(fields{1}) || isempty(fields{2})
            break;
        else
            genes{end+1} = fields{1};
            if doFlip
                seqs{end+1} = fliplr(strtrim(fields{2}));
            else
                seqs{end+1} = strtrim(fields{2});
            end
        end
    end    
end

fclose(f);

for i=1:numel(seqs)
    seqs{i} = Colorseq2Str(EncodeSOLID(seqs{i}));
end
seqToName = containers.Map(seqs, genes);
nameToSeq = containers.Map(genes, seqs);

