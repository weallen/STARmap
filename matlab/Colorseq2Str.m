function s = Colorseq2Str( cs )
s = '';
for i=1:numel(cs)
    s = [s num2str(cs(i))];
end

end

