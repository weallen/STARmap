function cs = Str2Colorseq( s )
%STR2COLORSEQ Summary of this function goes here
%   Detailed explanation goes here
cs = [];
for i=1:numel(s)
    cs = [cs str2num(s(i))];
end

end

