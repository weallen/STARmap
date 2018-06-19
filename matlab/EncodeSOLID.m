function colorSeq = EncodeSOLID( seq )

k = {'AT','CT','GT','TT',...
    'AG', 'CG', 'GG', 'TG',...
    'AC', 'CC', 'GC', 'TC',...
    'AA', 'CA', 'GA', 'TA'};
v = {4,3,2,1,3,4,1,2,2,1,4,3,1,2,3,4};
coding = containers.Map(k,v);
colorSeq = [coding(seq(1:2))];
for i=3:numel(seq)
    s = seq((i-1):i);
    colorSeq = [colorSeq coding(s)];
end
end

