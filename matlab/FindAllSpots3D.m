function points = FindAllSpots3D( imgs, varargin )
%FINDALLSPOTS3D Summary of this function goes here
%  imgs is cell array of each sequencing round
%  where each sequencing round is X x Y x Z x nCh

p = inputParser();
p.addParameter('thresholdRel', 0.2, @isnumeric);
p.addParameter('ksize', 5, @isnumeric);
p.addParameter('kdims', [10 10], @isnumeric);
p.parse(varargin{:});

thresholdRel = p.Results.thresholdRel;
ksize = p.Results.ksize;
kdims = p.Results.kdims;

points = struct();
points.centroids = [];
points.seqRound = [];
points.colorCh = [];

nrounds = numel(imgs);
tic;
for n=1:nrounds
    fprintf('%s\n', repmat('-',1,20));
    fprintf('Processing sequencing round %d [time=%03f]\n', n, toc);
    currImg = imgs{n};
    [X,Y,Z,nCh] = size(currImg);
    for c=1:nCh
        currColor = currImg(:,:,:,c);
        fprintf(' *> Processing color channel %d [time=%03f]\n', c, toc);
        currCentroids = FindSpots3D(currColor, ...
            'thresholdRel', thresholdRel, ...
            'ksize', ksize, 'kdims', kdims);
        nspots = size(currCentroids,1);
        round = repmat(n, 1, nspots);
        color = repmat(c, 1, nspots);
        points.centroids = [points.centroids; currCentroids];
        points.seqRound = [points.seqRound round];        
        points.colorCh = [points.colorCh color];
    end       
end
for i=1:nrounds
    fprintf(1, '%d spots total for round %d\n', sum(points.seqRound==i), i);
end

end

