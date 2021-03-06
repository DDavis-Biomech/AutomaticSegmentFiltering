function [SegFilter_Estimate] = join_sections(CPs, segment_estims)
% -----------
%
% Daniel J. Davis, The Pennsylvania State University, July 2019
%
% -----------
%
% Creates a signal estimate by joining together multiple segments which 
% have been filtered differently based on their frequency profiles. Uses 
% moving average function to smooth signal discrepancies at signal
% joins due to the use of different filter cut-off frequencies
%
% ----------
%
% INPUTS
% CPs        - list of sample indexes of change points
% segment_estims - segment estimates created by filtering signal at cut-off
%   frequencies defined by respective segment

[r,~] = size(CPs);
if r > 1
    CPs = CPs';
end

% number of change points and segments
numCPs = length(CPs);
numSegs = numCPs + 1;

% find shortest segment length to determine amount of "pad" on each side of
% join
evalCPs = [0, CPs, length(segment_estims)];
if numCPs > 1
shortestSmooth = min(diff(evalCPs));

% round the pad to be an integer
mult = 2;
roundDown = mult * floor((shortestSmooth-1)/mult) + 1;
roundUp = mult * ceil((shortestSmooth-1)/mult) + 1;

% potential residual after rounding
roundResi = [abs(roundDown - shortestSmooth), abs(roundUp - shortestSmooth)];

UporDown = find(max(roundResi));

if UporDown == 1
       smoothLength = roundDown;
else
       smoothLength = roundUp;
end

% if only 1 CP, smoothLength = 11. Number can be modified based on sampling
% frequency and data (IS NOT USED IN DOWLING DATA or in ASFP manuscript)
else 
    smoothLength = 11;
end


%  pads on each side of change point
pad = (smoothLength-1)/2;

%  row vector with weights to be applied to differences 
weights = linspace(1,0,smoothLength)';


% create joins where sections meet

%preallocate
overlaps = zeros(length(smoothLength), numCPs);
difs = zeros(length(smoothLength), numCPs);
weightedDifs = zeros(length(smoothLength), numCPs);
joins = zeros(length(smoothLength), numCPs);


for j = 1:numCPs   

% creates two overlaps per change point (one for each input into
% weighted-average function)
overlaps(1:smoothLength,(j*2)-1) = segment_estims(CPs(j) - pad : CPs(j) + pad, j);
overlaps(1:smoothLength,j*2) = segment_estims(CPs(j) - pad : CPs(j) + pad, j+1);


%  produce smoothed section "joins" from overlaps on each side of change
%  points
for r = 1:smoothLength
    
    difs(r,j) = overlaps(r,(j*2)-1) - overlaps(r,j*2);
    
    weightedDifs(r,j) = weights(r) * difs(r,j);
    
    joins(r,j) = overlaps(r,j*2) + weightedDifs(r,j);
    
end

end

%% Join segments %%

% first segment
% index from beginning of data up to the beginning of the pad (which is
% found at the change points - the pad +1)
first2pad_indx = 1:(CPs(1)-(pad+1));
SegFilter_Estimate(1,first2pad_indx) = segment_estims(first2pad_indx,1)';

for p = 2:numCPs
    
% index from beginning to ending of pad region
pad_indx = (CPs(p-1) - (pad+1))+1:(CPs(p-1)+pad);

% index from end of pad to beginning of next pad
pad2pad_indx = (CPs(p-1)+pad+1) : CPs(p) - (pad + 1);

SegFilter_Estimate(1,pad_indx) = joins(:,p-1)';
    
SegFilter_Estimate(1,pad2pad_indx) = segment_estims(pad2pad_indx,p)';

end

% last segment
% index from begin to end of last pad
lastpad_indx = (CPs(end) - (pad+1))+1:(CPs(end)+pad);
% index from last pad to end of data
pad2end_indx = (CPs(end)+pad+1):length(segment_estims(:,1));

SegFilter_Estimate(1,lastpad_indx) = joins(:,end)';
SegFilter_Estimate(1,pad2end_indx) = segment_estims(pad2end_indx,numSegs)'; 


% transpose to match input data
SegFilter_Estimate = SegFilter_Estimate';


end

%% The End %%

