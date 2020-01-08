function [ TKEO, MAD3, SegFilter_DisplEstim, SegFilter_AccelEstim, CPs, fco ] = ASFP(dt, forder, RawDisplace)
% -----------
%
% Daniel J. Davis and John H. Challis 
%   The Pennsylvania State University, July 2019
%
% -----------
%
% Function implements the Automatic Segment Filtering Procedure to process
% non-stationary signals. Segments the signal at change points defined by a
% median absolute deviation criterion. Applies an automatically determined
% filter to each section using the autocorrelation-based procedure (ABP;
% Challis, 1999).

% Details can be found in, "Davis D.J., Challis J.H. (in Press). Automatic segment 
%   filtering procedure for processing non-stationary signals. Journal of
%   Biomechanics."
%
% -------------
%
% INPUTS
% RawDisplace - noisy displacement data to be filtered and differentiated
%
% OUTPUTS
% TKEO        - Teager Kaiser Energy Operator of signal
% MAD3        - 3 x median absolute deviation of absolute value of TKEO, outlier criterion
% SegFilter_DisplEstim - Displacement estimate with segmental filtering
% SegFilter_AccelEstim - Acceleration estimate with segmental filtering
% CPs         - Sample indexes of change points as defined by the 3xMAD criterion
% fco        - Filter cut-off frequencies for each segment 
%
%
% NOTES
%  -----
%
%  The Teager Energy Operator is determined using,
%  x(n) = x^2(n) + x(n - 1)x(n + 1)
%
%  MAD is determined using,
%  MAD = median ( abs( xi - Mn))
%       xi = given data point
%       Mn = median of set of data points
%
% ------------------------------------------------------------------------

% single-filter using autocorrelation-based procedure
[ABP_DisplEstim, coABP] = filtmat_auto(dt,forder,RawDisplace);

% filter raw displacement data first as Teager-Kaiser Energy Operator (TKEO) 
% is sensitive to noise. Cut-off as selected by ABP on entire signal
[ TKEO ] = TeagerEnergy( ABP_DisplEstim ); 


% Compute 3* Median Absolute Deviation of absolute value of TKEO 
MAD3  = 3*mad(abs(TKEO));

% Determine number of change-points and each sample index (CPs)
[numCPs, CPs] = changePoints(TKEO, MAD3);
numSegs = numCPs+1;

% if only one segment identified, the single-filter approach is appropriate
if numSegs == 1
    SegFilter_DisplEstim = filtmat(dt, coABP, 2, RawDisplace);
    SegFilter_AccelEstim = fdiff2(dt, coABP, 2, RawDisplace);
    fco = coABP;
else
    
% initialize
fco = zeros(1,numSegs);

%  First Segment  %
% determines cut-off frequency based on raw displacement data in segment 1
% and filters raw data with said cut-off frequency
[ fco(1) ] = auto( RawDisplace(1:CPs(1)), dt, forder );

%  use frequency cut-off to double-differentiate the data, creating an acceleration
%  estimate for the full signal
segment_DisplEstims(:,1) =  filtmat(dt, fco(1), forder, RawDisplace); 
segment_AccelEstims(:,1) =  fdiff2(dt, fco(1), forder, RawDisplace);                 


% if procedure determines if multiple frequency cut-offs are necessary...
if length(fco) > 1

%  Middle Segment(s)  %
% creates displacement and acceleration estimates for each segment-specific
% cut-off frequency
for p = 2:numCPs
    
    [ fco(p) ] = auto( RawDisplace(CPs(p-1)+1:CPs(p)), dt, forder );
    segment_DisplEstims(:,p) =  filtmat(dt, fco(p), forder, RawDisplace); 
    segment_AccelEstims(:,p) =  fdiff2(dt, fco(p), forder, RawDisplace);         

end

%  Last segment %
[ fco(numSegs) ] = auto( RawDisplace(CPs(end)+1:end), dt, forder );
segment_DisplEstims(:,numSegs) =  filtmat(dt, fco(numSegs), forder, RawDisplace); 
segment_AccelEstims(:,numSegs) =  fdiff2(dt, fco(numSegs), forder, RawDisplace);         


%  Join segment estimates at the change-points using weighted-average joins. "pad" is
% simply a number of samples on either side of the change-points used for
% smoothing
[ SegFilter_DisplEstim ] = join(CPs, segment_DisplEstims);
[ SegFilter_AccelEstim ] = join(CPs, segment_AccelEstims);


elseif isemtpy(fco) == 1  % if the number of change points = 1 (ie, no change point necessary)
    
    SegFilter_DisplEstim = segment_DisplEstims;
    SegFilter_AccelEstim = segment_AccelEstims;
    
   
end

end


%% THE END %%
