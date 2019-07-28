function [ TKEO, MAD3, SegFilter_DEstimate, SegFilter_AEstimate, CPs, fco ] = segfilt_estim(RawDisplace)
% 
% Function implements a segment filtering procedure to process
% non-stationary signals. Segments the signal at change points defined by a
% median absolute deviation criterion. Applies an automatically determined
% filter to each section using the autocorrelation-based procedure (ABP;
% Challis, 1999).

% Details can be found in, "Davis DJ, Challis JH (2019). A filtering
% procedure to process non-stationary signals. Proceedings of the XXVII's
% Congress of the International Society of Biomechanics."

%-----------
%
% Daniel J. Davis, The Pennsylvania State University, July 2019
% John H. Challis, The Pennsylvania State University
%
% -----------
%
% INPUTS
% RawDisplace - noisy displacement data to be filtered and differentiated
% TKEOPlot    - option to plot TKEO
% AccelPlot   - option to plot the acceleration estimates
%
% OUTPUTS
% TKEO        - Teager Kaiser Energy Operator of signal
% MAD3        - 3 x median absolute deviation of TKEO, outlier criterion
% SegFilter_DEstimate - Displacement estimate with segmental filtering
% SegFilter_AEstimate = Acceleration estimate with segmental filtering
%
%
% NOTES
%  -----
%  ABP from: John H. Challis, The Penn. State University (December 4, 1997)
%
%  The Teager Energy Operator is determined using,
%  x(n) = x^2(n) + x(n - 1)x(n + 1)
%
%  MAD is determined using,
%  MAD = median ( abs( xi - Mn))
%       where b = normality constant = 1.4826
%       xi = given data point
%       Mn = median of set of data points

global forder dt

% single-filter using autocorrelation-based procedure
[ABP_DisplEstim, ~] = filtmat_auto(dt,forder,RawDisplace);
% [ABP_DisplEstim] = filtmat(dt, 2, forder, RawDisplace);


% filter raw displacement data first as Teager-Kaiser Energy Operator (TKEO) 
% is sensitive to noise. Cut-off as selected by ABP on entire signal
[ TeagEn ] = TeagerEnergy( ABP_DisplEstim );

% ABP is computed on TKEO to determine cut-off frequency as may still
% contain noise
[ TKEO, ~ ] = filtmat_auto( dt, forder, TeagEn);

% Compute 3* Median Absolute Deviation 
MAD3  = 3*mad(TKEO);

% Determine number of change-points and each sample index (CPs)
[numCPs, CPs] = changePoints(TKEO, MAD3);
numSegs = numCPs+1;

% initialize
fco = zeros(1,numSegs);

%  First Segment  %
% determines cut-off frequency based on raw displacement data in segment 1
% and filters raw data with said cut-off frequency
[ ~, fco(1) ] = filtmat_auto( dt, forder, RawDisplace(1:CPs(1)));

%  use frequency cut-off to double-differentiate the data, creating an acceleration
%  estimate
seg_Destims(:,1) =  filtmat(dt, fco(1), forder, RawDisplace); 
seg_Aestims(:,1) =  fdiff2(dt, fco(1), forder, RawDisplace);                 

% if procedure determines multiple frequency cut-offs are necessary...
if length(fco) > 1
    
%  Middle Segment(s)  %
for p = 2:numCPs
    
[ ~, fco(p) ] = filtmat_auto( dt, forder, RawDisplace(CPs(p-1)+1:CPs(p)));
seg_Destims(:,p) =  filtmat(dt, fco(p), forder, RawDisplace); 
seg_Aestims(:,p) =  fdiff2(dt, fco(p), forder, RawDisplace);         

end

%  Last segment %
[ ~, fco(numSegs) ] = filtmat_auto( dt, forder, RawDisplace(CPs(end)+1:end));
seg_Destims(:,numSegs) =  filtmat(dt, fco(numSegs), forder, RawDisplace); 
seg_Aestims(:,numSegs) =  fdiff2(dt, fco(numSegs), forder, RawDisplace);         


%  Join segment estimates at the change-points using weighted-average joins. "pad" is
% simply a number of samples on either side of the change-points used for
% smoothing
[ SegFilter_DEstimate ] = join(CPs, seg_Destims);
[ SegFilter_AEstimate ] = join(CPs, seg_Aestims);

else   % if the number of change points = 1 (ie, no change point necessary)
    
    SegFilter_DEstimate = seg_Destims;
    SegFilter_AEstimate = seg_Aestims;
    
end

end


%% THE END %%