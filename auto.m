function [ coopt ]= auto( nsignal, dt, forder)

%  purpose - procedure to work out optimal cut-off frequency using
%  -------   autocorrelation of residual of signal using iterative procedure
%
%  John H. Challis, The Penn. State University (December 4, 1997)
%
%  calling - [ coopt ]= auto( nsignal, dt, forder);
%  -------
%
%  input
%  -----
%  nsignal - array containing noisy signal
%  dt      - interval between samples
%  forder  - order of filter
%
%  output
%  ------
%  coopt   - optimum cut-off
%
%  calls
%  -----
%  filtmat - filters a matrix of data
%
%  notes
%  -----
%  1)  These computations are based on the algorithm described in
%      Challis, J.H. (1999)  A procedure for the automatic determination of filter cutoff
%      frequency for the processing of biomechanical data.  Journal of Applied Biomechanics  15:303-317.
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              %
%  set upper and lower bounds  %
%                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colow = 0.5;	       % lower bounds
coup  = 0.25 / dt;	   % upper bounds


%%%%%%%%%%%%%%%%%%%%%
%                   %
%  starting values  %
%                   %
%%%%%%%%%%%%%%%%%%%%%
ssignal = filtmat(dt, colow, forder, nsignal);  	        %  filter signal
resid = ssignal - nsignal;  								%  compute the residual
acorr = xcorr(resid,'coeff');								%  compute autocorrelation of residual of signal
fmin = sum(abs(acorr));										%  objective function
coopt = colow;												%  initial guess of optimum cut-off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                          %
%  loop with variable cut-off frequencies  %
%                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for co = colow+0.1:0.1:coup
%
    ssignal = filtmat(dt, co, forder, nsignal);                 %  filter signal
    resid = ssignal - nsignal;									%  compute the residual
    [acorr] = xcorr(resid,'coeff');								%  compute autocorrelation of residual of signal
    f = sum(abs(acorr));										%  objective function

%  is this a minimum?
%
    if f < fmin
        fmin = f;
    	coopt = co;
    end

%  loop back for next cut-off
%
end


%
%%
%%% The End %%%
%%
%
