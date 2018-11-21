function out = runmed(in, windowLength, edgepad)
%RUNMED - Smooth a time series using a running median filter.
%
% Syntax:  [out] = RUNMED(in, windowLength, edgepad)
%
% Performs a running median of length windowLength. The time series is
% mirror padded, nan padded or zero-order hold padded
%
% Inputs:
%    in - Time series
%
%    windowLength - Length of the running median. It must be odd.
%
%    edgepad - Describes how the filter will act at the edges. Options
%         are 'mirror', 'zeroorderhold' and 'nan'. Default is 'mirror'.
%
% Outputs:
%    out - Smoothed time series
%
% See also: RSKsmooth, RSKdespike.
%
% Author: RBR Ltd. Ottawa ON, Canada
% email: support@rbr-global.com
% Website: www.rbr-global.com
% Last revision: 2017-06-21

  if nargin == 2
      edgepad = 'mirror';
  end

  if mod(windowLength, 2) == 0
      error('windowLength must be odd');
  end

  padsize = (windowLength-1)/2;
  inpadded = padseries(in, padsize, edgepad);

  n = length(in);
  out = NaN*in;
  for ndx = 1:n
      out(ndx) = nanmedian(inpadded(ndx:ndx+(windowLength-1)));
  end

end
