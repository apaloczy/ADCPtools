function [udespiked, ispikes] = ztdespike(u, t, varargin)
% USAGE
% -----
% [udespiked, ispikes] = ztdespike(u, t, varargin)
%
% Identifies and treats spikes in a time series of vertical profiles of a variable 'u(z,t)'
% using a median filtering algorithm. A reference time series is created by filtering the input
% variable 'u' with a median filter of length 'windowLength'. A residual ("high-pass") series is
% formed by subtracting the reference series from the original signal. Data in the reference
% series lying outside of 'threshold' standard deviations are defined as spikes. Spikes are
% then treated by one of three methods (see below).
%
% Inputs:
% -------
%   [Required]   u - [nz x nt] matrix containing a time series of vertical
%                    profiles of a variable u(k, t), where 'k' is the
%                    vertical index.
%
%                t - [nt] Time vector associated with u.
%   [Optional]
%
%                threshold - Amount of standard deviations to use for the
%                      spike criterion. Default value is 2.
%
%                windowLength - Total size of the filter window. Must be
%                      odd. Default is 3.
%
%                action - Action to perform on a spike. The default is 'nan',
%                      whereby spikes are replaced with NaN.  Other options are
%                      'replace', whereby spikes are replaced with the
%                      corresponding reference value, and 'interp',
%                      whereby spikes are replaced with values calculated
%                      by linearly interpolating from the neighbouring
%                      points.
%
%                visualize - To give a diagnostic plot on specified profile
%                          number(s). Default is 0. If nonzero, it specifies
%                          the number of seconds to pause between bins.
%
% Ouputs:
% -------
% udespiked   - despiked time series.
% ispikes     - indices of the data points flagged as spikes.
%
% This code is a minor modification of the RSKdespike.m function included in the
% RSKtools package, developed by RBR Ltd. (https://bitbucket.org/rbr/rsktools).

  validActions = {'replace', 'interp', 'nan'};
  checkAction = @(u) any(validatestring(u,validActions));

  p = inputParser;
  addRequired(p, 'u', @isnumeric);
  addRequired(p, 't', @isnumeric);
  addOptional(p, 'threshold', 2, @isnumeric);
  addOptional(p, 'windowLength', 3, @isnumeric);
  addOptional(p, 'action', 'nan', checkAction);
  addOptional(p, 'visualize', 0, @isnumeric);
  parse(p, u, t, varargin{:})

  threshold = p.Results.threshold;
  windowLength = p.Results.windowLength;
  action = p.Results.action;
  visualize = p.Results.visualize;
  [nz, nt] = size(u);

  if visualize~=0
    f = figure; hold on; grid on;
  end

  udespiked = u.*nan;
  it = 1:nt; ispikes = [];
  for k=1:nz
      uk = u(k, :);
      [udspk, ispk] = despike(uk, t, threshold, windowLength, action);
      udespiked(k, :) = udspk;
      ispikes = [ispikes; ispk];
      if visualize~=0
        plot(it, uk, 'r');
        plot(it, udspk, 'k');
        scatter(it, uk(ispk), 'r*');
        xlabel('Time index'); ylabel('u');
        title('red=original, black=despiked, red stars=spike locations');
        pause(visualize)
      end
  end

  if visualize~=0
    close(f)
  end

  function [y, I] = despike(x, t, threshold, windowLength, action)
  % Replaces the values that are > threshold*standard deviation away from
  % the residual between the original time series and the running median
  % with the median, a NaN or interpolated value using the non-spike
  % values. The output is the x series with spikes fixed and I is the
  % index of the spikes.
  %
  % Source: 'RSKdespike.m' code in the RSKtools package (https://bitbucket.org/rbr/rsktools).

  y = x;
  ref = runmed(x, windowLength);
  dx = x - ref;
  sd = std(dx(isfinite(dx)));
  I = abs(dx) > threshold*sd;
  good = ~I;

  switch action
    case 'replace'
      y(I) = ref(I);
    case 'nan'
      y(I) = NaN;
    case 'interp'
      try
        y(I) = interp1(t(good), x(good), t(I));
      catch
        y(I) = NaN;
        warning(['Not enough good bins to interpolate over bad bins (' datestr(t(1)) ' - ' datestr(t(end)) ').']);
      end
  end
  end
end
