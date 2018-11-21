function out = nanpad(in, padsize)

%NANPAD - Pad a vector with NaNs.
%
% Syntax:  [out] = NANPAD(in, padsize)
% 
% Pads the 'in' vector with 'padsize' of values at the beginning and end of the
% 'in' vector. The added values are NaNs.
%
% Inputs:
%    in - Vector
%
%    padSize - Length of the padding on each side of the vector.
%         Must be <= length(in).
%
% Outputs:
%    out - Padded in vector of length length(in)+2*padSize
%
% Example: 
%    out = NANPAD([1:10], 3)
% 
%    out =
%
%    [NaN NaN NaN 1 2 3 4 5 6 7 8 9 10 NaN NaN NaN]
%
% See also: padseries, shiftarray.
%
% Author: RBR Ltd. Ottawa ON, Canada
% email: support@rbr-global.com
% Website: www.rbr-global.com
% Last revision: 2017-06-20


if isrow(in)
    pad = NaN(1,padsize);
    out = [pad in pad];
elseif iscolumn(in)
    pad = NaN(padsize,1);
    out = [pad; in; pad];
else
    error('in must be a row or column vector.')
end
end
