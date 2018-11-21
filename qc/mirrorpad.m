function out = mirrorpad(in, padsize)

%MIRRORPAD - Pad a vector by mirroring edge elements.
%
% Syntax:  [out] = MIRRORPAD(in, padsize)
% 
% Pads vector with 'padsize' of values at the beginning and end of the 
% 'in' vector. The added values are mirrors of 'in'. 
%
% Inputs:
%    in - a vector
%
%    padsize - Length of the padding on each side of the vector.
%         Must be <= length(in).
%
% Outputs:
%    out - Padded in vector of length length(in)+2*padSize
%
% Example: 
%    out = mirrorpad([1:10], 3)
%
% See also: shiftarray, padseries.
%
% Author: RBR Ltd. Ottawa ON, Canada
% email: support@rbr-global.com
% Website: www.rbr-global.com
% Last revision: 2017-06-28

if isrow(in)
    pre = fliplr(in(1:padsize));
    post = fliplr(in(end-padsize+1:end));
    out = [pre in post];
elseif iscolumn(in)
    pre = flipud(in(1:padsize));
    post = flipud(in(end-padsize+1:end));
    out = [pre; in; post];
else
    error('in must be a row or column vector.')
end

end
