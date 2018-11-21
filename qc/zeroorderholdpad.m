function out = zeroorderholdpad(in, padsize)

%ZEROORDERHOLDPAD - Pad a vector by entering zero-order hold values.
%
% Syntax:  [out] = ZEROORDERHOLDPAD(in, padsize)
% 
% Pads vector with 'padsize' of values at the beginning and end of the 'in'
% vector. The added values are repeats of the first and last entry. 
%
% Inputs:
%    in - Vector
%
%    padsize - Length of the padding on each side of the vector.
%         Must be <= length(in).
%
% Outputs:
%    out - Padded in vector of length length(in)+2*padsize.
%
% Example: 
%    out = zeroorderholdpad([1:10], 3)
% 
%    out =
%
%    [1 1 1 1 2 3 4 5 6 7 8 9 10 10 10 10]
%
% Author: RBR Ltd. Ottawa ON, Canada
% email: support@rbr-global.com
% Website: www.rbr-global.com
% Last revision: 2017-06-21

if isrow(in)
    pre = repmat(in(1), [1,padsize]);
    post = repmat(in(end), [1,padsize]);
    out = [pre in post];
elseif iscolumn(in)
    pre = repmat(in(1), [padsize,1]);
    post = repmat(in(end), [padsize,1]);
    out = [pre; in; post];
else
    error('in must be a row or column vector.')
end

end
