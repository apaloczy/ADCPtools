function inpadded = padseries(in, padsize, edgepad)

%PADSERIES - Add padsize amount of values of either side of the in vector.
%
% Syntax:  [inpadded] = PADSERIES(in, padsize, edgepad)
%
% Inputs:
%    in - Vector
%
%    padsize - Amount of values added to either edge of the vector.
%
%    edgepad - Describes what values will be added at the edges. Options
%         are 'mirror', 'zeroorderhold' and 'nan'.
%
% Outputs:
%    inpadded - Series with length 2*padsize+length(in).
%
% See also: runavg, runmed, runtriang.
%
% Author: RBR Ltd. Ottawa ON, Canada
% email: support@rbr-global.com
% Website: www.rbr-global.com
% Last revision: 2017-06-20

switch edgepad
    case 'mirror'
        inpadded = mirrorpad(in, padsize);
    case 'nan'
        inpadded = nanpad(in, padsize);
    case 'zeroorderhold'
        inpadded = zeroorderholdpad(in, padsize);
    otherwise
        error('edgepad argument is not recognized. Must be `mirror`, `nan` or `zeroorderhold`');
end

end