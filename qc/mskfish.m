function [b1, b2, b3, b4] = mskfish(b1, b2, b3, b4, amp1, amp2, amp3, amp4, varargin)
% USAGE
% -----
% [b1 b2 b3 b4] = mskfish(b1, b2, b3, b4, amp1, amp2, amp3, amp4)
%
%                 ---OR---
%
% [b1 b2 b3 b4] = mskfish(b1, b2, b3, b4, amp1, amp2, amp3, amp4,
%                        'Threshold', thresh, 'use3BeamSol', true|false)
%
% REFERENCE
% ---------
% ADCP Coordinate transformation: Formulas and calculations (2010), p 22-23
% P/N 951-6079-00 (January 2010), Teledyne RD Instruments
% Available at:
% http://www.teledynemarine.com/Documents/Brand%20Support/RD%20INSTRUMENTS/
% Technical%20Resources/Manuals%20and%20Guides/General%20Interest/Coordinate_Transformation.pdf
options = struct('Threshold', 50, 'use3BeamSol', false);
optionNames = fieldnames(options);

if any(strcmp(varargin, 'Threshold'))
  thresh = varargin{find(strcmp(varargin, 'Threshold'))+1};
else
  thresh = options.Threshold; % Default value.
end

if any(strcmp(varargin, 'use3BeamSol'))
  use3BeamSol = varargin{find(strcmp(varargin, 'use3BeamSol'))+1};
else
  use3BeamSol = options.use3BeamSol; % Default value.
end

[nz, nt] = size(b1);
for i=1:nt
  Bi = [b1(:,i) b2(:,i) b3(:,i) b4(:,i)];
  for k=1:nz
    ampk = [amp1(k,i) amp2(k,i) amp3(k,i) amp4(k,i)];
    ampmax = max(ampk);
    [ampmin fm] = mink(ampk, 2); % Weakest and second-weakest echo.
    damp = ampmax - ampmin;

    if damp(1)>thresh  % Fish in at least 1 beam.
      Bi(fm(1)) = NaN; % Can still calculate 3-beam solutions.
      if use3BeamSol==true
        Bi = sol3b(Bi);
      end
      if k<nz                % Also mark cells k+1 as bad, because
        Bi(k+1,fm(1)) = NaN; % echo is measured at the end of the cells.
      end

      if damp(2)>thresh % Fish in at least 2 beams.
        Bi(k,:) = NaN;  % Mark all beams as bad.
        if k<nz
          Bi(k+1,:) = NaN;
        end
      end
    else % No fish detected on any beam at cell k.
      ;
    end
  end
  b1(:,i) = Bi(:,1);
  b2(:,i) = Bi(:,2);
  b3(:,i) = Bi(:,3);
  b4(:,i) = Bi(:,4);
end

end

function bki = sol3b(bki)
  fbad = isnan(bki);   % b1 + b2 = b3 + b4; Solve for bad beam.
  if sum(fbad)==1      % Only one bad beam allowed for 3-beam solutions.
    fbad = find(fbad);
    if fbad==1     % Beam 1 is bad.
      bki(1) = bki(3) + bki(4) - bki(2);
    elseif fbad==2 % Beam 2 is bad.
      bki(2) = bki(3) + bki(4) - bki(1);
    elseif fbad==3 % Beam 3 is bad.
      bki(3) = bki(1) + bki(2) - bki(4);
    elseif fbad==4 % Beam 4 is bad.
      bki(4) = bki(1) + bki(2) - bki(3);
    end
  end

end
