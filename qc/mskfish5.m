function [b1, b2, b3, b4, b5] = mskfish5(b1, b2, b3, b4, b5, amp1, amp2, amp3, amp4, amp5, varargin)
% USAGE
% -----
% [b1 b2 b3 b4 b5] = mskfish5(b1, b2, b3, b4, b5, amp1, amp2, amp3, amp4, amp5, 'Threshold', thresh)
%
% REFERENCE
% ---------
% Adapted from:
% ADCP Coordinate transformation: Formulas and calculations (2010), p 22-23
% P/N 951-6079-00 (January 2010), Teledyne RD Instruments
% Available at:
% http://www.teledynemarine.com/Documents/Brand%20Support/RD%20INSTRUMENTS/
% Technical%20Resources/Manuals%20and%20Guides/General%20Interest/Coordinate_Transformation.pdf
options = struct('Threshold', 50);
optionNames = fieldnames(options);

if any(strcmp(varargin, 'Threshold'))
  thresh = varargin{find(strcmp(varargin, 'Threshold'))+1};
else
  thresh = options.Threshold; % Default value.
end

[nz, nt] = size(b1);
for i=1:nt
  Bi = [b1(:,i) b2(:,i) b3(:,i) b4(:,i) b5(:,i)];
  for k=1:nz
    ampk = [amp1(k,i) amp2(k,i) amp3(k,i) amp4(k,i), amp5(k,i)];
    ampmax = max(ampk);
    [ampmin fm] = mink(ampk, 3); % Weakest, second-weakest and third-weakest echo.
    damp = ampmax - ampmin;

    if damp(1)>thresh  % Fish in at least 1 beam.
      Bi(fm(1)) = NaN;
      if k<nz                % Also mark cells k+1 as bad, because
        Bi(k+1,fm(1)) = NaN; % echo is measured at the end of the cells.
      end

        if damp(2)>thresh % Fish in at least 2 beams.
          Bi(fm(2)) = NaN;
          if k<nz                % Also mark cells k+1 as bad, because
            Bi(k+1,fm(2)) = NaN; % echo is measured at the end of the cells.
          end

          if damp(3)>thresh % Fish in at least 3 beams.
            Bi(k,:) = NaN;  % Mark all beams as bad.
            if k<nz
              Bi(k+1,:) = NaN;
            end
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
  b5(:,i) = Bi(:,5);
end
