function [b1, b2, b3, b4] = janus3beamsol(b1, b2, b3, b4)
% Usage
% -----
% [b1, b2, b3, b4] = janus3beamsol(b1, b2, b3, b4)
%
% Calculates a three-beam solution for a bad beam when the other three Janus beams have good data.
[Nz Nt] = size(b1);

for nt=1:Nt
  for nz=1:Nz            % Set error velocity to zero: const*(b1 + b2 -b3 -b4) = 0.
    bki = [b1(nz,nt) b2(nz,nt) b3(nz,nt) b4(nz,nt)];
    fbad = isnan(bki);   % b1 + b2 = b3 + b4; Solve for bad beam.
    if sum(fbad)==1      % Only one bad beam allowed for 3-beam solutions.
      fbad = find(fbad);
      if fbad==1     % Beam 1 is bad.
        b1(nz,nt) = bki(3) + bki(4) - bki(2);
      elseif fbad==2 % Beam 2 is bad.
        b2(nz,nt) = bki(3) + bki(4) - bki(1);
      elseif fbad==3 % Beam 3 is bad.
        b3(nz,nt) = bki(1) + bki(2) - bki(4);
      elseif fbad==4 % Beam 4 is bad.
        b4(nz,nt) = bki(1) + bki(2) - bki(3);
      end
    end
end

end
