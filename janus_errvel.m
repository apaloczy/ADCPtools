function Verr = janus_errvel(b1, b2, b3, b4, theta)
% Usage
% -----
% Verr = janus_errvel(b1, b2, b3, b4, theta)
%
% Calculates the error velocity, defined as the scaled difference between
% the two independent measurements of vertical velocity (w_12 and w_34),
% Verr = (b1 + b2 -b3 -b4)./(2.*sin(theta))./sqrt(2).
[Nz Nt] = size(b1);

B = cat(3, b1, b2, b3, b4);

%     b1 b2 b3 b4
A = [+1 +1 -1 -1];

Verr = zeros(Nz, Nt);
for nz=1:Nz
  disp(['Calculating error velocity at bin ', num2str(nz), '/', num2str(Nz)])
  for nt=1:Nt
    Verr(nz,nt) = A*squeeze(B(nz,nt,:));
  end
end

d2r = pi/180;
theta = theta.*d2r;

% Following notation on p. 11 of TRDI's ADCP coordinate tramsformation
% primer (ADCP Coordinate Transformation - formulas and calculations,
% P/N 951-6079-00, January 2008);
a = 1./(2.*sin(theta));
d = a./sqrt(2);
Verr = Verr.*d;

end
