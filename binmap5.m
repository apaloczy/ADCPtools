function [b1m, b2m, b3m, b4m, b5m] = binmap5(b1, b2, b3, b4, b5, r, r5, theta, ptch, roll, how)
%USAGE
%-----
%[b1m, b2m, b3m, b4m, b5m] = binmap5(b1, b2, b3, b4, b5, r, r5, theta, ptch, roll, how)
%
% theta, ptch and roll in RADIANS.
%
%Interpolate beam-coordinate velocities to fixed horizontal planes based on tilt angles
%(pitch and roll).
Sth = sin(theta);
Cth = cos(theta);

Sph2 = sin(ptch);
Cph2 = cos(ptch);
Sph3 = sin(roll);
Cph3 = cos(roll);

Z = r.*Cth;
Zbot = Z(1);
Ztop = Z(end);
dZu = Z(end) - Z(end-1);
rmax = r(end);
z00 = [0 0 1]';

[nz, nt] = size(b1);
for k=1:nt
  PRk = [Cph3(k)             0       Sph3(k);
         Sph2(k).*Sph3(k)  Cph2(k)  -Sph2(k).*Cph3(k);
        -Sph3(k).*Cph2(k)  Sph2(k)   Cph2(k).*Cph3(k)];
  PR(:,:,k) = PRk;
end

%      b1    b2    b3    b4   b5
E = [-Sth  +Sth    0     0    0;
       0     0   -Sth  +Sth   0;
     -Cth  -Cth  -Cth  -Cth  -1];

Bo = cat(3, b1, b2, b3, b4);

for i=1:5
  Ei = E(:,i);

  if i==5
    Boi = b5;
    r = r5;
  else
    Boi = Bo(:,:,i); % z, t, bi.
  end

  for k=1:nt
    zi = abs((PR(:,:,k)*Ei)'*z00).*r; % Actual bin height, dot product of tilt matrix with along-beam distance vector.
    boi = Boi(:,k);

    for J=1:nz
      Zj = Z(J);
      if strcmp(how, 'linear')                         % Linear interpolation.
        j = nearfl(zi, Zj);
        jj = j + 1;
        if jj>nz
          jj = nz;
          j = nz - 1;
        end
        zij = zi(j);
        zijj = zi(jj);
        if Zbot>zij
          bmi(J,k) = boi(1);
        elseif Ztop<zijj
          bmi(J,k) = boi(end);
        else
          dzj = zijj - zij;
          bmi(J,k) = ((Zj - zij)./dzj).*boi(j) + ((zijj - Zj)./dzj).*boi(jj);
        end
      elseif strcmp(how, 'nn')                         % Nearest-neighbor interpolation.
        j = near(zi, Zj);
        if j>nz
          j = nz;
        end
        bmi(J,k) = boi(j);
      end
    end
  end

  if i==5
    b5m = bmi;
  else
    Bm(:,:,i) = bmi;
  end

end

b1m = Bm(:,:,1);
b2m = Bm(:,:,2);
b3m = Bm(:,:,3);
b4m = Bm(:,:,4);

end
