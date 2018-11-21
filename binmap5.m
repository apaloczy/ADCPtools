function [b1m, b2m, b3m, b4m, b5m] = binmap5(b1, b2, b3, b4, b5, r, r5, theta, ptch, roll, how)
%USAGE
%-----
%[b1m, b2m, b3m, b4m, b5m] = binmap5(b1, b2, b3, b4, b5, r, r5, theta, ptch, roll, how)
%
%Interpolate beam-coordinate velocities to fixed horizontal planes based on tilt angles
%(pitch and roll).
d2r = pi./180;
theta = theta.*d2r;
Sth = sin(theta);
Cth = cos(theta);

ptch = ptch.*d2r;
roll = roll.*d2r;
Sph2 = sin(ptch);
Cph2 = cos(ptch);
Sph3 = sin(roll);
Cph3 = cos(roll);

Z = r.*Cth;
dZu = Z(end) - Z(end-1);
rmax = r(end);
nz = length(Z);
z00 = [0 0 1]';

nt = length(ptch);
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
          j = j - 1;
          jj = jj - 1;
        end
        zij = zi(j);
        zijj = zi(jj);
        dzj = zijj - zij;
        bmi(J,k) = ((Zj - zij)./dzj).*boi(j) + ((zijj - Zj)./dzj).*boi(jj);
      elseif strcmp(how, 'nn')                         % Nearest-neighbor interpolation.
        bmi(J,k) = boi(near(zi, Zj));
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
