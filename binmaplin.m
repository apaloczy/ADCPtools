function [b1m, b2m, b3m, b4m, b5m] = binmaplin(b1, b2, b3, b4, b5, r, theta, ptch, roll)
%USAGE
%-----
%[b1m, b2m, b3m, b4m, b5m] = binmaplin(b1, b2, b3, b4, b5, r, theta, ptch, roll)
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
while Z(end)<(rmax-dZu)
  Z = [Z Z(end)+dZu]; % Complete vertical axis to map to.
end
z00 = [0 0 1]';

nt = length(ptch);
for k=1:nt
  PRk = [Cph3(k)             0       Sph3(k);
         Sph2(k).*Sph3(k)  Cph2(k)  -Sph2(k).*Cph3(k);
        -Sph3(k).*Cph2(k)  Sph2(k)   Cph2(k).*Cph3(k)];
  PR(:,:,k) = PRk;
end

%     b1    b2    b3    b4   b5
E = [-Sth  +Sth    0     0    0;
       0     0   -Sth  +Sth   0;
     -Cth  -Cth  -Cth  -Cth  -1];

Bo = cat(3, b1, b2, b3, b4, b5);
for i=1:5
  Ei = E(:,i);
  Boi = Bo(:,:,i); % z, t, bi.

  for k=1:nt
    zi = abs((PR(:,:,k)*Ei)'*z00).*r; % Actual bin height, dot product of tilt matrix with along-beam distance vector.
    boi = Boi(:,k);

    for J=1:nz-1
      Zj = Z(J);
      j = nearfl(zi, Zj); jj = j + 1;
      zij = zi(j); zijj = zi(jj); dzj = zijj - zij;
      bmi(J,k) = ((Zj - zij)./dzj).*boi(j) + ((zijj - Zj)./dzj).*boi(jj);
    end
  end
  Bm(:,:,i) = bmi;

end

b1m = Bm(:,:,1);
b2m = Bm(:,:,2);
b3m = Bm(:,:,3);
b4m = Bm(:,:,4);
b5m = Bm(:,:,5);

end

function idxfl = nearfl(x, x0)
  % Get the lowest index of the two points in
  % vector x that bound the number x0.
  dx = abs(x - x0);
  flr = mink(dx, 2);
  il = find(dx==flr(1));
  ir = find(dx==flr(2));
  idxfl = min([il ir]);
end
