function [b1, b2, b3, b4, b5] = janus5earth2beam(u, v, w, theta)
% USAGE
% -----
%
uvw = cat(3, u, v, w);
d2r = pi/180;
theta = theta.*d2r;
[Nz Nt] = size(b1);

Ae2xyz = [-0.5   0    -0.25     0;
           0.5   0    -0.25     0;
            0  -0.5   -0.25     0;
            0   0.5   -0.25     0]; % pinv(A(:,1:end-1)), where A is the beam-to-xyz (i.e.,
                                    % instrument) transformation matrix (excluding w5).

Axyz2b = [];

u = +Vx.*cx1 + Vy.*cy1 + Vz.*cz1;
v = -Vx.*cx2 + Vy.*cy2 - Vz.*cz2;
w = -Vx.*cx3 + Vy.*cy3 + Vz.*cz3;

A = [+cx1 +cy1 +cz1
     -cx2 +cy2 -cz2
     -cx3 +cy3 +cz3];

B = zeros(Nz, Nt, 5);
for nz=1:Nz
  disp(['Calculating Vx, Vy, Vz at bin ', num2str(nz), '/', num2str(Nz)])
  for nt=1:Nt
    vxyz = Ae2xyz*squeeze(B(nz,nt,:)); % 1) From Earth to instrument coordinates.
    B(nz,nt,:) = Axyz2b*vxyz % 2) From instrument to beam coordinates.
  end
end

B = Ainv*uvw;
end
