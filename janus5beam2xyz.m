function [Vx, Vy, Vz, Vz5] = janus5beam2xyz(b1, b2, b3, b4, b5, theta)
% Usage
% -----
% [Vx, Vy, Vz, Vz5] = janus5beam2xyz(b1, b2, b3, b4, b5, theta)
%
% INPUTS
% -------
% b1, b2, b3, b4, b5  [nz -by- nt] matrices of along-beam velocity components.
% theta                            angle between local vertical and the direction
%                                  of the Janus beams.
%
% OUTPUTS
% -------
% [Vx, Vy, Vz, Vz5]   [x, y, z]-ward components of instrument-referenced velocity vector.
%  Vz                    ** Vz5 is the z velocity component from beam 5 only.**
%                        Vz is the z velocity component considering the conventional form
%                        of w (average of the estimates of w from planes 1-2 and 3-4).
%
%                     x-axis:   Increases in beam 1's direction, away from instrument.
%                     y-axis:   Increases in beam 3's direction, away from instrument.
%                     z-axis:   Increases in beam 5's direction, away from instrument.
B = cat(3, b1, b2, b3, b4, b5);
d2r = pi/180;
theta = theta.*d2r;
uvfac = 1./(2.*sin(theta));
wfac = 1./(4.*cos(theta)); % For w derived from beams 1-4.
[Nz Nt] = size(b1);

% 3rd row: w from the average of the 4 Janus beams.
% 4rd row: w from the 5th beam only.
A = [-1  1  0  0  0;
      0  0 -1  1  0;
     -1 -1 -1 -1  0;
      0  0  0  0 -1];

uvw = zeros(Nz, Nt, 4);
for nz=1:Nz
  disp(['Calculating Vx, Vy, Vz at bin ', num2str(nz), '/', num2str(Nz)])
  for nt=1:Nt
    uvw(nz,nt,:) = A*squeeze(B(nz,nt,:));
  end
end

Vx = uvw(:,:,1)*uvfac;
Vy = uvw(:,:,2)*uvfac;
Vz = uvw(:,:,3)*wfac;
Vz5 = uvw(:,:,4);

end
