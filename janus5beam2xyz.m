function [Vx, Vy, Vz, Vz5] = janus5beam2xyz(b1, b2, b3, b4, b5, theta, varargin)
% Usage
% -----
% [Vx, Vy, Vz, Vz5] = janus5beam2xyz(b1, b2, b3, b4, b5, theta, ...
%                                     r, ptch, roll, BinmapType, use3BeamSol)
%
% Where 'BinmapType' is the type of bin mapping to perform on the raw beam velocities.
%
% theta, ptch and roll in RADIANS.
%
% INPUTS
% -------
% b1, b2, b3, b4, b5  [nz -by- nt] matrices of along-beam velocity components.
%
% theta                            angle between local vertical and the direction
%                                  of the Janus beams.
%
% BinmapType ['linear' or 'nearest' or 'none']
%                                  Whether to map the beam velocities to fixed horizontal
%                                  planes with linear interpolation ('linear') or
%                                  nearest-neighbor interpolation ('nearest') prior to converting
%                                  to instrument coordinates (Ott, 2002; Dewey & Stringer, 2007).
%
% OUTPUTS
% -------
% [Vx, Vy, Vz, Vz5]   [x, y, z]-ward components of instrument-referenced velocity vector.
%                        ** Vz5 is the z velocity component from beam 5 only.**
%                        Vz is the z velocity component considering the conventional form
%                        of w (average of the estimates of w from planes 1-2 and 3-4).
%
%                     x-axis:   Increases in beam 1's direction, away from instrument.
%                     y-axis:   Increases in beam 3's direction, away from instrument.
%                     z-axis:   Increases in beam 5's direction, away from instrument.
if length(varargin)>0
  r = varargin{1};
  ptch = varargin{2};
  roll = varargin{3};
  BinmapType = varargin{4};
  use3BeamSol = varargin{5};
else
  BinmapType = 'none';
  use3BeamSol = false;
end

if strcmp(BinmapType, 'none')
  disp('Bin mapping NOT applied.')
else
  disp(['Mapping bins to horizontal planes with *',BinmapType,'* interpolation.'])
  [b1, b2, b3, b4, b5] = binmap5(b1, b2, b3, b4, b5, r, theta, ptch, roll, BinmapType);
end

if use3BeamSol==true
  [b1, b2, b3, b4] = janus3beamsol(b1, b2, b3, b4);
end

[Nz Nt] = size(b1);
B = cat(3, b1, b2, b3, b4, b5);
uvfac = 1./(2.*sin(theta));
wfac = 1./(4.*cos(theta)); % For w derived from beams 1-4.

% 3rd row: w from the average of the 4 Janus beams.
% 4rd row: w from the 5th beam only.
A = [-1  1  0  0  0;
      0  0 -1  1  0;
     -1 -1 -1 -1  0;
      0  0  0  0 -1];

vxyz = zeros(Nz, Nt, 4);
for nz=1:Nz
  disp(['Calculating Vx, Vy, Vz at bin ', num2str(nz), '/', num2str(Nz)])
  for nt=1:Nt
    vxyz(nz,nt,:) = A*squeeze(B(nz,nt,:));
  end
end

Vx = vxyz(:,:,1)*uvfac;
Vy = vxyz(:,:,2)*uvfac;
Vz = vxyz(:,:,3)*wfac;
Vz5 = vxyz(:,:,4);

end
