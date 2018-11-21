function [Vx, Vy, Vz] = janus2xyz(b1, b2, b3, b4, theta, varargin)
% Usage
% -----
% [Vx, Vy, Vz] = janus2xyz(b1, b2, b3, b4, theta)
%
%    OR
%
% [Vx, Vy, Vz] = janus2xyz(b1, b2, b3, b4, theta, 'Binmap', BinmappingType, ...
%                          ptch, roll, r)
%
% Where 'BinmappingType' is the type of bin mapping to perform on the raw beam velocities.
%
% INPUTS
% -------
% b1, b2, b3, b4                   [nz -by- nt] matrices of along-beam velocity components.
%
% theta                            angle between local vertical and the direction
%                                  of the Janus beams.
%
% Binmap ['linear' or 'nearest']
%                                  Whether to map the beam velocities to fixed horizontal
%                                  planes with linear interpolation ('linear') or
%                                  nearest-neighbor interpolation ('nearest') prior to converting
%                                  to instrument coordinates (Ott, 2002; Dewey & Stringer, 2007).
%                                  *The default is to NOT perform any bin mapping.
%
% OUTPUTS
% -------
% [Vx, Vy, Vz]        [x, y, z]-ward components of instrument-referenced velocity vector.
%                               Vz is the z velocity component considering the conventional form
%                               of w (average of the estimates of w from planes 1-2 and 3-4, often
%                               called "error velocity").
%
%                     x-axis:   Increases in beam 1's direction, away from instrument.
%                     y-axis:   Increases in beam 3's direction, away from instrument.
%                     z-axis:   Increases upward direction, away from instrument.
options = struct('Binmap', 'none');
optionNames = fieldnames(options); % read the acceptable names.
nArgs = length(varargin);          % count arguments.

if ~isempty(varargin)
  BinmapType = varargin{find(strcmp(varargin, 'Binmap'))+1};
  if ~strcmp(BinmapType, 'none')

    if nArgs<4
       error('Need pitch and roll angles and along-beam coordinate for bin-mapping.')
    end

    ptch = varargin{3};
    roll = varargin{4};
    r = varargin{5};

  end
else
  BinmapType = 'none';
end

switch BinmapType
case 'linear'
  disp('Mapping beams to horizontal planes with *linear* interpolation.')
  [b1, b2, b3, b4] = binmap(b1, b2, b3, b4, r, theta, ptch, roll, 'linear');
case 'nn'
  disp('Mapping beams to horizontal planes with *nearest-neighbor* interpolation.')
  [b1, b2, b3, b4] = binmap(b1, b2, b3, b4, r, theta, ptch, roll, 'nn');
case 'none'
  disp('Bin mapping NOT applied.')
otherwise
  error(['Invalid bin mapping method: ' BinmapType '.'])
end

B = cat(3, b1, b2, b3, b4);
d2r = pi/180;
theta = theta.*d2r;
uvfac = 1./(2.*sin(theta));
wfac = 1./(4.*cos(theta)); % For w derived from beams 1-4.
[Nz Nt] = size(b1);

% 3rd row: w from the average of the 4 Janus beams.
A = [-1  1  0  0;
      0  0 -1  1;
     -1 -1 -1 -1];

vxyz = zeros(Nz, Nt, 3);
for nz=1:Nz
  disp(['Calculating Vx, Vy, Vz at bin ', num2str(nz), '/', num2str(Nz)])
  for nt=1:Nt
    vxyz(nz,nt,:) = A*squeeze(B(nz,nt,:));
  end
end

Vx = vxyz(:,:,1)*uvfac;
Vy = vxyz(:,:,2)*uvfac;
Vz = vxyz(:,:,3)*wfac;

end
