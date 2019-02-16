function [Vx, Vy, Vz] = janus2xyz(b1, b2, b3, b4, theta, varargin)
% Usage
% -----
% [Vx, Vy, Vz] = janus2xyz(b1, b2, b3, b4, theta)
%
%    OR
%
% [Vx, Vy, Vz] = janus2xyz(b1, b2, b3, b4, theta, 'use3BeamSol', true|false, ...
%                         'Binmap', BinmappingType, ptch, roll, r)
%
% Where 'BinmappingType' is the type of bin mapping to apply on the raw beam-coordinate velocities.
%
% INPUTS
% -------
% b1, b2, b3, b4                   [nz -by- nt] matrices of along-beam velocity components.
%
% theta                            angle between local vertical and the direction
%                                  of the Janus beams.
%
% Binmap ['linear' or 'nearest']   Whether to map the beam velocities to fixed horizontal
%                                  planes with linear interpolation ('linear') or
%                                  nearest-neighbor interpolation ('nearest') prior to converting
%                                  to instrument coordinates (Ott, 2002; Dewey & Stringer, 2007).
%                                  *The default is to NOT perform any bin mapping.
%
% use3BeamSol [true or false]      Whether to use three-beam solutions when exactly one beam has
%                                  no data in one cell.
%
% OUTPUTS
% -------
% [Vx, Vy, Vz]        [x, y, z]-ward components of instrument-referenced velocity vector.
%                               Vz is the z velocity component considering the conventional form
%                               of w (average of the estimates of w from planes 1-2 and 3-4").
%
%                     x-axis:   Increases in beam 1's direction, away from instrument.
%                     y-axis:   Increases in beam 3's direction, away from instrument.
%                     z-axis:   Increases upward direction, away from instrument.
options = struct('Binmap', 'none', 'use3BeamSol', false);
optionNames = fieldnames(options); % read the acceptable names.
nArgs = length(varargin);          % count arguments.

if ~isempty(varargin)
  use3BeamSol = varargin{find(strcmp(varargin, 'use3BeamSol'))+1};
  BinmapType = varargin{find(strcmp(varargin, 'Binmap'))+1};
  if ~strcmp(BinmapType, 'none')

    if nArgs<4
       error('Need pitch/roll angles and along-beam coordinate for bin-mapping.')
    end

    if isempty(use3BeamSol)
      ptch = varargin{3};
      roll = varargin{4};
      r = varargin{5};
      use3BeamSol = options.use3BeamSol;
    else
      ptch = varargin{5};
      roll = varargin{6};
      r = varargin{7};
  end

  end
else
  BinmapType = options.BinmapType;
  use3BeamSol = options.use3BeamSol;
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

if use3BeamSol==true
  [b1, b2, b3, b4] = janus3beamsol(b1, b2, b3, b4);
end

[Nz Nt] = size(b1);
B = cat(3, b1, b2, b3, b4);
d2r = pi/180;
theta = theta.*d2r;
uvfac = 1./(2.*sin(theta));
wfac = 1./(4.*cos(theta)); % For w derived from beams 1-4.

% 3rd row: w from the average of the 4 Janus beams.
%     b1 b2 b3 b4
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
