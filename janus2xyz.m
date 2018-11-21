function [Vx, Vy, Vz] = janus2xyz(b1, b2, b3, b4, theta, varargin)
% Usage
% -----
% [Vx, Vy, Vz] = janus2xyz(b1, b2, b3, b4, theta)
%
%    OR
%
% [Vx, Vy, Vz] = janus2xyz(b1, b2, b3, b4, theta, 'Binmap', BinmappingType, ...
%                          ptch, roll, rz)
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
  if ~strcmp(varargin{1}, 'none')

    if nArgs<4
       error('Need pitch and roll angles and along-beam coordinate for bin-mapping.')
    end

    for pair = varargin(1:2)'
      inpName = pair{1};
      if any(strcmp(inpName,optionNames))
        options.(inpName) = pair{2};
      else
        error('option %s is not recognized',inpName)
      end
    end
    ptch = varargin{2};
    roll = varargin{3};
    rz = varargin{4};

  end
end

switch options.Binmap
case 'linear'
  disp('Mapping beams to horizontal planes with *linear* interpolation.')
  [b1, b2, b3, b4] = binmaplin(b1, b2, b3, b4, rz, theta, ptch, roll);
case 'nearest'
  disp('Mapping beams to horizontal planes with *nearest-neighbor* interpolation.')
  % Implement nn.
case 'none'
  disp('Bin mapping NOT applied.')
otherwise
  error(['Invalid bin mapping method: ' option.Binmap '.'])
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
