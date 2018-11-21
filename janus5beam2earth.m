function [u, v, w, w5] = janus5beam2earth(head, ptch, roll, theta, b1, b2, b3, b4, b5, varargin)
% USAGE
% -----
% [u, v, w, w5] = janus5beam2earth(head, ptch, roll, theta, b1, b2, b3, b4, b5,
%                                 'uvwBeam5', true, 'Gimbaled', true)
%
%                              OR
%
% [u, v, w, w5] = janus5beam2earth(head, ptch, roll, theta, b1, b2, b3, b4, b5,
%                                 'uvwBeam5', true, 'Gimbaled', true, 'Binmap', BinmappingType, rz)
%
% Calculates Earth velocities (u,v,w) = (east,north,up) from beam-referenced velocity time series
% from a 5-beam Janus ADCP, (e.g., Appendix A of Dewey & Stringer (2007), Equations A3-A11).
%
% nz, nt, nb = number of vertical bins, data records, beams.
%
%    TRDI CONVENTION:
%    ================
%
% * Velocity toward transducers' faces: POSITIVE
% * Counter-clockwise PITCH (tilt about x-AXIS): POSITIVE (beam 4 higher than beam 3)
% * Counter-clockwise ROLL (tilt about y-AXIS):  POSITIVE (beam 1 higher than 2)
%
%       ^ positive y axis, psi = 0
%       |
%       3
%       |
%       |
%       |
% 2 --- O --- 1 ---> positive x axis, psi = -90
%       |
%       |
%       |
%       4
%
%    NORTEK CONVENTION:
%    ==================
%
% * Velocity toward transducers' faces: NEGATIVE
% * Counter-clockwise PITCH (tilt about y-AXIS): POSITIVE (beam 1 higher than beam 3)
% * Counter-clockwise ROLL (tilt about x-AXIS):  NEGATIVE (beam 2 higher than beam 4)
%
%       ^ positive y axis, psi = 0
%       |
%       4
%       |
%       |
%       |
% 3 --- O --- 1 ---> positive x axis, psi = -90
%       |
%       |
%       |
%       2
%
% TRDI convention for beam numbering:   psi = [psi1 psi2 psi3 psi4] = [-90 90 0 180].
% Nortek convention for beam numbering: psi = [psi1 psi2 psi3 psi4] = [-90 180 90 0].
%
% INPUTS
% ------
% b1, b2, b3, b4, b5    [nz -by- nt] matrices of along-beam velocity components.
% head, ptch, roll      [nt]         vectors with (time-dependent) heading, pitch
%                                    and roll angles, following D&S2007's notation.
%
% theta                              Beam angle measured from the vertical.
%                                    *For RDI Sentinel V and Nortek Signature: 25.
%
% uvwBeam5     [true or false]   whether to calculate [u, v, w] using the independent information
%                                from beam 5 (defaults true). If false, the usual four-beam
%                                solution using w derived from beams 1-4 is calculated.
%
% Gimbaled     [true or false]    Whether the ADCP was deployed with a gimbaled roll sensor
%                                 (default true). Applies the correction to the raw pitch angle
%                                 if the pitch/roll sensors were mounted rigidly to the
%                                 instrument ('Gimbaled'==false), or the correction to the raw
%                                 heading angle if the ADCP was mounted on a gimbal (Dewey &
%                                 Stringer, 2007; Lohrmann et al., 1990).
%
% Binmap ['none' or 'linear' or 'nearest']
%                                  Whether to map the beam velocities to fixed horizontal
%                                  planes with linear interpolation ('linear') or nearest-neighbor
%                                  interpolation ('nearest') prior to converting
%                                  to instrument coordinates (Ott, 2002; Dewey & Stringer, 2007).
%                                  *The default is to NOT perform any bin mapping.
%
% Code for parsing named options as inputs from: https://stackoverflow.com/questions/2775263/
% how-to-deal-with-name-value-pairs-of-function-arguments-in-matlab
%
% OUTPUTS
% -------
% [u, v, w, w5]           [east, north, up, up-(from vertical beam only)] components
%                         of Earth-referenced velocity vector.
%
% For TRDI instruments, call function like this:
% [u, v, w] = janus2earth(head, ptch, roll, theta, b1, b2, b3, b4, b5)
%
% For Nortek instruments, call function like this:
% [u, v, w] = janus2earth(head, -roll, ptch, theta, -b1, -b3, -b4, -b2, -b5)
options = struct('uvwBeam5', true,'Gimbaled', true, 'Binmap', 'none', 'rz', NaN);
optionNames = fieldnames(options); % read the acceptable names.
nArgs = length(varargin);          % count arguments.

if ~isempty(varargin) && any(strcmp(varargin, 'Binmap'))
  BinmapType = varargin{find(strcmp(varargin, 'Binmap'))+1};
  if ~strcmp(BinmapType, 'none')

    if isnan(varargin{end})
       error('Need along-beam coordinate for bin-mapping.')
    else
      rz = varargin{end};
      varargin = varargin(1:end-1);
      nArgs = nArgs - 1;
    end

  end
end

if round(nArgs/2)~=nArgs/2
   error('Need propertyName/propertyValue pairs')
end

for pair = reshape(varargin,2,[])
  inpName = pair{1};
  if any(strcmp(inpName,optionNames))
    options.(inpName) = pair{2};
  else
    error('option %s is not recognized',inpName)
  end
end

nz = size(b1, 1);              % Number of vertical bins.
nt = size(b1, 2);              % Number of records in the time series.

d2r = pi/180;
head = head.*d2r; ptch = ptch.*d2r; roll = roll.*d2r;

% Time-dependent angles (heading, pitch and roll).
Sph1 = sin(head);
Sph2 = sin(ptch);
Sph3 = sin(roll);
Cph1 = cos(head);
Cph2 = cos(ptch);
Cph3 = cos(roll);

if options.Gimbaled==true % Correct heading (D&S 2007, eq. A2).
  disp('Gimbaled instrument case.')
  Sph2Sph3 = Sph2.*Sph3;
  head = head + asin( Sph2Sph3./sqrt(Cph2.^2 + Sph2Sph3.^2) );
else                      % Correct pitch (D&S 2007, eq. A1; Lohrmann et al. 1990, eq. A1).
  disp('Fixed instrument case.')
  ptch = asin( (Sph2.*Cph3)./sqrt(1 - (Sph2.*Sph3).^2) );
end

% Convert instrument-referenced velocities
% to Earth-referenced velocities.
% Option 1: Classic four-beam solution.
% Option 2: five-beam solution for [u, v, w].
cx1 = Cph1.*Cph3 + Sph1.*Sph2.*Sph3;
cx2 = Sph1.*Cph3 - Cph1.*Sph2.*Sph3;
cx3 = Cph2.*Sph3;
cy1 = Sph1.*Cph2;
cy2 = Cph1.*Cph2;
cy3 = Sph2;
cz1 = Cph1.*Sph3 - Sph1.*Sph2.*Cph3;
cz2 = Sph1.*Sph3 + Cph1.*Sph2.*Cph3;
cz3 = Cph2.*Cph3;

% Convert beam-referenced velocities to instrument-referenced velocities.
% NOTE: The convention used here (positive x axis = horizontally away from beam 1) and
%                                 positive y axis = horizontally away from beam 3) is not
%                                 the same as the one used by the instrument's firmware if
%                                 the coordinate transformation mode is set to "instrument
%                                 coordinates" before deployment.
if ~strcmp(options.Binmap, 'none')
  [Vx, Vy, Vz, Vz5] = janus5beam2xyz(b1, b2, b3, b4, b5, theta, 'Binmap', options.Binmap, ptch./d2r, roll./d2r, rz);
else
  [Vx, Vy, Vz, Vz5] = janus5beam2xyz(b1, b2, b3, b4, b5, theta);
end
w5 = Vz5.*cz3; % w from beam 5 only.

if options.uvwBeam5==true
  disp('Using vertical beam for [u, v, w].')
  u = +Vx.*cx1 + Vy.*cy1 + Vz5.*cz1;
  v = -Vx.*cx2 + Vy.*cy2 - Vz5.*cz2;
  w = -Vx.*cx3 + Vy.*cy3 + w5;
else
  disp('Using only beams 1-4 for [u, v, w].')
  u = +Vx.*cx1 + Vy.*cy1 + Vz.*cz1;
  v = -Vx.*cx2 + Vy.*cy2 - Vz.*cz2;
  w = -Vx.*cx3 + Vy.*cy3 + Vz.*cz3;
end
