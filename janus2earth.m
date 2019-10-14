function [u, v, w] = janus2earth(head, ptch, roll, theta, b1, b2, b3, b4, varargin)
% USAGE
% -----
% [u, v, w] = janus2earth(head, ptch, roll, theta, b1, b2, b3, b4)
%
%                              OR
%
% [u, v, w] = janus2earth(head, ptch, roll, theta, b1, b2, b3, b4, r, ...
%                         Gimbaled, BinmapType, use3BeamSol)
%
% Calculates Earth velocities (u,v,w) = (east,north,up) from beam-referenced velocity time series
% from a 4-beam Janus ADCP, (e.g., Appendix A of Dewey & Stringer (2007), Equations A3-A11).
%
% nz, nt, nb = number of vertical bins, data records, beams.
%
%============================================================================
% For TRDI instruments, call function like this:
% [u, v, w] = janus2earth(head, ptch, roll, theta, b1, b2, b3, b4)
%
% For Nortek instruments, call function like this:
% [u, v, w] = janus2earth(head-90, roll, -ptch, theta, -b1, -b3, -b4, -b2)
%============================================================================
%
%    TRDI CONVENTION:
%    ================
%
% * Velocity toward transducers' faces: POSITIVE
% * Clockwise PITCH (tilt about x-AXIS): POSITIVE (beam 3 higher than beam 4)
% * Clockwise ROLL (tilt about y-AXIS):  POSITIVE (beam 2 higher than beam 1)
%
% * Heading increases CLOCKWISE from the *Y-AXIS*.
%
%       ^ positive y axis, psi = 0
%       |
%       3
%       |
%       |
%       |
% 2 --- O --- 1 ---> positive x axis, psi = +90
%       |
%       |
%       |
%       4
%
%    NORTEK CONVENTION:
%    ==================
%
% * Velocity toward transducers' faces: NEGATIVE
% * Counter-clockwise PITCH (tilt about y-AXIS, equivalent to -ROLL in the TRDI convention): POSITIVE (beam 1 higher than beam 3)
% * Clockwise ROLL (tilt about x-AXIS, equivalent to PITCH in the TRDI convention):  POSITIVE (beam 4 higher than beam 2)
%
% Heading increases CLOCKWISE from the *X-AXIS*.
%
%       ^ positive y axis, psi = -90
%       |
%       4
%       |
%       |
%       |
% 3 --- O --- 1 ---> positive x axis, psi = 0
%       |
%       |
%       |
%       2
%
% INPUTS
% ------
% b1, b2, b3, b4        [nz -by- nt] matrices of along-beam velocity components.
% head, ptch, roll      [nt]         vectors with (time-dependent) heading, pitch
%                                    and roll angles, following D&S2007's notation.
%
% theta                              Beam angle measured from the vertical.
%                                    *For RDI Sentinel V and Nortek Signature: 25.
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
% use3BeamSol [true or false]      Whether to use three-beam solutions when exactly one beam has
%                                  no data in one cell.
%
% Code for parsing named options as inputs from: https://stackoverflow.com/questions/2775263/
% how-to-deal-with-name-value-pairs-of-function-arguments-in-matlab
%
% OUTPUTS
% -------
% [u, v, w]           [east, north, up, up-(from vertical beam only)] components
%                      of Earth-referenced velocity vector.
if length(varargin)>0
  r = varargin{1};
  Gimbaled = varargin{2};
  BinmapType = varargin{3};
  use3BeamSol = varargin{4};
else
  Gimbaled = true;
  BinmapType = 'none';
  use3BeamSol = false;
end

nz = size(b1, 1);              % Number of vertical bins.
nt = size(b1, 2);              % Number of records in the time series.

d2r = pi/180;
head = head.*d2r; ptch = ptch.*d2r; roll = roll.*d2r;
theta = theta.*d2r;

% Time-dependent angles (heading, pitch and roll).
Sph1 = sin(head);
Sph2 = sin(ptch);
Sph3 = sin(roll);
Cph1 = cos(head);
Cph2 = cos(ptch);
Cph3 = cos(roll);

if Gimbaled==true % Correct heading (D&S 2007, eq. A2).
  disp('Gimbaled instrument case.')
  Sph2Sph3 = Sph2.*Sph3;
  head = head + asin( Sph2Sph3./sqrt(Cph2.^2 + Sph2Sph3.^2) );
  Sph1 = sin(head);
  Cph1 = cos(head);
else                      % Correct pitch (D&S 2007, eq. A1; Lohrmann et al. 1990, eq. A1).
  disp('Fixed instrument case.')
  ptch = asin( (Sph2.*Cph3)./sqrt(1 - (Sph2.*Sph3).^2) );
  Sph2 = sin(ptch);
  Cph2 = cos(ptch);
end

% Convert instrument-referenced velocities
% to Earth-referenced velocities.
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
if strcmp(BinmapType, 'none')
  [Vx, Vy, Vz] = janus2xyz(b1, b2, b3, b4, theta);
else
  [Vx, Vy, Vz] = janus2xyz(b1, b2, b3, b4, theta, r, ptch, roll, BinmapType, use3BeamSol);
end

u = +Vx.*cx1 + Vy.*cy1 + Vz.*cz1;
v = -Vx.*cx2 + Vy.*cy2 - Vz.*cz2;
w = -Vx.*cx3 + Vy.*cy3 + Vz.*cz3;

end
