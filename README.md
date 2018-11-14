# ADCPtools
Beam to Earth coordinate transformations and other utilities for Acoustic Doppler Current Profiler (ADCP) data.

## Convert profiles of along-beam velocity to Earth coordinates (east-north-up).
Convert raw along-beam velocity (positive toward transducer) time series (ntimes x nbins) to Earth-referenced velocity:
```
[u, v, w, w5] = janus5beam2earth(head, ptch, roll, theta, b1, b2, b3, b4, b5)
```

or, with all optional arguments specified,
```
[u, v, w, w5] = janus5beam2earth(head, ptch, roll, theta, b1, b2, b3, b4, b5, ...
                                 'uvwBeam5', true, 'Gimbaled', true, 'Binmap', true, rz)
```

Where

* ```b[1-5]``` are the along-beam velocity time series [ntimes x nbins] (positive toward transducer face).

* ```u```, ```v```, ```w``` are the Earth-referenced (eastward, northward, upward) velocity time series [ntimes x nbins].

* ```w5``` is the vertical velocity calculated from the vertical beam only.

* ```head```, ```pitch``` and ```roll``` are the heading (rotation about the z-axis (beam 5), positive clockwise), pitch (rotation about the x-axis (beam 1), positive clockwise) and roll (rotation about the y-axis (beam 3), positive clockwise) angles measured by the ADCP, in degrees.

* ```theta``` is the beam angle measured from beam 5's direction (vertical if the instrument is level), in degrees.

* ```'uvwBeam5'``` (```true```/```false```) indicates whether to use the vertical beam data to calculate ```u```, ```v``` and ```w```.

* ```'Gimbaled'``` (```true```/```false```) indicates whether the pitch and roll sensors were gimbaled (mounted on a free-swiveling platform as opposed to a rigid frame).

* ```'Binmap'``` (```'none'```/```'linear'```) indicates whether or not to perform bin-mapping on the raw beam velocities. Bin-mapping here means interpolating the velocities from all beams to the same horizontal plane before converting to instrument and then Earth coordinates. Defaults to ```'none'```, _i.e._, the assumption that pitch and roll are zero at all times.

* ```rz``` is a vector with the along-beam positions of the center of each bin. Required only if ```'Binmap'``` is ```true```.
