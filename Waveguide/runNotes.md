# Run notes

| UID     | Parameters/Description | Results comments |
--------- | ---------------------- | ---------------- |
138881 | Lx=5 um, Lz=500 um. dx=3nm, dz=3nm | Looks OK with hdfview |
765719 | Lx=5 um, Lz=500 um. dx=3nm, dz=3nm, sparse save | File deleted |
22380 | Lx=5 um, Lz=500 um. dx=3nm, dz=3nm, sparse save, higher threshold than 765716 | File deleted |
446080 | Lx=5 um, Lz=500 um. dx=3nm, dz=3nm R=80mm | File deleted |
457131 | Lx=5 um, Lz=500 um. dx=3nm, dz=3nm R=80mm | HDF5 file corrupted. File is deleted |
964284 | Lx=5 um, Lz=500 um. dx=3nm, dz=10nm R=80mm | Seg. fault. |
452330 | Lx=5 um, Lz=500 um. dx=3nm, dz=10nm R=80mm | Looks good |
297428 | Lx=5 um, Lz=500 um. dx=3nm, dz=100nm R=80mm | Looks good, plottable in python :) |
876401 | Lx=5 um, Lz=500 um. dx=3nm, dz=100nm R=40mm | Wave is lost at approx. z = 0.4 mm. Files deleted. |
205961 | Lx=5 um, Lz=500 um. dx=3nm, dz=50nm R=40mm | Wave is lost at approx. z=0.4 mm. Files deleted. |
651953 | Lx=2.5 um, Lz=500 um. dx=3nm, dz=50nm R=80mm | Looks very good |
392705 | Lx=2.5 um, Lz=500 um. dx=3nm, dz=50nm R=40mm, delta=4.49E-5 | Wrong parameters. Files deleted. |
852045 | Lx=5.0 um, Lz=500 um. dx=3nm, dz=50nm R=40mm, delta=4.49E-5 | The beam is lost around 0.4 mm. Files deleted |
385496 | Lx=5.0 um, Lz=500 um. dx=3nm, dz=50nm R=40mm, delta=4.49E-5, try with full circle formula | The beam is lost around 0.4 mm. Files deleted |
263117 | Lx=5.0 um, Lz=500 um. dx=1nm, dz=50nm R=40mm, delta=4.49E-5 | Looks very, very good |
972076 | Lx=5.0 um, Lz=500 um. dx=1nm, dz=50nm R=40mm, delta=4.49E-5 | Forgot to save the transmission. Files deleted. |
715235 | Lx=5.0 um, Lz=500 um. dx=1nm, dz=50nm R=40mm, delta=4.49E-5 | Errors in the control file. Files are deleted |
516123 | Lx=5.0 um, Lz=500 um. dx=1nm, dz=100nm R=40mm, delta=4.49E-5 | Looks good, Damping: 2.6 mm |
751007 | Lx=5.0 um, Lz=200 um. dx=1nm, dz=50nm R=20mm, delta=4.49E-5 | Forgot to copy executable to server. Files are deleted. |
799571 | Lx=5.0 um, Lz=200 um. dx=1nm, dz=50nm R=20mm, delta=4.49E-5 | Too low stepsize when computing transmission. Filed deleted  |
776143 | Lx=2.0 um, Lz=200 um. dx=1nm, dz=50nm R=20mm, delta=4.49E-5 | Too low stepsize when computing transmission. Files deleted |
519783 | Lx=2.0 um, Lz=200 um. dx=1nm, dz=50nm R=20mm, delta=4.49E-5 | Looks good. Damping: 1.18 mm |
821665 | Lx=2.5 um, Lz=300 um. dx=1nm, dz=50nm R=30mm, delta=4.49E-5 | Looks good. Damping: 1.83 mm |
628369 | Lx=5.0 um, Lz=500 um. dx=1nm, dz=100nm R=40mm, delta=4.14E-5 | Error with HDF5 and Armadillo |
272786 | Lx=5.0 um, Lz=500 um. dx=1nm, dz=100nm R=40mm, delta=4.14E-5 | T=0.79. L=2.13 mm. Think T. Salditt uses delta=4.49E-5 |
181476 | Lx=5.0 um, Lz=500 um. dx=1nm, dz=100nm R=40mm, delta=4.49E-5 | Not analysed |
%%% | __Parameter sweep over R__ | %%% |
117646 | Lz/R=0.01 max(dx) =1.0, dz=calculated, delta=4.45E-5, R=10 mm | Looks Ok |
559170 | R = 20 mm | Looks OK |
240280 | R = 30 mm | Looks OK |
647340 | R = 40 mm | Looks OK |
920626 | R = 60 mm | Looks OK |
992235 | R = 80 mm | Looks OK |
605795 | R = 100 mm | To large stepsize in z? |
187870 | R = 150 mm | Too large stepsize in z? |
109596 | R = 200 mm | Too large stepsize in z? |
%%% | __Parameter sweep over R__ | %%% |
248375 |  Lz/R=0.01 max(dx) = 1 nm, max(dz)=100.0 nm, delta=4.45E-5, R=10 mm | Looks Ok |
55170 | R = 20 mm | Looks OK |
493127 | R = 30 mm | Looks OK |
253932 | R = 40 mm | Looks OK |
859966 | R = 60 mm | Looks OK |
91570 | R = 80 mm | Looks OK |
184063 | R = 100 mm | Looks OK |
866239 | R = 150 mm | Looks OK |
%%% | %%% | %%%
293663 | Lz/R=0.01 max(dx) =1.0, dz=calculated, delta=4.45E-5, Straight waveguide | L=5.35 mm, way too wide domain in the transverse direction |
702004 | Lz/R=0.01 max(dx) =1.0, dz=calculated, delta=4.45E-5, Straight waveguide | L = 5.32 mm, looks good |
966449 | Lz/R=0.01 max(dx) =1.0, dz=calculated, delta=4.45E-5, Straight waveguide, with real part | L = 5.32 mm, looks good |
%%% | __Parameter sweep over R saving the real part of the solution in addition__ | %%% |
694357 | Lz/R=0.01 max(dx) = 1 nm, max(dz)=100.0 nm, delta=4.45E-5, R=10 mm | Looks Ok |
725810 | R = 20 mm | Looks OK |
519089 | R = 30 mm | Looks OK |
824320 | R = 40 mm | Looks OK |
269826 | R = 60 mm | Looks OK |
262354 | R = 80 mm | Looks OK |
945091 | R = 100 mm | Looks OK |
521840 | R = 150 mm | Looks OK |
%%% | %%% | %%%
281636 | 1D eigenmode run. R=40mm. delta=4.49E-5 | Lookds good|
754680 | Coupled. R=80mm. delta=4.49E-5, sep=50 nm, couplerStart=1um | Difficult to see the second wg on log scale |
363907 | Coupled. R=80mm. delta=4.49E-5, sep=150 nm, couplerStart=1um | Good, need to close the upper for a longer time |
305415 | Coupled. R=80mm. delta=4.49E-5, sep=150 nm, couplerStart=50um | Ok, difficult to see |
379863 | Coupled. R=80mm, R2=2R1. delta=4.49E-5, sep=150 nm, couplerStart=50um | Some signals are catched |
%%% | %%% | %%%
747532 | Lz/R=0.01 max(dx) =1.0, dz=calculated, delta=4.45E-5, R=80mm, with far field | Ok. Some useful attributes where misisn. Files deleted. |
805696 | Lz/R=0.01 max(dx) =1.0, dz=calculated, delta=4.45E-5, R=80mm, with far field | Running |
255950 | 1D eigenmode run. R=40mm. delta=4.49E-5, 20 modes | Lookds good|
192712 | R=40mm. delta=4.49E-5 Gaussian beam. z0=-0.um, x0=0.0, w0=10nm | Looks OK |
254910 | R=40mm. delta=4.49E-5 Gaussian beam. z0=-0.um, x0=50nm, w0=10nm | Looks OK |
360525 | R=40mm. delta=4.49E-5 Plane wave. AngleZ=0.2 | Looks OK |
384176 | R=40mm. delta=4.49E-5 Gaussian beam. z0=-0.um, x0=50nm, w0=1 nm | Looks OK |
