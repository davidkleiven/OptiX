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
651953 | Lx=2.5 um, Lz=500 um. dx=3nm, dz=50nm R=80mm | Looks very good. Files deleted |
392705 | Lx=2.5 um, Lz=500 um. dx=3nm, dz=50nm R=40mm, delta=4.49E-5 | Wrong parameters. Files deleted. |
852045 | Lx=5.0 um, Lz=500 um. dx=3nm, dz=50nm R=40mm, delta=4.49E-5 | The beam is lost around 0.4 mm. Files deleted |
385496 | Lx=5.0 um, Lz=500 um. dx=3nm, dz=50nm R=40mm, delta=4.49E-5, try with full circle formula | The beam is lost around 0.4 mm. Files deleted |
263117 | Lx=5.0 um, Lz=500 um. dx=1nm, dz=50nm R=40mm, delta=4.49E-5 | Looks very, very good. Files deleted |
972076 | Lx=5.0 um, Lz=500 um. dx=1nm, dz=50nm R=40mm, delta=4.49E-5 | Forgot to save the transmission. Files deleted. |
715235 | Lx=5.0 um, Lz=500 um. dx=1nm, dz=50nm R=40mm, delta=4.49E-5 | Errors in the control file. Files are deleted |
516123 | Lx=5.0 um, Lz=500 um. dx=1nm, dz=100nm R=40mm, delta=4.49E-5 | Looks good, Damping: 2.6 mm. Files deleted |
751007 | Lx=5.0 um, Lz=200 um. dx=1nm, dz=50nm R=20mm, delta=4.49E-5 | Forgot to copy executable to server. Files are deleted. |
799571 | Lx=5.0 um, Lz=200 um. dx=1nm, dz=50nm R=20mm, delta=4.49E-5 | Too low stepsize when computing transmission. Filed deleted  |
776143 | Lx=2.0 um, Lz=200 um. dx=1nm, dz=50nm R=20mm, delta=4.49E-5 | Too low stepsize when computing transmission. Files deleted |
519783 | Lx=2.0 um, Lz=200 um. dx=1nm, dz=50nm R=20mm, delta=4.49E-5 | Looks good. Damping: 1.18 mm |
821665 | Lx=2.5 um, Lz=300 um. dx=1nm, dz=50nm R=30mm, delta=4.49E-5 | Looks good. Damping: 1.83 mm |
628369 | Lx=5.0 um, Lz=500 um. dx=1nm, dz=100nm R=40mm, delta=4.14E-5 | Error with HDF5 and Armadillo |
272786 | Lx=5.0 um, Lz=500 um. dx=1nm, dz=100nm R=40mm, delta=4.14E-5 | T=0.79. L=2.13 mm. Think T. Salditt uses delta=4.49E-5 |
181476 | Lx=5.0 um, Lz=500 um. dx=1nm, dz=100nm R=40mm, delta=4.49E-5 | Not analysed |
%%% | __Parameter sweep over R__ (All files are deleted) | %%% |
117646 | Lz/R=0.01 max(dx) =1.0, dz=calculated, delta=4.45E-5, R=10 mm | Looks Ok |
559170 | R = 20 mm | Looks OK |
240280 | R = 30 mm | Looks OK |
647340 | R = 40 mm | Looks OK |
920626 | R = 60 mm | Looks OK |
992235 | R = 80 mm | Looks OK |
605795 | R = 100 mm | To large stepsize in z? |
187870 | R = 150 mm | Too large stepsize in z? |
109596 | R = 200 mm | Too large stepsize in z? |
%%% | __Parameter sweep over R__ (All files are deleted) | %%% |
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
%%% | __Parameter sweep over R saving the real part of the solution in addition__ (All files deleted)| %%% |
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
754680 | Coupled. R=80mm. delta=4.49E-5, sep=50 nm, couplerStart=1um | Difficult to see the second wg on log scale. Files deleted |
363907 | Coupled. R=80mm. delta=4.49E-5, sep=150 nm, couplerStart=1um | Good, need to close the upper for a longer time. Files deleted |
305415 | Coupled. R=80mm. delta=4.49E-5, sep=150 nm, couplerStart=50um | Ok, difficult to see. Files deleted |
379863 | Coupled. R=80mm, R2=2R1. delta=4.49E-5, sep=150 nm, couplerStart=50um | Some signals are catched. Files deleted |
%%% | %%% | %%%
747532 | Lz/R=0.01 max(dx) =1.0, dz=calculated, delta=4.45E-5, R=80mm, with far field | Ok. Some useful attributes where misisn. Files deleted. |
805696 | Lz/R=0.01 max(dx) =1.0, dz=calculated, delta=4.45E-5, R=80mm, with far field | Not analysed |
255950 | 1D eigenmode run. R=40mm. delta=4.49E-5, 20 modes | Lookds good|
192712 | R=40mm. delta=4.49E-5 Gaussian beam. z0=-0.um, x0=0.0, w0=10nm | Looks OK |
254910 | R=40mm. delta=4.49E-5 Gaussian beam. z0=-0.um, x0=50nm, w0=10nm | Looks OK |
360525 | R=40mm. delta=4.49E-5 Plane wave. AngleZ=0.2 | Looks OK |
384176 | R=40mm. delta=4.49E-5 Gaussian beam. z0=-0.um, x0=50nm, w0=1 nm | Looks OK |
384533 | R=40mm. delta=4.49E-5 Plane cylindrical coordinates z0=-0.um, x0=0.0 | Looks very good OK |
416996 | R=40mm. delta=4.49E-5 Plane cylindrical coordinates z0=-0.um, x0=0.0, proper axis labels | Looks very good OK |
968548 | R=10mm. delta=4.49E-5 Plane cylindrical coordinates z0=-0.um, x0=0.0 | Looks very good OK |
649695 | R=40mm. delta=4.49E-5 Plane cartesian coordinates z0=-0.um, x0=0.0, simple averagins | Still hairs, simple averaging has no effect, error with transmission |
%%% | __Parameter sweep over R saving the real part of the solution in addition (Cylindrical crd)__ | %%% |
241802 | Lz=500 um max(dx) = 1 nm, max(dz)=100.0 nm, delta=4.45E-5, R=10 mm, inc angle=0.2 deg | Looks Ok |
929539 | R = 20 mm | Looks OK |
735044 | R = 30 mm | Looks OK |
873746 | R = 40 mm | Looks OK |
510399 | R = 60 mm | Looks OK |
668643 | R = 80 mm | Looks OK |
729962 | R = 100 mm | Looks OK |
788708 | R = 150 mm | Looks OK |
%%% | %%% | %%%
550454 | Beam split angle 0.015 | Looks OK |
918592 | Beam split angle 0.015, delta=1E-5, beta=1E-7 | Looks OK |
609003 | R=40mm delta=4.45E-5, beta=3.45E-6 gaussian WG| No hairs |
237703 | R=40mm delta=4.45E-5, beta=3.45E-6 linear ramp WG, f=1.0 | 3 hairs |
158240 | R=40mm delta=4.45E-5, beta=3.45E-6 linear ramp WG, f=2.0 | 3 hairs |
118543 | R=40mm delta=4.45E-5, beta=3.45E-6 linear ramp WG, f=5.0 | 3 hairs |
440622 | R=40mm delta=4.45E-5, beta=3.45E-6 step, w=5nm | No hairs |
577660 | R=40mm delta=4.45E-5, beta=3.45E-6 step, w=100nm, nz=nx=1000 | Hairy |
645627 | R=40mm delta=4.45E-5, beta=3.45E-6 step, w=100nm, nz=nx=100 | Smoothed out |
462501 | R=40mm delta=4.45E-5, beta=3.45E-6 step, w=100nm, nz=nx=500 | Hairs |
476450 | R=40mm delta=4.45E-5, beta=3.45E-6 step, w=100nm, nz=nx=2000 | Hairs |
869908 | R=40mm delta=4.45E-5, beta=3.45E-6 step, w=100nm, nz=nx=3000 | Hairs |
558362 | R=40mm delta=4.45E-5, beta=3.45E-6 step, w=100nm, nz=3000 nx=100 | Some interesting staircase |
272786 | R=40mm delta=4.45E-5, beta=3.45E-6 step, w=100nm, nz=3000 nx=200 | Interasting staircase features |
133595 | R=40mm delta=4.45E-5, beta=3.45E-6 step, w=100nm, nz=3000 nx=200, cylindrical, proper transmission | Looks OK |
352671 | R=40mm delta=4.45E-5, beta=3.45E-6 step, w=100nm, nz=1000 nx=1000, cartesian | Too high losses |
558175 | R=40mm delta=4.45E-5, beta=3.45E-6 step, w=100nm, nz=5000 nx=5000, cartesian | Looks good |
460643 | R=40mm delta=4.45E-5, beta=3.45E-6 step, w=100nm, nz=5000 nx=5000, cartesian, proper real part | Looks good |
121216 | R=40mm delta=4.45E-5, beta=3.45E-6 step, w=100nm, nz=5000 nx=5000, cylindrical, proper far field | Looks good |
268220 | R=40mm delta=4.45E-5, beta=3.45E-6 step, w=100nm, nz=5000 nx=5000, cartesian, proper far field | Looks good |
450776 | Inc angle sweep Nz=10000, Nx=1000 Lz=5mm | Looks very good |
772477 | Inc angle sweep Nz=10000, Nx=1000 Lz=5mm with alcohol | Looks very good |
509750 | Cylindrical proper far field Lz=500 um | Looks OK |
722020 | Cartesian proper far field Lz=400 um | Looks OK |
459397 | Cylindrical proper far field Lz=400 um | Looks OK |
682477 | Cartesian with b-track Lz=400 um | Looks OK, more equal to the other cartesian |
768375 | Cartesiann Nz=10000, Nx=2000, Lz=400um | Far field looks very good |
352700 | Eigenmode run, saving more parameters in the hdf5 file | Looks good |
174797 | Eigenmode run, 15 modes | Looks good |
172060 | Eigenmode run | Looks good |
363100 | Eigenmode run, 70 modes | Looks good |
577579 | Eigenmodes straight waveguide | Looks good |
953999 | Cylindrical, distance is widht on both side | Looks good |
563999 | Straight waveguide | Wrong boundary conditions |
478288 | Straight waveguide | Too short waveguide |
579186 | Straight waveguide | Looks good |
%%% | __Parameter sweep over R__ | %%% |
265156 |  Lz/R=0.01, delta=4.45E-5, R=10 mm, Nx=2000, Nz=10000 | Looks Ok |
963983 | R = 20 mm | Looks OK |
681603 | R = 30 mm | Looks OK |
827383 | R = 40 mm | Looks OK |
116355 | R = 60 mm | Looks OK |
507568| R = 80 mm | Looks OK |
680743 | R = 100 mm | Looks OK |
326866 | R = 150 mm | Looks OK |
%%% | %%% | %%% |
796988 | R = 40 mm, L = 1200 um, width=1200 nm, Nx=2000, Nz=3000 | Looks good, solution mimics the Eikonal equation |
%%% | __Parameter sweep over width___ | %%% |
114123  | w = 10 nm | No propagating modes
301088 | w = 20 nm | Look OK |
433090 | w = 30 nm | Looks OK |
857974 | w = 40 nm | Looks OK |
348856 | w = 60 nm | Looks OK |
539750 | w = 80 nm  | Looks OK |
763629 | w = 100 nm | Looks OK |
397296 | w = 150 nm | Looks OK |
%%% | __Parameter sweep over width__ | %%% |
326678 | w = 40 nm | Looks OK |
479615 | w = 60 nm | Looks OK |
943087 | w = 80 nm | Looks OK |
963353 | w = 100 nm | Looks OK |
524050 | w = 200 nm | Looks OK |
473836 | w = 400 nm | Looks OK |
903590 | w = 800 nm | Looks OK |
272447 | w = 1200 nm  | Looks OK |
%%% | %%% | %%% |
237813 | Incident angle sweep straight, no alcohol. width=69.8 nm | Looks good |
441159 | Incident angle sweep straight, alcohol, width=69.8 nm | Looks good |
269064 | Incident angle sweep, w=69.8 nm, E=10keV | Looks good |
936436 | Incident angle sweep, w=699.8 m, E=10keV | Looks good |
4289323 | Fresnel simulation of run 857974 |
497528 | Beam split | Look good |
9262170 | Fresnel two curved WG | Looks very good, can not be used have to use the cartesian |
629884 | Curved WG Lz=400 um, R=40mm with border tracker | Looks good
6043124 | Fresnel of two curved WG, using code with border tracker. Based on 629884. | Looks very good |
570321 | Eigenmode SiO2 core vaccuum, E=10 keV | Looks good |
304440 | Eigenmode SiO2 core C2H6O2 E=10keV | Looks good |
125047 | Lengths + angle sweep L=2.5 mm to 3.5 mm | Far field is very dependent on the length |
4439827 | Fresnel two curved WG, using code with border tracker. Based 629884. With padding to avoid periodic effects | Looks good |
899810 | Cylidnrial Nx=100, Nz=3000 | Looks OK |
477129 | Cylindrical Nx=50, Nz=3000 | Looks OK |
191298 | Cylindrical Nx=200, Nz=3000 | Not analysed |
581287 | Cylindrical Nx=400, Nz=3000 | Not analysed |
563966 | Cylindrical Nx=800, Nz=3000 | Not analysed |
504230 | Cartesian Nx=100, brdtrack, Nz=3000 | Not analysed |
943809 | Cartesian Nx=200, brdtrack, Nz=3000 | Not analysed |
171918 | Cartesian Nx=400, brdtrack, Nz=3000 | Not analysed |
224190 | Cartesian Nx=800, brdtrack, Nz=3000 | Not analysed |
868127 | Cartesian Nx=100, brdtrack, Nz=3000, with transmission | Not analysed |
631690 | Cartesian Nx=200, brdtrack, Nz=3000, with transmission | Not analysed |
976491 | Cartesian Nx=400, brdtrack, Nz=3000, with transmission | Not analysed |
496524 | Cartesian Nx=800, brdtrack, Nz=3000, with transmission | Not analysed |
