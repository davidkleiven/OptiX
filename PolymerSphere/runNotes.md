# Description of runs with different UID

UID    | Description                   | Result Comment                                                  |
-------|:-----------------------------:|-----------------------------------------------------------------|
848166 | kR = 5, eps = 1-1E-5 + 1E-6, Nodes=4099   | The results deviate from the exact solution         |
613242 | kR = 10, eps=1-1E-5 + 1E-6, Nodes=4099    | The results deviate from exact, probably discretisation |
321847 | kR = 5, eps=1-1E-5+1E-6, Nodes = 4099     | Forgot to do git pull before launching, so it is equation to run 848166 |
735305 | kR = 5, eps=1E-4+1E-6, Nodes = 4099 | Aborted, network error, screen were used so why? |
35754  | kR = 5. eps=1-1E-3 + 1E-6, Nodes = 4099 | Aborted, network error, screen were used  |
284528 | kR=5, eps=1E-4 + 1E-6i, Nodes = 4099, periodic only spheres | Out of memory while assembling why? |
200001 | lambda=0.7, eps=1-1E-5+1E-6i, Nodes=23302, MoM prec=0, VCFIE prec=0 | Puma run, the results show lobes" |
627234 | kR=5, eps=0.81+1E-6i, Nodes = 1027 | Very good results, short runtime |
176935 | kR=5, eps=0.98+1E-6i, Nodes = 1027 | Very good results, short runtime |
131925 | kR=5, eps=0.998 + 1E-6i, Nodes = 1027  | Good, deviates more than run 627234 and 176935 |
501846 | kR=10, eps=0.9999+1E-6i, Nodes = 771, cube | Deviates quite a bit |
161737 | kR=10, eps=0.99998+1E-6, Nodes = 771, cube | Deviates quite a bit |
562343 | kR=10, eps=0.998+1E-6i, Nodes = 771, cube | Very good results |
358059 | kR=5, eps=0.998+1E-6i, Nodes = 771, cube | Very good results |
625837 | kR=10, eps=0.998+1E-6i, Nodes = 3075, cube | Was apperently launched with kR=0 |
| 685739 | kR=5, eps=0.99998+1E-6i, Nodes = 3075, cube | Better, still it does not follow good near the tails |
| 272259 | kR=10, eps=0.99998+1E-6i, Nodes = 3075, cube | Results did not fit good |
| 345269 | kR=5, eps=0.99998+1E-6i, Nodes=4452, sphere | Not analysed, the results will probably be poor |
|200758 | kR=5, eps=0.995+0.005i, Nodes=739, sphere | Good results, can proceed with these parameters |
| 216156 | kR=5, eps=0.99500.005i, Nsphere=739, Nsubst=835, Grazing angle=90 | Induced currents seems OK |
| 16522 | kR=5, eps=0.995+0.005i, Nsphere=1016, Nsubst=1382, Grazing angle=90, Gaussian beam waist=3.0  | Forgot to do git pull, same as 216156 but with better resolution |
| 542786 | kR=5, eps=0.995+0.005i, Nsphere=1016, Nsubst=1382, Grazing angle=90, Gaussian beam waist=3.0 | Pattern mush narrower than expected |
| 184327 | kR=5, eps=0.995+0.005i, Nsphere=1016, Nsubst=1382, Grazing angle=90, Gaussian beam waist=1.0 | Pattern seems narrower than expected |
| 75347 | kR=5, eps=0.995+0.001i, N=739, Gaussian waist = 3.0 | Good results  |  
| 356437 | kR=5, eps=0.995+0.001i, Nsphere=739, Ncube=3248. Gaussian. waist=3.0 | Pattern is very narrow |
| 87373 | kR=5, eps=0.995+0.001i, Nsphere=739, Ncube=3248. Gaussian. waist=2.0 | Same as 356437 |
| 498264 | kR=5, eps=0.995+0.001i, Nsphere=739, Ncube=2305. Gaussian. waist=2.0 Thickness=0.2 | Observe an effect of the film |
| 961865 | kR=5, eps=0.995+0.001i, Nsphere=739, Ncube=2305. Plane wave | Looks like the film scatters a lot |
| 621734 | kR=5, eps=0.995+0.001i, Nsphere=739, Ncube=2305. Gaussian. waist=1.0 Thickness=0.2. Note: this beam diverges a lot | Running |
| 156499 | kR=5, eps=0.995+0.001i, Nsphere=739, Ncube=4409. Plane wave. Plane width=30 | Running |
