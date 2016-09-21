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
625837 | kR=10, eps=0.998+1E-6i, Nodes = 3075, cube | Running |
