import numpy as np

z = 100.0
fname = "r_obs.txt"

x = np.linspace(-z,z,101)
y = np.linspace(-z,z,101)

out = open( fname, 'w' )
for i in range(0,len(x)):
    for j in range(0,len(y)):
        out.write("%.2f %.2f %.2f\n"%(x[i],y[j],z))
out.close()
