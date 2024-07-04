from matplotlib import pyplot as plt
import numpy as np

with open("./atlas_1007.dat", "r") as f:
    d = f.read()

d = d.split("\n")
keys = d[0].split()
values = [i.split() for i in d[1:-1]]
values = [i for i in values if i[5] == 'orange']

X = [float(i[0]) for i in values]
Y = [float(i[3]) for i in values]
plt.scatter(X, Y)

g1 = [i for i in values if float(i[0]) > 59470 and float(i[0]) < 59475]
g2 = [i for i in values if float(i[0]) > 59475 and float(i[0]) < 59477]
g3 = [i for i in values if float(i[0]) > 59477 and float(i[0]) < 59479]

print([[i[0] for i in g] for g in [g1, g2, g3]])

g1_binned = []
g2_binned = []
g3_binned = []

for i in range(3):
    group = [g1, g2, g3][i] 
    binned = [g1_binned, g2_binned, g3_binned][i]
    
    for j in range(len(group[0])):
        try:
            l = np.mean([float(g[j]) for g in group])
        except:
            l = group[0][j]
        binned.append(l)

print([g[0] for g in [g1_binned, g2_binned, g3_binned]])

x1 = [float(g[0]) for g in g1]
y1 = [float(g[3]) for g in g1]
dy1 = [float(g[4]) for g in g1]
x2 = [float(g[0]) for g in g2]
y2 = [float(g[3]) for g in g2]
dy2 = [float(g[4]) for g in g2]
x3 = [float(g[0]) for g in g3]
y3 = [float(g[3]) for g in g3]
dy3 = [float(g[4]) for g in g3]

x = x1 + x2 + x3
y = y1 + y2 + y3
dy = dy1 + dy2 + dy3
print(x, y, dy)

xx1 = g1_binned[0] 
yy1 = g1_binned[3] 
dyy1 = g1_binned[4] 
xx2 = g2_binned[0] 
yy2 = g2_binned[3] 
dyy2 = g2_binned[4] 
xx3 = g3_binned[0] 
yy3 = g3_binned[3]
dyy3 = g3_binned[4] 

xx = [xx1, xx2, xx3]
yy = [yy1, yy2, yy3]
dyy = [dyy1, dyy2, dyy3]
print(xx, yy, dyy)

plt.scatter(x, y)
plt.errorbar(x, y, yerr=dy, linestyle="none")
plt.scatter(xx, yy)
plt.errorbar(xx, yy, yerr=dyy, linestyle="none")
plt.show()

data = [keys, g1_binned, g2_binned, g3_binned]
data = [[str(i) for i in j] for j in data]
data = [" ".join(i) for i in data]
data = "\n".join(data)
print(data)

with open("./atlas_1007_binned.dat", "w") as f:
    f.write(data)
