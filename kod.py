import numpy
import matplotlib.pyplot as plt

'''
Lagrange
'''

data = numpy.array([[1., 3.],[2., 1.],[3.5, 4.],[5., 0.],[6., .5],[9., -2.],[9.5, -3.]])
N = data.shape[0]
x = numpy.linspace(0, 10, 1000)


#Znajdujemy wielomiany bazy
lagrange_poly = numpy.ones((N, x.shape[0]))
for i in range(N):
    for j in range(N):
        if j!=i:
            lagrange_poly[i,:] *= (x-data[j,0])/(data[i,0]-data[j,0])
            print(i, data[i, 0], data[j, 0])
            
#Znajdujemy wzór Lagrange'a
lagrange = numpy.zeros(x.shape[0]) 
for n in range(N):
    lagrange += lagrange_poly[n, :] * data[n, 1]


'''
Sklejana 1 stopnia
'''
a = numpy.zeros(data.shape[0])

for i in range(data.shape[0] - 1):
    a[i] = (data[i+1, 1] - data[i, 1]) / (data[i+1, 0] - data[i, 0])


sklejana1st = numpy.zeros(x.shape[0])
for i in range(data.shape[0] - 1):
    sklejana1st += ((data[i+1, 1] - data[i, 1])/(data[i+1, 0] - data[i, 0]) * (x - data[i, 0]) + data[i, 1]) * (x < data[i+1, 0]) * (x >= data[i, 0]) 

sklejana1st += ((data[1, 1] - data[0, 1])/(data[1, 0] - data[0, 0]) * (x - data[0, 0]) + data[0, 1]) * (x <= data[0, 0]) * (x >= 0)
sklejana1st += ((data[-1, 1] - data[-2, 1])/(data[-1, 0] - data[-2, 0]) * (x - data[-2, 0]) + data[-2, 1]) * (x <= 10) * (x >= data[-1, 0])

'''
Sklejana 3 stopnia
'''

n = data.shape[0] - 1
h_i = numpy.zeros(n)
b_i = numpy.zeros(n)
u = numpy.zeros(n - 1)
v = numpy.zeros(n - 1)
z = numpy.zeros(n + 1)
a = numpy.zeros(n)
b = numpy.zeros(n)
c = numpy.zeros(n)

for i in range(n):
    h_i[i] = data[i+1, 0] - data[i, 0]
    b_i[i] = (6/h_i[i]) * (data[i+1, 1] - data[i, 1])


u1 = 2*(h_i[0] + h_i[1])

u[0] = u1

for i in range(1, n - 1):
    u[i] = 2*(h_i[i] + h_i[i+1]) - (h_i[i]**2)/u[i-1]

v1 = b_i[1] - b_i[0]
v[0] = v1

for i in range(1, n - 1):
    v[i] = b_i[i + 1] - b_i[i] - (h_i[i] * v[i-1])/u[i-1]

for i in range(n - 1, 0, -1):
    z[i] = (v[i-1] - (h_i[i] * z[i+1]))/u[i-1]

sklejana3st = numpy.zeros(x.shape[0])

for i in range(n):
    a[i] = (1/(6 * h_i[i])) * (z[i+1] - z[i])
    b[i] = z[i]/2
    c[i] = ((-1 * h_i[i])/6) * (z[i+1] + (2 * z[i])) + 1/h_i[i] * (data[i+1, 1] - data[i, 1])

for i in range(n):
    sklejana3st += (data[i, 1] + (x - data[i, 0]) * (c[i] + (x - data[i, 0]) * (b[i] + (x - data[i, 0]) * a[i]))) * ((x < data[i + 1, 0]) * (x >= data[i, 0]))

sklejana3st += (data[0, 1] + (x - data[0, 0]) * (c[0] + (x - data[0, 0]) * (b[0] + (x - data[0, 0]) * a[0]))) * ((x <= data[0, 0]) * (x >= 0))
sklejana3st += (data[-1, 1] + (x - data[-1, 0]) * (c[-1] + (x - data[-1, 0]) * (b[-1] + (x - data[-1, 0]) * a[-1]))) * ((x <= 10) * (x >= data[-1, 0]))



plt.title('Różne funkje Interpolacji')
plt.plot(x, lagrange, label="Interpolacja Lagrange'a")
plt.plot(x, sklejana1st, label='Interpolacja Sklejana 1 stopnia')
plt.plot(x, sklejana3st, label='Interpolacja Sklejana 3 stopnia')
plt.legend()
plt.show()



from scipy.interpolate import CubicSpline

spline = CubicSpline(data[:, 0], data[:, 1])
spline = spline(x)

plt.title('Porównanie CubicSpline i skl 3 stopnia')
plt.plot(x, sklejana3st, label='Sklejana 3 stopnia')
plt.plot(x, spline, label='Cubic Spline')
plt.legend()
plt.show()
