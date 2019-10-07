from timeit import default_timer as timer
from Code.Func_Class.Func_dynamical_friction import *


r = 0
m = 96
totolm = [96]
time = []
centerx, centery, centerz = [], [], []
for i in range(0, 7680):
    start = timer()
    data1 = data(str(i))
    points = data1[0]
    center = Center(points, points[0].x, points[0].y, points[0].z, r, m)
    centerx.append(center.x)
    centery.append(center.y)
    centerz.append(center.z)
    r = radius(points, center)
    m = m
    totolm.append(m)
    time.append(i)
    end = timer()
    print("in progress!! ", i, " time: ", end - start)

f = open("mass_file.txt", "w")
f2 = open("center_file", "w")

print('Center is ', center.x, ' ', center.y, ' ', center.y, ' with radius ', r)
for i in range(len(time)):
    f.write(str(totolm[i]) + "\n")
    f2.write("{} {} {}\n".format(centerx[i], centery[i], centerz[i]))
    print("In timestep {} the center has mass is {}".format(time[i], totolm[i]))
