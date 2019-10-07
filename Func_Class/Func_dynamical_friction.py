#some function used for dynamical friction research procject with Prof.Newberg in 2017
#updated
#Write in Python 3 by JiaZhao Lin :)  9/17/2018

from Code.Func_Class.Point_Class import Point
import numpy as np
from scipy.stats import gaussian_kde
from numpy import sqrt, log, pi, e
from scipy.special import erf



def isotropic_v_dispersion(a, v0, r):
    return (v0 ** 2) / 2 * (4 * r * (a + r) ** 2 * log(1 + r / a) - a * r * (5 * a + 4 * r)) / (a ** 2 * (2 * a + r))


def anisotropic_v_dispersion(a, v0, r):
    return (v0 ** 2) / 2 * (3 * a + 2 * r) / (2 * a + r)


def dynamical_friction(point):
    k = 1.428
    r = point.r()
    Lambda = r / 1.6 / k
    a = 12
    v0 = 73

    sigma = isotropic_v_dispersion(a, v0, r)
    X = point.v() / sqrt(2 * sigma)

    density = 5329 / (2 * pi) * ((3 + r ** 2 / 144) / 144 / (1 + r ** 2 / 144) ** 2)

    F = -4 * pi * point.m ** 2 * log(Lambda) * density / point.v() ** 2 * (
                erf(X) - 2 * X / pi ** (1 / 2) * e ** (-X ** 2))
    return (F * point.vx / point.v(), F * point.vy / point.v(), F * point.vz / point.v(), abs(F))


def acceleration(point, DF_option = 'no'):
    """
    return list that contain accelerations in each axis caused by MW
    :param point:
    :param DF_option:
    :return:
    """
    FB = 1.52954402e5 / ((point.r() + 0.7) ** 2)
    FDR = 4.45865888e5 / ((point.R() ** 2 + (6.5 + (point.z ** 2 + 0.0676) ** 0.5) ** 2) ** (3 / 2)) * point.R()
    FDZ = 4.45865888e5 / ((point.R() ** 2 + (6.5 + (point.z ** 2 + 0.0676) ** 0.5) ** 2) ** (3 / 2)) * (
            6.5 + (point.z ** 2 + 0.0676) ** 0.5) / (point.z ** 2 + 0.0676) ** 0.5 * point.z
    FH = 5329 / (1 + point.r() ** 2 / 144) * point.r() * 2 / 144
    a = (FB + FH)
    ax = a * point.x / point.r() + FDR * point.x / point.R()
    ay = a * point.y / point.r() + FDR * point.y / point.R()
    az = a * point.z / point.r() + FDZ
    DF = dynamical_friction(point)

    if DF_option == 'yes':
        return (-ax + DF[0] / point.m, -ay + DF[1] / point.m, -az + DF[2] / point.m, DF[3])
    else:
        return (-ax + DF[0], -ay, -az, DF[3])


def Center(points, x, y, z, r, m):
    """
    function that create the center with center of mass velocity of
    all the mass within the radius r
    :param points: point object
    :param x: position x
    :param y: position y
    :param z: position z
    :param r: radius that we wish to include the mass
    :param m: mass of the center
    :return: point
    """
    p = (0, 0, 0)
    tm = 0
    fake = Point(x, y, z, 0, 0, 0, 0)
    for i in range(len(points)): #find the center of mass momentum
        if (Point.distance(points[i], fake) <= r):
            p[0] += points[i].p()[0]
            p[1] += points[i].p()[1]
            p[2] += points[i].p()[2]
            tm += points[i].m
        else:
            break
    return Point(x, y, z, p[0] / tm, p[1] / tm, p[2] / tm, m)


def mass(point, center, r):
    """
    sum all the mass with radius r, the points is sorted by distance to the center
    :param point: list of point
    :param center: center point
    :param r: radius
    :return: total mass
    """
    j = 0
    tm = 0
    points = point
    for i in range(len(points)):
        if Point.distance(points[i], center) <= r:
            j = i
        else:
            break
    for i in range(j + 1):
        tm += points[i].m
    return tm


def radius(points, center, saves = 5):
    """
    find radius by checking if a point is bounded to the center
    :param points: sorted list of point by distance to the center
    :param center: center point
    :param saves: number of encountering consecutive outlier before stop
    :return: raduis
    """
    tm = 0
    x = 0
    save = saves
    for i in range(1, len(points)):
        tm += points[i].m
        if Point.bound(points[i], center):
            save = saves
        elif save == 0:
            x = tm
            for j in range(i, 0, -1): #stop when save = 0
                x -= points[j].m
                if x <= tm * 0.98:  #take %98 of mass
                    return Point.distance(points[j], center)
        else:
            save -= 1
    x = tm
    for j in range(len(points) - 1, -1, -1):
        x -= points[j].m
        if x <= tm * 0.98:
            return Point.distance(points[j], center)


def data(file,colored = 'no'):
    """
    import data
    :param file: file name
    :param colored: choice if the output points needed to be colored,it will be much faster if not
    :return: if colored return (points,d,x,y,z); if not colored (points,x,y,z)
    """
    f = open(file)
    for i in range(1, 5):
        f.readline()
    i = f.readline()

    x, y, z, vx, vy, vz, m = [], [], [], [], [], [], []
    for i in f:
        l = i.split(",")
        x.append(float(l[1]))
        y.append(float(l[2]))
        z.append(float(l[3]))
        vx.append(float(l[7]))
        vy.append(float(l[8]))
        vz.append(float(l[9]))
        m.append(float(l[10]))

    points = []
    if colored:
        xyz = np.vstack([x, y, z])
        kde = gaussian_kde(xyz)
        density = kde(xyz)
        idx = density.argsort()
        d = []
        x1, y1, z1, vx1, vy1, vz1, m1, d = [], [], [], [], [], [], [], []
        for i in range(len(x)):
            j = idx[i]
            x1.append(x[j])
            y1.append(y[j])
            z1.append(z[j])
            vx1.append(vx[j])
            vy1.append(vy[j])
            vz1.append(vz[j])
            m1.append(m[j])
            d.append(density[j])
        for i in range(len(x)): #making points
            points.append(Point(x1[i], y1[i], z1[i], vx1[i], vy1[i], vz1[i], m1[i]))
        return points, d, x1, y1, z1
    else:
        for i in range(len(x)): #making points
            points.append(Point(x[i], y[i], z[i], vx[i], vy[i], vz[i], m[i]))
        return points, x, y, z
