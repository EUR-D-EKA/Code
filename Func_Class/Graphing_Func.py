import numpy as np
import matplotlib.pyplot as plt
from math import sqrt
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import gaussian_kde
from Code.Func_Class.Point_Class import Point


def plot(file, title = "", xlabel = "", ylabel=""):
    x = []
    for line in open(file):
        x.append(float(line))

    plt.figure()
    plt.plot(x)
    plt.xlabel('{}'.format(xlabel))
    plt.xticks([0,2000,4000,6000,8000],[0,1,2,3,4] )
    plt.ylabel('{}'.format(ylabel))
    plt.title('{}'.format(title))
    plt.show()

def parse_date(file):
    x, y, z, vx, vy, vz, m = [], [], [], [], [], [], []
    f = open(file)
    for i in range(1, 6):
        f.readline()
    f.readline()

    for i in f:
        l = i.split(",")
        x.append(float(l[2]))
        y.append(float(l[3]))
        z.append(float(l[4]))
        vx.append(float(l[8]))
        vy.append(float(l[9]))
        vz.append(float(l[10]))
        m.append(float(l[11]))
    return x, y, z, vx, vy, vz, m

def find_COM(file):
    parsed = parse_date(file)
    m1 = parsed[6][-1]

    P1 = Point(0,0,0,0,0,0,0)
    P2 = Point(0, 0, 0, 0, 0, 0, 0)
    for i in range(len(parsed[0])):
        if m1 != parsed[6][i]:
            P1.combine_point(Point(parsed[0][i],parsed[1][i],parsed[2][i],parsed[3][i],parsed[4][i],parsed[5][i],parsed[6][i]))
        else:
            P2.combine_point(Point(parsed[0][i],parsed[1][i],parsed[2][i],parsed[3][i],parsed[4][i],parsed[5][i],parsed[6][i]))

    print(P1)
    print(P2)

def scatter_plot(file, title = "", xlabel = "", ylabel="", axis1 = 0, axis2 = 1):
    parsed = parse_date(file)
    x, y = parsed[axis1], parsed[axis2]

    points = [ [-39.9749, 3.92602, 33.9475],[-33.9665, -3.18828, 31.178], [-25.8759, -8.33567, 25.7046],[ -21.8274, -10.8032, 23.7138],[-16.5659, -11.7899, 19.6945],[-12.8173, -12.8845, 16.3934 ],[-7.7442, -14.6548, 11.4513]]

    plt.figure()
    plt.scatter(x,y,s=1,c='k')
    for i in range(len(points)):
        plt.scatter(points[i][axis1],points[i][axis2],s=5,c='r')

    plt.xlabel('{}'.format(xlabel))
    plt.ylabel('{}'.format(ylabel))
    plt.title('{}'.format(title))
    # plt.xlim(-80, 60)
    # plt.ylim(-75, 125)
    plt.show()

def plot_heat_map(file):
    """
    plot the heat map of the point. the more red means more density
    :param file: name
    :return: graph
    """
    parsed = parse_date(file)
    x, y, z, vx, vy, vz, m = parsed

    xyz = np.vstack([x, y, z])
    kde = gaussian_kde(xyz)
    density = kde(xyz)

    idx = density.argsort()
    x1, y1, z1, vx1, vy1, vz1, m1, d = [], [], [], [], [], [], [], []
    for i in range(len(x)):
        j = idx[i]
        x1.append(x[j])
        y1.append(y[j])
        z1.append(z[j])
        vx1.append(vx[j])
        vy1.append(vy[j])
        vz1.append(vz[idx[i]])
        m1.append(m[j])
        d.append(density[j])

    print("Center is", x1[-1], y1[-1], z1[-1])
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x1, y1, z1, c=d, edgecolor='', s=1)
    ax.set_xlabel('Xcd nb Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    plt.show()


def plot_two_body(file, file1, skip = 1):
    """
    plot the position of two body
    :param file: name
    :param file1: name
    :param skip: skipping some particle to reduce computation
    :return: graph
    """
    f = open(file)
    f1 = open(file1)
    x, y, z, x1, y1, z1 = [], [], [], [], [], []
    for line in f:
        position = line.split()
        position1 = f1.readline().split()
        x.append(float(position[0]))
        y.append(float(position[1]))
        z.append(float(position[2]))
        x1.append(float(position1[0]))
        y1.append(float(position1[1]))
        z1.append(float(position1[2]))

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x[::skip], y[::skip], z[::skip], c='b', edgecolor='')
    ax.scatter(x1[::skip], y1[::skip], z1[::skip], c='r', edgecolor='')
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    plt.show()

def plot_three_body(file, file1, file2, skip = 1):
    """
    plot the position of three body
    :param file: name
    :param file1: name
    :param file2: name
    :param skip: skipping some particle to reduce computation
    :return: graph
    """
    f = open(file)
    f1 = open(file1)
    f2 = open(file2)
    x, y, z, x1, y1, z1, x2, y2, z2 = [], [], [], [], [], [], [], [], []
    for line in f:
        position = line.split()
        position1 = f1.readline().split()
        position2 = f2.readline().split()
        x.append(float(position[0]))
        y.append(float(position[1]))
        z.append(float(position[2]))
        x1.append(float(position1[0]))
        y1.append(float(position1[1]))
        z1.append(float(position1[2]))
        x2.append(float(position2[0]))
        y2.append(float(position2[1]))
        z2.append(float(position2[2]))

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x[::skip], y[::skip], z[::skip], c='b', edgecolor='')
    ax.scatter(x1[::skip], y1[::skip], z1[::skip], c='r', edgecolor='')
    ax.scatter(x2[::skip], y2[::skip], z2[::skip], c='g', edgecolor='')
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    plt.show()

def plot_hist(file):
    f = open(file)
    x = []
    y = []
    for i in f:
        if not i.startswith("1"):
            continue
        z = i.split()
        x.append(float(z[1]))
        y.append(float(z[3]))

    fig, axs = plt.subplots()
    axs.bar(x, y, align='center')
    plt.ylim(0, 0.35)
    plt.show()

def format_shift(file):
    x = []
    f = open(file)
    for i in f:
        x.append(i.split('\n')[0])
    f.close()

    x = x[::-1]
    f = open("../Output/shift.txt",'w')
    for i in x:
        f.write("{}\n".format(i))
    f.close()

if __name__ == '__main__':
    # print('yes')
    # plot_heat_map("../Output/some_output_LMC20000.out")
    # scatter_plot("../Output/7685_OG",'','x','y')
    # scatter_plot("../Output/7685_new_LMC_shift",'','x','y')
    # plot_hist("../Output/hist_noLMC_noshift")
    # plot_hist("../Output/path_to_input_hist")
    # format_shift("../Output/shift_r.txt")
    # find_COM('../Output/some_output_LMC20000.out')
    # find_COM('../Output/some_output_LMC1_fixed.out')
    # scatter_plot("../Output/some_output_LMC1.out", 'LMC has 1 particle ran for 3.8 Gyr', 'x', 'y',0,1)
    scatter_plot("../Output/some_output.out", 'LMC has 1 particle', 'x', 'z',0,2)
    # scatter_plot("../Output/7685_OG", 'Original orphan stream', 'z', 'y', 0, 1)