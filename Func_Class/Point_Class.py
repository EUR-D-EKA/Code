#Point class used for dynamical friction research procject with Prof.Newberg in 2017
#updated With more stuff
#Write in Python 3 by JiaZhao Lin :)  9/17/2018

from numpy import sqrt, log, e, pi
from scipy.special import erf

class Point(object):
    """
    Point class represent a single body
    """
    def __init__(self, x0, y0, z0, vx0, vy0, vz0, m0):
        self.x = x0
        self.y = y0
        self.z = z0
        self.vx = vx0
        self.vy = vy0
        self.vz = vz0
        self.m = m0

    def __str__(self):
        s = "Center of mass Point Class Contains mass of {} Having position and velocity of\n{}".format(self.m,str((self.x,self.y,self.z,self.vx,self.vy,self.vz)))
        return s

    def r(self):
        """
        :return: distance from the center
        """
        return sqrt(self.x ** 2 + self.y ** 2 + self.z ** 2)

    def v(self):
        """
        :return: velocity relative to center
        """
        return sqrt(self.vx ** 2 + self.vy ** 2 + self.vz ** 2)

    def R(self):
        """
        :return:
        """
        return sqrt(self.x ** 2 + self.y ** 2)

    def p(self):
        """
        :return: momentum relative to center
        """
        return (self.m * self.vx, self.m * self.vy, self.m * self.vz)

    def abs_p(self):
        p1 = self.p()
        return sqrt(p1[0] ** 2 + p1[1] ** 2 + p1[2] ** 2)

    def RKE(self, o):
        """
        :param o: point object
        :return: relative kinetic energy
        """
        return 0.5 * self.m * ((self.vx - o.vx) ** 2 + (self.vy - o.vy) ** 2 + (self.vz - o.vz) ** 2)

    def RV(self, o):
        """
        :param o: point object
        :return: relative speed
        """
        return sqrt((self.vx - o.vx) ** 2 + (self.vy - o.vy) ** 2 + (self.vz - o.vz) ** 2)

    def distance(self, o):
        """
        :param o: point object
        :return: relative distance
        """
        return sqrt((self.x - o.x) ** 2 + (self.y - o.y) ** 2 + (self.z - o.z) ** 2)

    def bound(self, o):
        """
        check if a object is bounded to it by checking if total energy is less than 0
        :param o: point object
        :return: bool
        """
        return (self.RKE(o) - self.m * o.m / self.distance(o)) < 0

    def dist_from_com(self, o, PM = True):
        if PM:
            d = self.distance(o)
            return o.m * d / (self.m + o.m)
        else:
            return

    def combine_point(self,o):
        x = (self.x*self.m+o.x*o.m)/(self.m+o.m)
        y = (self.y * self.m + o.y * o.m) / (self.m + o.m)
        z = (self.z * self.m + o.z * o.m) / (self.m + o.m)
        vx = (self.vx*self.m+o.vx*o.m)/(self.m+o.m)
        vy = (self.vy * self.m + o.vy * o.m) / (self.m + o.m)
        vz = (self.vz * self.m + o.vz * o.m) / (self.m + o.m)
        m = (self.m + o.m)
        self.x = x
        self.y = y
        self.z = z
        self.vx = vx
        self.vy = vy
        self.vz = vz
        self.m = m


    def advance(self, ax, ay, az, time, m0):
        """
        classical method of advancing this point given the acceleration
        :param ax:
        :param ay:
        :param az:
        :param time:
        :param m0:
        :return:
        """
        self.x += self.vx * time + 0.5 * ax * time ** 2
        self.y += self.vy * time + 0.5 * ay * time ** 2
        self.z += self.vz * time + 0.5 * az * time ** 2
        self.vx += ax * time
        self.vy += ay * time
        self.vz += az * time
        self.m = m0

    def r_advance(self, ax, ay, az, time, m0):
        """
        classical method of reverse advancing this point given the acceleration
        :param ax:
        :param ay:
        :param az:
        :param time:
        :param m0:
        :return:
        """
        self.m = m0
        self.vx -= ax * time
        self.vy -= ay * time
        self.vz -= az * time
        self.x -= self.vx * time + 0.5 * ax * time ** 2
        self.y -= self.vy * time + 0.5 * ay * time ** 2
        self.z -= self.vz * time + 0.5 * az * time ** 2

    def probe_advance(self, ax, ay, az, time, m0):
        x = self.x + self.vx * time + 0.5 * ax * time ** 2
        y = self.y + self.vy * time + 0.5 * ay * time ** 2
        z = self.z + self.vz * time + 0.5 * az * time ** 2
        self.x = (self.x + x) / 2
        self.y = (self.y + y) / 2
        self.z = (self.z + z) / 2
        self.vx += ax * time
        self.vy += ay * time
        self.vz += az * time
        self.m = m0

    def probe_r_advance(self, ax, ay, az, time, m0):
        self.m = m0
        self.vx -= ax * time
        self.vy -= ay * time
        self.vz -= az * time
        x = self.x - self.vx * time - 0.5 * ax * time ** 2
        y = self.y - self.vy * time - 0.5 * ay * time ** 2
        z = self.z - self.vz * time - 0.5 * az * time ** 2
        self.x = (self.x + x) / 2
        self.y = (self.y + y) / 2
        self.z = (self.z + z) / 2

    def r_advance_vel(self, ax, ay, az, time):
        self.vx -= ax * time
        self.vy -= ay * time
        self.vz -= az * time

    def r_advance_pos(self, time):
        self.x -= self.vx * time
        self.y -= self.vy * time
        self.z -= self.vz * time

    def two_body_acc(self, o):
        """
        acceleration cause by another body
        :param o: another body
        :return: list of acceleration
        """
        d = self.distance(o)
        ax = o.m/d**2 * (self.x - o.x)/d
        ay = o.m/d**2 * (self.y - o.y)/d
        az = o.m/d**2 * (self.z - o.z)/d
        return -ax, -ay, -az

    def isotropic_v_dispersion(self, a, v0, r):
        return (v0 ** 2) / 2 * (4 * r * (a + r) ** 2 * log(1 + r / a) - a * r * (5 * a + 4 * r)) / (
                    a ** 2 * (2 * a + r))

    def anisotropic_v_dispersion(self, a, v0, r):
        return (v0 ** 2) / 2 * (3 * a + 2 * r) / (2 * a + r)

    def dynamical_friction(self, o , a, v0):
        r = self.distance(o)
        k = 1.428  #softening length
        Lambda = r / 1.6 / k  #Coulomb logarithm
        rv = self.RV(o)

        sigma = self.isotropic_v_dispersion(a, v0, r)  #one-dimensional velocity dispersion of the dark matter halo
        X = rv / sqrt(2 * sigma)

        density = 5329 / (2 * pi) * ((3 + r ** 2 / 144) / 144 / (1 + r ** 2 / 144) ** 2)

        F = -4 * pi * self.m ** 2 * log(Lambda) * density / rv ** 2 * (
                erf(X) - 2 * X / pi ** (1 / 2) * e ** (-X ** 2))
        return (F * (self.vx - o.vx) / rv, F * (self.vy - o.vy) / rv, F * (self.vz - o.vz) / rv, abs(F))

    def MiyamotoNagaiDiskAccel(self, b, c, R, x, y,z, MD):
        ADR = MD / ((R ** 2 + (b + (z ** 2 + c ** 2) ** 0.5) ** 2) ** (3 / 2))
        ADZ = MD / ((R ** 2 + (b + (z ** 2 + c ** 2) ** 0.5) ** 2) ** (3 / 2)) * \
              ( b + (z ** 2 + c ** 2) ** 0.5) / (z ** 2 + c ** 2) ** 0.5
        return -ADR * x, -ADR * y, -ADZ *z

    def LogHaloAccel(self, v0, a, x, y, z, gamma = 1):
        ALH = 2 * v0**2 * gamma**2 /(gamma**2 * (x**2 + y**2 + a**2) + z**2)
        return -ALH*x, -ALH*y, -ALH * z /gamma**2

    def SphericalAccel(self, d, x, y, z, r, MB):
        SA =  MB / (r * (r + d) ** 2)
        return -SA * x, -SA * y , -SA * z

    def MW_acc(self, o, DF_option = False, NFW_halo = False, scale_v = False):
        """
        acceleration cause by three component mw potential
        :param o: mw body
        :return: acceleration
        """
        a = 12  #scale parameter
        b = 6.5  #scale parameters
        c = 0.26  #scale parameters
        d = 0.7  #scale length
        MB = 1.52954402e5  #mass of bulge 3.4*10**10
        MD = 4.45865888e5  #mass of disk 10*11
        v0 = 74.61  # scale velocity

        if NFW_halo: #update the params if we are using NFW halo
            b = 3.5
            c = 0.35
            d = 0.7
            MD = 5.5*10**10 /222288.47
            MB = 10**10 /222288.47


        #----------------------Accerleration--------------------------------
        #relative position
        r = self.distance(o)
        x = self.x - o.x
        y = self.y - o.y
        z = self.z - o.z
        R = sqrt(x ** 2 + y ** 2)

        #acc from bulge
        AB = self.SphericalAccel(d,x,y,z,r,MB)
        #acc from disk
        AD = self.MiyamotoNagaiDiskAccel(b,c,R,x,y,z,MD)
        #acc from halo
        if NFW_halo:

            c1 = 2999
            rs = 31.27
            M_viral = 1.5 * 10 ** 12 / 222288.47

            H = 80000
            rho = 3*H**2/8/pi
            r200 = (M_viral/200/rho/4/pi*3)**(1/3)
            c1 = rs/r200

            cons = log(1+c1) - c1/(1+c1)
            AH = M_viral/r**2/cons*log(1+r/rs) + M_viral/r/cons/(1+r/rs)/rs

        else:
            if scale_v: #you want to change the scale velocity?
                v0 = ((o.m / r ** 2 - sqrt(AD[0]**2+AD[1]**2 + AD[2]**2) - AB) * ((1 + r ** 2 / a**2) * r * 2 / a**2))**0.5

            AH = self.LogHaloAccel(v0,a,x,y,z)
        #--------------------END---------------------------------------------

        DF = self.dynamical_friction(o, a, v0)

        acc_x = AB[0] + AD[0] + AH[0]
        acc_y = AB[1] + AD[1] + AH[1]
        acc_z = AB[2] + AD[2] + AH[2]

        if DF_option:
            return (acc_x + DF[0] / self.m, acc_y + DF[1] / self.m, acc_z + DF[2] / self.m, DF[3])
        else:
            return (acc_x, acc_y, acc_z, DF[3])

    def Plummer_acc(self, o, a=15):
        """
        acceleration cause by Plummer potential
        :param o: body with Plummer potential
        :param a: scale radius
        :return: acceleration
        """

        r = self.distance(o)

        A = o.m * r / (r ** 2 + a**2)**(3/2)

        return -A * (self.x - o.x)/r, -A * (self.y - o.y)/r, -A * (self.z - o.z)/r

    def Point_acc(self, o):

        r = self.distance(o)

        return - o.m * (self.x - o.x)/r**3, - o.m * (self.y - o.y)/r**3, - o.m * (self.z - o.z)/r**3

    def forward_simulation(self,o,t,file):
        f = open(file)
        for i in f:
            x=i.split('\n')[0].split(' ')
            ax,ay,az,ax1,ay1,az1=float(x[0]),float(x[1]),\
                                 float(x[2]),float(x[3]),\
                                 float(x[4]),float(x[5])
            a,b,c,d = o.MW_acc(self)
            o.r_advance_vel(a-ax1, b-ay1, c-az1, -t / 2)
            o.r_advance_pos(-t)
            # o.x -= 0.5 *ax1*t**2
            # o.y -= 0.5 * ay1 * t ** 2
            # o.z -= 0.5 * az1 * t ** 2
            a, b, c, d = o.MW_acc(self)
            o.r_advance_vel(a-ax,b-ay,c-az, -t / 2)
        print('forward simulated',o.x,o.y,o.z,o.vx,o.vy,o.vz)
