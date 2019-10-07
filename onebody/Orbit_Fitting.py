from Code.Func_Class.Point_Class import Point
from math import cos, sin, asin,pi,atan2
import matplotlib.pyplot as plt

def lbr_to_xyz(l, b, r):
    x = r* cos(l * pi/180)*cos(b* pi/180) - 8
    y = r* sin(l * pi/180)*cos(b* pi/180)
    z = r* sin(b * pi/180)
    return x,y,z

def xyz_to_lbr(x,y,z):
    xp = x + 8
    l = atan2(y , xp) * 180/pi
    b = atan2(z , (xp**2 + y**2)**0.5) * 180/pi
    if l < 0:
        l += 360
    r = (xp**2 + y**2 + z**2)**0.5
    return l,b,r

def v_los(x,y,z,vx,vy,vz):
    r = ((x+8)**2+y**2+z**2)**0.5
    x_ = (x+8) /r
    y_ = y / r
    z_ = z / r
    v_ = x_*vx + y_*vy + z_*vz
    return v_

def record(path,i,chi_2,Q):
    f = open(path, 'a')
    f.write("TimeStep:{}, Chi_2:{}, Q:{}\n".format(i, chi_2, Q))

def print_orbit(Q):
    x,y,z,vx,vy,vz = Model_cal(Q)
    l,b,R,rv = [],[],[],[]
    for i in range(len(x)):
        l1,b1,r1 = xyz_to_lbr(x[i],y[i],z[i])
        rv.append(v_los(x[i],y[i],z[i], vx[i],vy[i],vz[i]))
        l.append(l1)
        b.append(b1)
        R.append(((x[i]+8) ** 2 + y[i] ** 2 + z[i] ** 2) ** 0.5)

    fig, ax = plt.subplots()
    l_data = [173, 187, 205, 218, 234, 249.5, 271]
    b_data = [46.5, 50, 52.5, 53.5, 53.5, 50, 38]
    rv_data = [115.5, 119.7, 139.8, 131.5, 111.3, 101.4, 38.4]
    R_data = [46.8, 40.7, 32.4, 29.5, 24.4, 21.4, 18.6]

    plt.scatter(l, b, s=1, c='b', label="l vs b fit")
    plt.scatter(l, rv, s=1, c='g', label="l vs rv fit")
    plt.scatter(l, R, s=1, c='k', label="l vs R fit")
    plt.scatter(l_data, b_data, s=3, c='r')
    plt.scatter(l_data, rv_data, s=3, c='r')
    plt.scatter(l_data, R_data, s=3, c='r')

    plt.title("Moving MW and Massive LMC Fit")

    ax.legend()
    plt.show()

def print_fit(l,b,R,rv):
    l_data = [173, 187, 205, 218, 234, 249.5, 271]
    b_data = [46.5, 50, 52.5, 53.5, 53.5, 50, 38]
    rv_data = [115.5, 119.7, 139.8, 131.5, 111.3, 101.4, 38.4]
    R_data = [46.8, 40.7, 32.4, 29.5, 24.4, 21.4, 18.6]

    plt.scatter(l, b, s=1, c='b')
    plt.scatter(l, rv, s=1, c='g')
    plt.scatter(l, R, s=1, c='k')
    plt.scatter(l_data, b_data, s=5, c='r')
    plt.scatter(l_data, rv_data, s=5, c='r')
    plt.scatter(l_data, R_data, s=5, c='r')
    plt.show()

def one_bodies(x,y,z,vx,vy,vz):

    l_lim = 0

    mw_mass = 10 ** 12 / 222288.47
    mw_position_v = (0, 0, 0, 0, 0, 0)

    lmc_mass = 1 * 10 ** 11 / 222288.47
    lmc_position_v = (-1.1, -41.1, -27.9, -57, -226, 221)

    MW = Point(*mw_position_v, mw_mass)
    LMC = Point(*lmc_position_v, lmc_mass)

    t = 0.00051398136723478

    x1,y1,z1,vx1,vy1,vz1 = [],[],[],[],[],[]
    l1,b1,r1,rv1,R1 = [],[],[],[],[]
    orbit_point = Point(x,y,z,vx,vy,vz,0)

    while l_lim < 300:
        orbit_point.r_advance_vel(orbit_point.MW_acc(MW)[0]  ,\
                          orbit_point.MW_acc(MW)[1] ,\
                          orbit_point.MW_acc(MW)[2] , t / 2)
        orbit_point.r_advance_pos(t)
        orbit_point.r_advance_vel(orbit_point.MW_acc(MW)[0] ,\
                          orbit_point.MW_acc(MW)[1] ,\
                          orbit_point.MW_acc(MW)[2] , t / 2)

        LMC.r_advance_vel(*LMC.MW_acc(MW)[:-1], t / 2)
        LMC.r_advance_pos(t)
        LMC.r_advance_vel(*LMC.MW_acc(MW)[:-1], t / 2)

        MW.r_advance_vel(*MW.Point_acc(LMC), t / 2)
        # MW.r_advance_pos(t)
        MW.r_advance_vel(*MW.Point_acc(LMC), t / 2)
        x1.append(orbit_point.x)
        y1.append(orbit_point.y)
        z1.append(orbit_point.z)
        ll,bb,rr = xyz_to_lbr(orbit_point.x, orbit_point.y, orbit_point.z)
        rv1.append(v_los(orbit_point.x, orbit_point.y, orbit_point.z, orbit_point.vx, orbit_point.vy, orbit_point.vz))
        R1.append(((orbit_point.x+8)**2+ orbit_point.y**2+ orbit_point.z**2)**0.5)
        l1.append(ll)
        b1.append(bb)
        r1.append(rr)
        vx1.append(orbit_point.vx)
        vy1.append(orbit_point.vy)
        vz1.append(orbit_point.vz)

        l_lim = ll

    MW = Point(*mw_position_v, mw_mass)
    LMC = Point(*lmc_position_v, lmc_mass)
    t = -0.00051398136723478
    orbit_point = Point(x, y, z, vx, vy, vz, 0)

    l_lim = 200
    while l_lim > 170:
        orbit_point.r_advance_vel(orbit_point.MW_acc(MW)[0],\
                          orbit_point.MW_acc(MW)[1] ,\
                          orbit_point.MW_acc(MW)[2] , t / 2)
        orbit_point.r_advance_pos(t)
        orbit_point.r_advance_vel(orbit_point.MW_acc(MW)[0],\
                          orbit_point.MW_acc(MW)[1] ,\
                          orbit_point.MW_acc(MW)[2] , t / 2)

        LMC.r_advance_vel(*LMC.MW_acc(MW)[:-1], t / 2)
        LMC.r_advance_pos(t)
        LMC.r_advance_vel(*LMC.MW_acc(MW)[:-1], t / 2)

        MW.r_advance_vel(*MW.Point_acc(LMC), t / 2)
        # MW.r_advance_pos(t)
        MW.r_advance_vel(*MW.Point_acc(LMC), t / 2)
        x1.append(orbit_point.x)
        y1.append(orbit_point.y)
        z1.append(orbit_point.z)
        ll,bb,rr = xyz_to_lbr(orbit_point.x, orbit_point.y, orbit_point.z)
        rv1.append(v_los(orbit_point.x, orbit_point.y, orbit_point.z,orbit_point.vx, orbit_point.vy, orbit_point.vz))
        R1.append(((orbit_point.x+8) ** 2 + orbit_point.y ** 2 + orbit_point.z ** 2) ** 0.5)
        l1.append(ll)
        b1.append(bb)
        r1.append(rr)
        vx1.append(orbit_point.vx)
        vy1.append(orbit_point.vy)
        vz1.append(orbit_point.vz)

        l_lim = ll

    # print_fit(l1,b1,R1,rv1)

    return x1,y1,z1,vx1,vy1,vz1

def match_point(x_l,y_l,z_l,vx_l,vy_l,vz_l,l_0):
    for i in range(len(x_l)):
        l1,b1,r1 = xyz_to_lbr(x_l[i],y_l[i],z_l[i])
        # qq, ww, ee = lbr_to_xyz(l1, b1, r1)
        # print(qq)
        # print(l1)
        if abs(l1-l_0) < 0.4:
            return b1, v_los(x_l[i],y_l[i],z_l[i],vx_l[i],vy_l[i],vz_l[i]), ((x_l[i]+8)**2+y_l[i]**2+z_l[i]**2)**0.5
    # print(l_0)
    return 0

def Model_cal(Q):
    # I = (173, 46.5)
    I = (218, 53.5)
    x,y,z,vx,vy,vz = Q_to_poiont(Q,I)
    x_l, y_l, z_l, vx_l, vy_l, vz_l = one_bodies(x, y, z, vx, vy, vz)
    return x_l,y_l,z_l,vx_l,vy_l,vz_l

def Chi2(Q):

    l_data =  [173,   187,   205,   218,   234,   249.5, 271]
    b_data =  [46.5,  50,    52.5,  53.5,  53.5,  50,    38]
    rv_data = [115.5, 119.7, 139.8, 131.5, 111.3, 101.4, 38.4]
    R_data =  [46.8,  40.7,  32.4,  29.5,  24.4,  21.4,  18.6]
    b_model,rv_model,R_model = [], [], []
    sigma_b, sigma_rv, sigma_R = [0.7,1,0.7,1,0.7,0.7,3.5], [6.7,6.9,4.6,3.1,11.1,2.9,1.7], [4.5,1.9,1.5,1.4,1.2,1,0.9]

    chi2_b = 0
    chi2_rv = 0
    chi2_R = 0

    eta = len(b_data)*3-4-1

    x_l,y_l,z_l,vx_l,vy_l,vz_l = Model_cal(Q)

    for i in range(len(l_data)):
        b_, rv_, R_ = match_point(x_l,y_l,z_l,vx_l,vy_l,vz_l,l_data[i])
        b_model.append(b_)
        rv_model.append(rv_)
        R_model.append(R_)
    for i in range(len(b_data)):
        chi2_b += ((b_model[i] - b_data[i])/sigma_b[i])**2
        chi2_rv += ((rv_model[i] - rv_data[i])/sigma_rv[i])**2
        chi2_R += ((R_model[i] - R_data[i])/sigma_R[i])**2
    print(rv_model)
    print(rv_data)
    print("chi_b",chi2_b/eta,"chi_rv",chi2_rv,"chi_R",chi2_R/eta)
    chi2 = (chi2_b + chi2_rv + chi2_R)/eta

    return chi2

def Cal_Delta(Q,i, h):
    Qa = list(Q)
    Qb = list(Q)
    Qa[i] += h
    Qb[i] -= h
    return (Chi2(Qa) - Chi2(Qb))/2/h

def Gradien_descent(Q, Lambda):
    Q_new = []
    h = [0.1,1,1,1]

    for i in range(len(Q)):
        Q_new.append(Q[i] - h[i] * Lambda * Cal_Delta(Q,i, h[i]))

    return Q_new

def Q_to_poiont(Q,I):
    x, y, z = lbr_to_xyz(*I, Q[0])
    vx, vy, vz = Q[1:]

    return x,y,z,vx,vy,vz

def main():
    path = "../Output/Orbit_fitting_record.txt"

    # Q = [50.88341798196363, -217.47791991915517, 126.084990798251, 155.7335104594512]
    Q = [22.97220260367188, -180.90399959605722, 165.36251739264196, 78.54737409624518]
    Q=[28.6,-156,79,107]

    chi2_old = Chi2(Q)
    Lambda = 1
    for i in range(100):
        record(path,i,chi2_old,Q)
        print("timestep:",i,"chi_2",chi2_old,"Q ",Q,Lambda)
        Q = Gradien_descent(Q,Lambda)
        chi2_new = Chi2(Q)
        if chi2_new < chi2_old:
            Lambda *= 1.03
        else:
            Lambda *= 0.8
        chi2_old = chi2_new

    return Q




if __name__ == '__main__':
    # main()
    # print_orbit([27.636802120830918, -167.7393095375379, 103.63014042710199, 110.38789579909194])
    print(Chi2([28.6,-156,79,107]))
    print_orbit([28.6,-156,79,107])
    # print(Q_to_poiont([27.636802120830918, -167.7393095375379, 103.63014042710199, 110.38789579909194],(218, 53.5)))