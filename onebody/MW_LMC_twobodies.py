from Code.Func_Class.Point_Class import Point
from Code.Func_Class.Graphing_Func import *
import copy



#-------------------------set prarms---------------------------------------
path = "../Output/mw_position_velocity"
path1 = "../Output/lmc_position_velocity"
path2 = "../Output/shift_r.txt"
path3 = "../Output/dw_position_velocity"
f = open(path,'w')
f1 = open(path1,'w')
f2 = open(path2,'w')
f3 = open(path3,'w')

mw_reverse = True
probe = False

r_Gyr = 4  #number of Gyr to reverse
step = 10000

mw_mass = 10**12 / 222288.47
mw_position_v = (0,0,0,0,0,0)

lmc_mass = 1*10**11 / 222288.47
lmc_position_v = (-1.1,-41.1,-27.9,-57,-226,221)

dwarf_mass = 96
# dwarf_position_v = (-21.405585, -10.473591, 22.990306, -156, 79, 107)
dwarf_position_v = (-20.954108506210716, -10.120858780110828, 22.216032990350396, -167.7393095375379, 103.63014042710199, 110.38789579909194)
#---------------------------------------------------------------------------


MW = Point(*mw_position_v,mw_mass)
LMC = Point(*lmc_position_v,lmc_mass)
DW = Point(*dwarf_position_v,dwarf_mass)
f.write("{} {} {}\n".format(MW.x, MW.y, MW.z))
f1.write("{} {} {}\n".format(LMC.x, LMC.y, LMC.z))
f3.write("{} {} {}\n".format(DW.x, DW.y, DW.z))

t = 1 / step
t = 0.00051398136723478
closest = DW.distance(LMC)
print(closest)
for i in range(7685):
    ox,oy,oz,ovx,ovy,ovz = MW.x, MW.y, MW.z, MW.vx, MW.vy, MW.vz
    mw_ax,mw_ay,mw_az = MW.Point_acc(LMC)
    lmc_ax, lmc_ay, lmc_az , DF = LMC.MW_acc(MW)
    dw_ax, dw_ay, dw_az = DW.MW_acc(MW)[0] + DW.Point_acc(LMC)[0],\
                          DW.MW_acc(MW)[1] + DW.Point_acc(LMC)[1],\
                          DW.MW_acc(MW)[2] + DW.Point_acc(LMC)[2]

#--------------------probing-------------------------
    if probe:
        f_MW = copy.copy(MW)
        f_LMC = copy.copy(LMC)
        f_MW.probe_r_advance(mw_ax,mw_ay,mw_az,t/2,MW.m)
        f_LMC.probe_r_advance(lmc_ax, lmc_ay, lmc_az,t/2,LMC.m)
        mw_ax, mw_ay, mw_az = f_MW.Point_acc(f_LMC)
        lmc_ax, lmc_ay, lmc_az, DF = f_LMC.MW_acc(f_MW)
#----------------------------------------------------

#---------------------MW reverse---------------------
    if mw_reverse:
        LMC.r_advance_vel(lmc_ax, lmc_ay, lmc_az, t/2)
        LMC.r_advance_pos(t)
        f1.write("{} {} {} {} {} {}\n".format(LMC.x, LMC.y, LMC.z, LMC.vx, LMC.vy, LMC.vz))
        LMC.r_advance_vel(*LMC.MW_acc(MW)[:-1],t/2)

        DW.r_advance_vel(dw_ax, dw_ay, dw_az, t / 2)
        DW.r_advance_pos(t)
        f3.write("{} {} {} {} {} {}\n".format(DW.x, DW.y, DW.z, DW.vx, DW.vy, DW.vz))
        DW.r_advance_vel(DW.MW_acc(MW)[0] + DW.Point_acc(LMC)[0],\
                          DW.MW_acc(MW)[1] + DW.Point_acc(LMC)[1],\
                          DW.MW_acc(MW)[2] + DW.Point_acc(LMC)[2], t / 2)

        MW.r_advance_vel(mw_ax,mw_ay,mw_az,t/2)
        MW.r_advance_pos(t)
        f2.write("{} {} {} ".format(mw_ax, mw_ay, mw_az))
        mw_ax, mw_ay, mw_az = MW.Point_acc(LMC)
        f2.write("{} {} {}\n".format(mw_ax, mw_ay, mw_az))
        MW.r_advance_vel(mw_ax, mw_ay, mw_az, t/2)
        f.write("{} {} {} {} {} {}\n".format(MW.x, MW.y, MW.z, MW.vx, MW.vy, MW.vz))
        if closest > DW.distance(LMC):
            closest = DW.distance(LMC)

#----------------------------------------------------
    else:
        MW.r_advance(mw_ax, mw_ay, mw_az, t, MW.m)
        LMC.r_advance(lmc_ax, lmc_ay, lmc_az, t, LMC.m)
        DW.r_advance(dw_ax, dw_ay, dw_az, t, DW.m)
        f.write("{} {} {} {} {} {}\n".format(MW.x, MW.y, MW.z, MW.vx, MW.vy, MW.vz))
        f1.write("{} {} {} {} {} {}\n".format(LMC.x, LMC.y, LMC.z, LMC.vx, LMC.vy, LMC.vz))
        f2.write("{} {} {}\n".format(mw_ax, mw_ay, mw_az))
        f3.write("{} {} {} {} {} {}\n".format(DW.x, DW.y, DW.z, DW.vx, DW.vy, DW.vz))


# for i in range(8000):
#     mw_ax,mw_ay,mw_az = MW.two_body_acc(LMC)
#     lmc_ax, lmc_ay, lmc_az = LMC.two_body_acc(MW)
#     MW.advance(mw_ax,mw_ay,mw_az,t,MW.m)
#     LMC.advance(lmc_ax, lmc_ay, lmc_az,t,LMC.m)

print(MW.x,MW.y,MW.z,MW.vx,MW.vy,MW.vz)
print(LMC.x,LMC.y,LMC.z,LMC.vx,LMC.vy,LMC.vz)
print(DW.x, DW.y, DW.z, DW.vx, DW.vy, DW.vz)
print("LMC converted",LMC.x-MW.x,LMC.y-MW.y,LMC.z-MW.z,LMC.vx-MW.vx,LMC.vy-MW.vy,LMC.vz-MW.vz)
print("DW converted",DW.x-MW.x,DW.y-MW.y,DW.z-MW.z,DW.vx-MW.vx,DW.vy-MW.vy,DW.vz-MW.vz)
f.close()
f1.close()
f2.close()
f3.close()
format_shift(path2)
print(closest)

plot_three_body(path,path1,path3,100)
#plot(path2)


fake_MW = Point(*mw_position_v,mw_mass)
fake_LMC = Point(LMC.x-MW.x,LMC.y-MW.y,LMC.z-MW.z,LMC.vx-MW.vx,LMC.vy-MW.vy,LMC.vz-MW.vz,lmc_mass)
fake_MW.forward_simulation(fake_LMC,t,"../Output/shift.txt")