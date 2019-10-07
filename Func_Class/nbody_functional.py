import os
import subprocess
from subprocess import call
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.patches as mpatches

# # # # # # # # # # # # # # # # # # # # # #
#        USEFUL CLASSES                   #
# # # # # # # # # # # # # # # # # # # # # #
class nbody_outputs:#a class that takes in data from nbody output files and makes them available
    def __init__(self, file_name):
        self.file_name = file_name
        self.get_data()
        
    def get_data(self):
        self.xs = []; self.ys = []; self.zs = []
        self.ls = []; self.bs = []; self.rs = []
        self.vxs = []; self.vys = []; self.vzs = []
        self.vls = []; self.ms = []; self.tps = []
        
        f = open(self.file_name, 'r')
        read_data = False
        
        for line in f:
            if (line.startswith("# ignore")):
                read_data = True
                continue
            if(line.startswith("</bodies>")):
                break
            if(read_data):
                ss = line.split(', ')
                ty = int(ss[0])
                x = float(ss[1])
                y = float(ss[2])
                z = float(ss[3])
                
                l = float(ss[4])
                b = float(ss[5])
                r = float(ss[6])
                
                vx = float(ss[7])
                vy = float(ss[8])
                vz = float(ss[9])
                
                m  = float(ss[10])
                vl = float(ss[11])
                self.xs.append(x); self.ys.append(y); self.zs.append(z)
                self.ls.append(l); self.bs.append(b); self.rs.append(r)
                self.vxs.append(vx); self.vys.append(vy); self.vzs.append(vz)
                self.tps.append(ty); self.ms.append(m); self.vls.append(vl)
            
        f.close()
        
    def dark_light_split(self):#splits the data between light and dark
        self.light_x , self.light_y , self.light_z    = ([] for i in range(3))
        self.light_l , self.light_b , self.light_r    = ([] for i in range(3))
        self.light_vx , self.light_vy , self.light_vz = ([] for i in range(3))
        self.light_vl, self.light_m                   = ([] for i in range(2))
        
        self.dark_x , self.dark_y , self.dark_z       = ([] for i in range(3))
        self.dark_l , self.dark_b , self.dark_r       = ([] for i in range(3))
        self.dark_vx , self.dark_vy , self.dark_vz    = ([] for i in range(3))
        self.dark_vl, self.dark_m                     = ([] for i in range(2))
        
        for i in range(0, len(self.xs)):
            if(self.tps[i] == 0):
                self.light_x.append(self.xs[i])
                self.light_y.append(self.ys[i])
                self.light_z.append(self.zs[i])
                
                self.light_l.append(self.ls[i])
                self.light_b.append(self.bs[i])
                self.light_r.append(self.rs[i])
                
                self.light_vx.append(self.vxs[i])
                self.light_vy.append(self.vys[i])
                self.light_vz.append(self.vzs[i])
                
                self.light_vl.append(self.vls[i])
                self.light_m.append(self.ms[i])
                
            if(self.tps[i] == 1):
                self.dark_x.append(self.xs[i])
                self.dark_y.append(self.ys[i])
                self.dark_z.append(self.zs[i])
                
                self.dark_l.append(self.ls[i])
                self.dark_b.append(self.bs[i])
                self.dark_r.append(self.rs[i])
                
                self.dark_vx.append(self.vxs[i])
                self.dark_vy.append(self.vys[i])
                self.dark_vz.append(self.vzs[i])
                
                self.dark_vl.append(self.vls[i])
                self.dark_m.append(self.ms[i])    
                
    def rescale_l(self):#to change l range from [0:360] to [-180:180]
        for i in range(0, len(self.ls)):
            if(self.ls[i] > 180.0):
                self.ls[i] = self.ls[i] - 360.0
    
    def convert_lambda_beta(self, split):#to convert l,b to lambda, beta
        if(split):#if the data was split between light and dark
            self.light_lambda = []; self.light_beta = []
            self.dark_lambda  = []; self.dark_beta = []
        else:
            self.betas   = []
            self.lambdas = []
        
        for i in range(0, len(self.ls)):
            lmbda_tmp, beta_tmp = self.convert_to_Lambda_Beta(self.ls[i], self.bs[i], self.rs[i], False)
            if(split):
                if(self.tps[i] == 0):
                    self.light_lambda.append(lmbda_tmp)
                    self.light_beta.append(beta_tmp)
                    
                if(self.tps[i] == 1):
                    self.dark_lambda.append(lmbda_tmp)
                    self.dark_beta.append(beta_tmp)
            else:
                self.lambdas.append(lmbda_tmp)
                self.betas.append(beta_tmp)
            
        
        
    def convert_to_Lambda_Beta(self, x1, x2, x3, cartesian):#can convert l,b or x,y,z to lambda beta
        phi   = mt.radians(128.79)
        theta = mt.radians(54.39)
        psi   = mt.radians(90.70)
        
        if(cartesian):
            x_coor = x1
            y_coor = x2
            z_coor = x3
            x_coor += 8.0 #convert to solar centric
        else:
            l = mt.radians(x1)
            b = mt.radians(x2)
            r = x3
            
            x_coor = r * mt.cos(l) * mt.cos(b) #this is solar centered x
            y_coor = r * mt.sin(l) * mt.cos(b)
            z_coor = r * mt.sin(b)
        
        #A = MB
        B = [x_coor, y_coor, z_coor]
        M_row1 = [mt.cos(psi) * mt.cos(phi) - mt.cos(theta) * mt.sin(phi) * mt.sin(psi),
                mt.cos(psi) * mt.sin(phi) + mt.cos(theta) * mt.cos(phi) * mt.sin(psi),
                mt.sin(psi) * mt.sin(theta)]
        
        M_row2 = [-mt.sin(psi) * mt.cos(phi) - mt.cos(theta) * mt.sin(phi) * mt.cos(psi),
                -mt.sin(psi) * mt.sin(phi) + mt.cos(theta) * mt.cos(phi) * mt.cos(psi),
                mt.cos(psi) * mt.sin(theta)]
        
        M_row3 = [mt.sin(theta) * mt.sin(phi), 
                -mt.sin(theta) * mt.cos(phi),
                mt.cos(theta)]
        
        A1 = M_row1[0] * B[0] + M_row1[1] * B[1] + M_row1[2] * B[2]
        A2 = M_row2[0] * B[0] + M_row2[1] * B[1] + M_row2[2] * B[2]
        A3 = M_row3[0] * B[0] + M_row3[1] * B[1] + M_row3[2] * B[2]
        
        beta = mt.asin(-A3 / mt.sqrt(A1 * A1 + A2 * A2 + A3 * A3))
        lamb = mt.atan2(A2, A1)
        
        beta = mt.degrees(beta)
        lamb = mt.degrees(lamb)
        
        return lamb, beta
    
    
    def binner_vlos(self, angle_cuttoffs):
        self.binned_vlos = []
        self.which_bin   = []
        
        bin_size = abs(angle_cuttoffs[0] - angle_cuttoffs[1]) / angle_cuttoffs[2]
        mid_bins = []
        which_lambda = []
        which_beta = []
        #setting up the mid bin coordinates
        for i in range(0, angle_cuttoffs[2]):
            mid_bin = angle_cuttoffs[0] + i * bin_size + bin_size / 2.0
            mid_bins.append(mid_bin)
            
        #transform to lambda beta coordinates from lbr
        self.convert_lambda_beta(False)
        
        for i in range(0, len(self.lambdas)):#go through the different lambda coordinates
            if(self.betas[i] >= angle_cuttoffs[3] and self.betas[i] <= angle_cuttoffs[4]):#if it is between the beta cuttoffs
                for j in range(0, len(mid_bins)):#go through the bin coordinates
                    left_edge  = mid_bins[j] - bin_size / 2.0 #edges of the bin
                    right_edge = mid_bins[j] + bin_size / 2.0 #edges of the bin
                    
                    if(self.lambdas[i] >= left_edge and self.lambdas[i] <= right_edge):#check if the lambda coor falls in the bin
                        self.which_bin.append(mid_bins[j])#which mid bin it should be 
                        self.binned_vlos.append(self.vls[i])#the line of sight vel
                        
                        which_lambda.append(self.lambdas[i])#the coordinate that was put there
                        which_beta.append(self.betas[i])
                        break 

class nbody_histograms:#a class that takes in data from nbody histogram files and makes them available
    def __init__(self, file_name):
        self.file_name = file_name
        self.get_data()
        
    def get_data(self):
        self.lbins = []; self.counts = []; self.count_err = []; self.vd = []; self.vd_error = []
        read_data = False

        lines = open(self.file_name, 'r')
        for line in lines:
            if (line.startswith("betaBins")):
                read_data = True
                continue
            if(line.startswith("</histogram>")):
                break
            
            if(read_data):
                if(line.startswith("\n")):
                    continue
                else:
                    ss = line.split(' ')
                    self.lbins.append(    float(ss[1]))
                    self.counts.append(   float(ss[3]))
                    self.count_err.append(float(ss[4]))
                    self.vd.append(       float(ss[5]))
                    self.vd_error.append( float(ss[6]))
                    
        lines.close()
        

        
class nbody_running_env:
    def __init__(self, lua_file, version, path):
        self.lua_file      = lua_file
        self.version       = version
        self.path          = path
    
    def build(self):
        os.chdir("./")
        #-DCMAKE_C_COMPILER=/usr/bin/cc 
        #os.system("rm -r nbody_test")
        #os.system("mkdir nbody_test")
        os.chdir("nbody_test")
        os.system("cmake -DCMAKE_BUILD_TYPE=Release -DNBODY_DEV_OPTIONS=ON -DBOINC_RELEASE_NAMES=OFF -DNBODY_GL=ON -DNBODY_STATIC=ON -DBOINC_APPLICATION=OFF -DSEPARATION=OFF -DNBODY_OPENMP=ON    " + self.path + "milkywayathome_client/")
        os.system("make -j ")
        os.chdir("../")
    
    
    def single_run(self, parameters, simulation_hist):
        ft    = str(parameters[0])
        bt    = str(parameters[1])
        rl    = str(parameters[2])
        rr    = str(parameters[3])
        ml    = str(parameters[4])
        mr    = str(parameters[5])
        print('running nbody')
        os.chdir("nbody_test/bin/")
        os.system("./milkyway_nbody" + self.version + " \
            -f " + self.path + "lua/" + self.lua_file + " \
            -z " + self.path + "quick_plots/hists/" + simulation_hist + ".hist \
            -o " + self.path + "quick_plots/outputs/" + simulation_hist + ".out \
            -n 10 -b -u -i " + (ft) + " " + bt + " " + rl + " " + rr + " " + ml + " " + mr)
        
        #os.system("./milkyway_nbody" + self.version + " \
            #-f " + self.path + "lua/" + self.lua_file + " \
            #-z " + self.path + "quick_plots/hists/" + simulation_hist + ".hist \
            #-o " + self.path + "quick_plots/outputs/" + simulation_hist + ".out \
            #-n 10 -b -u -i 100.0 ~/Desktop/research/test2.out" )
        
    
    def run_and_compare(self, parameters, correctans_hist, comparison_hist):
        ft    = str(parameters[0])
        bt    = str(parameters[1])
        rl    = str(parameters[2])
        rr    = str(parameters[3])
        ml    = str(parameters[4])
        mr    = str(parameters[5])
        print('running nbody 2')
        os.system(" " + self.path + "nbody_test/bin/milkyway_nbody" + self.version + " \
            -f " + self.path + "lua/" + self.lua_file + " \
            -h " + self.path + "quick_plots/hists/" + correctans_hist + ".hist \
            -z " + self.path + "quick_plots/hists/" + comparison_hist + ".hist \
            -o " + self.path + "quick_plots/outputs/" + comparison_hist + ".out \
            -n 10 -b -P -i  " + (ft) + " " + bt + " " + rl + " " + rr + " " + ml + " " + mr )
        
        
    def match_hists(self, hist1, hist2):
        print "matching histograms: "
        #using call here instead so the format of using it is on record
        call([" " + self.path + "nbody_test/bin/milkyway_nbody" + self.version  
            + " -h " + self.path + "quick_plots/hists/" + hist1 + '.hist'
            + " -S " + self.path + "quick_plots/hists/" + hist2 + '.hist'], shell=True)
        print hist1, "\n", hist2
        print "\n"
        return 0
    
    
    
    def match_hists_pipe_output(hist1, hist2, pipe_name):
        print "matching histograms: "
        #using call here instead so the format of using it is on record
        call([" " + self.path + "nbody_test/bin/milkyway_nbody" + self.version  
            + " -h " + self.path + "" + hist1 + '.hist'
            + " -S " + self.path + "" + hist2 + '.hist' + " 2>>" + pipe_name], shell=True)
        print hist1, "\n", hist2
        print "\n"
        return 0
    

    def single_piping_output(self, parameters, simulation_hist, output):
        ft    = str(parameters[0])
        bt    = str(parameters[1])
        rl    = str(parameters[2])
        rr    = str(parameters[3])
        ml    = str(parameters[4])
        mr    = str(parameters[5])
        print('running nbody. piping output to file')
        os.chdir("nbody_test/bin/")
        os.system("./milkyway_nbody" + self.version + " \
            -f " + self.path + "lua/" + self.lua_file + " \
            -z " + self.path + "quick_plots/hists/" + simulation_hist + ".hist \
            -o " + self.path + "quick_plots/outputs/" + output + ".out \
            -n 12 -b -P  -i " + (ft) + " " + bt + " " + rl + " " + rr + " " + ml + " " + mr + " \
        2>> " + out + "_piped.out")
        
        
    def run_from_checkpoint(self, parameters, simulation_hist, output):
        ft    = str(parameters[0])
        bt    = str(parameters[1])
        rl    = str(parameters[2])
        rr    = str(parameters[3])
        ml    = str(parameters[4])
        mr    = str(parameters[5])
        print('running nbody from checkpoint')
        os.chdir("nbody_test/bin/")
        os.system("./milkyway_nbody" + self.version + " \
            -f " + path + "lua/" + self.lua_file + " \
            -z " + path + "quick_plots/hists/" + simulation_hist + ".hist \
            -o " + path + "quick_plots/outputs/" + output + ".out \
            -n 10 -b  -P --no-clean-checkpoint --checkpoint=CHECKPOINT_NAME " + (ft) + " " + bt + " " + rl + " " + rr + " " + ml + " " + mr)
        
        
    
    def run_and_compare_from_checkpoint(self, parameters, comparison_hist, simulation_hist, output):
        ft    = str(parameters[0])
        bt    = str(parameters[1])
        rl    = str(parameters[2])
        rr    = str(parameters[3])
        ml    = str(parameters[4])
        mr    = str(parameters[5])
        print('running nbody 2')
        os.system(" " + self.path + "nbody_test/bin/milkyway_nbody_1.66_x86_64-pc-linux-gnu__mt" + self.version + " \
            -f " + self.path + "lua/" + lua_file + " \
            -h " + self.path + "quick_plots/hists/" + comparison_hist + ".hist \
            -z " + self.path + "quick_plots/hists/" + simulation_hist + ".hist \
            -o " + self.path + "quick_plots/outputs/" + output + ".out \
            -n 10 -b -P -i --no-clean-checkpoint --checkpoint=CHECKPOINT_NAME " + (ft) + " " + bt + " " + rl + " " + rr + " " + ml + " " + mr )

     
        
        #os.chdir("../")
    # #     



# # # # # # # # # # # # # # # # # # # # # #
#        histogram plot                   #
# # # # # # # # # # # # # # # # # # # # # #
# # 

def plot(hist1, hist2, name, label1, label2):
    ylimit = 1.0
    xlower = 180 
    xupper = -180
    w_overlap = 2.5
    w_adjacent = 1.5
    folder = 'quick_plots/hists/'
    #folder = ''
    #folder = 'like_surface/'
    save_folder_ove = 'quick_plots/comp_hist_plots/overlap/'
    save_folder_adj = 'quick_plots/comp_hist_plots/adj/'
    #os.system("" + path + "scripts/plot_matching_hist.py " + hist1 + " " + hist2)
    print "plot histogram 1: ", hist1
    print "plot histogram 2: ", hist2
    plot_hist1 = hist1 + ".hist"
    plot_hist2 = hist2 + ".hist"

    
    print("plotting histograms\n")
    hist1 = nbody_histograms(folder + plot_hist1)
    hist2 = nbody_histograms(folder + plot_hist2)
            
            
    # plot overlapping #
    #f, (f1, f2) = plt.subplots(2, sharex = True, sharey = True)
    #plt.subplot(211)
    plt.bar(hist1.lbins, hist1.counts, width = w_overlap, color='k', alpha=1,    label= label1)
    plt.bar(hist2.lbins, hist2.counts, width = w_overlap, color='r', alpha=0.75, label= label2)
    plt.title('Histogram of Light Matter Distribution After 4 Gy')
    plt.xlim((xlower, xupper))
    plt.ylim((0.0, ylimit))
    plt.ylabel('counts')
    plt.xlabel('Lambda')
    plt.legend()
    plt.savefig(save_folder_ove + name + '_overlapping.png', format='png')
    plt.clf()
    #plt.show()
        
    # plot_adjacent #
    plt.subplot(211)
    #f, (f1, f2) = plt.subplots(2, sharex = True, sharey = True)
    plt.bar(hist1.lbins, hist1.counts, width = w_adjacent, color='b')
    plt.legend(handles=[mpatches.Patch(color='b', label= plot_hist1)])
    plt.title('Histogram of Light Matter Distribution After 4 Gy')
    plt.xlim((xlower, xupper))
    plt.ylim((0.0, ylimit))
    plt.ylabel('counts')
    plt.xlabel('Lambda')

    plt.subplot(212)
    plt.bar(hist2.lbins, hist2.counts, width = w_adjacent, color='k')
    plt.legend(handles=[mpatches.Patch(color='k', label= plot_hist2)])
    plt.xlim((xlower, xupper))
    plt.ylim((0.0, ylimit))
    plt.xlabel('l')
    plt.ylabel('counts')
    #f.subplots_adjust(hspace=0)
    plt.savefig(save_folder_adj + name + '.png', format='png')
    plt.clf()
    #plt.show()
    return 1
# # 
def plot_veldisp(hist1, hist2, name, label1, label2):
    ylimit = 100
    xlower = 180 
    xupper = -180
    w_overlap = 2.5
    w_adjacent = 1.5
    folder = 'quick_plots/hists/'
    #folder = 'like_surface/'
    save_folder_ove = 'quick_plots/comp_hist_plots/overlap/'
    save_folder_adj = 'quick_plots/comp_hist_plots/adj/'
    #os.system("" + path + "scripts/plot_matching_hist.py " + hist1 + " " + hist2)
    print "plot histogram 1: ", hist1
    print "plot histogram 2: ", hist2
    plot_hist1 = hist1 + ".hist"
    plot_hist2 = hist2 + ".hist"

    print("plotting histograms\n")
    hist1 = nbody_histograms(folder + plot_hist1)
    hist2 = nbody_histograms(folder + plot_hist2)
            
    # plot overlapping #    
    #f, (f1, f2) = plt.subplots(2, sharex = True, sharey = True)
    #plt.subplot(211)
    plt.bar(hist1.lbins, hist1.vd, width = w_overlap, color='k', alpha=1,    label= label1)
    plt.bar(hist2.lbins, hist2.vd, width = w_overlap, color='r', alpha=0.75, label= label2)
    #plt.bar(hist1.lbins, hist1.count_err, width = w_overlap, color='black', alpha=0.75, label= label2)
    #plt.bar(hist2.lbins, hist2.count_err, width = w_overlap, color='b', alpha=0.75, label= label2)
    plt.title('Line of Sight Vel Disp Distribution')
    plt.xlim((xlower, xupper))
    plt.ylim((0.0, ylimit))
    plt.ylabel('vel disp')
    plt.legend()
    plt.savefig(save_folder_ove + name + '_overlapping.png', format='png')
    plt.clf()
    #plt.show()
        
    # plot_adjacent #
    plt.subplot(211)
    #f, (f1, f2) = plt.subplots(2, sharex = True, sharey = True)
    plt.bar(hist1.lbins, hist1.vd, width = w_adjacent, color='b')
    plt.legend(handles=[mpatches.Patch(color='b', label= plot_hist1)])
    plt.title('Line of Sight Vel Disp Distribution')
    plt.xlim((xlower, xupper))
    plt.ylim((0.0, ylimit))
    plt.ylabel('counts')
    plt.xlabel('Lambda')

    plt.subplot(212)
    plt.bar(hist2.lbins, hist2.vd, width = w_adjacent, color='k')
    plt.legend(handles=[mpatches.Patch(color='k', label= plot_hist2)])
    plt.xlim((xlower, xupper))
    plt.ylim((0.0, ylimit))
    plt.xlabel('l')
    plt.ylabel('vel disp')
    plt.xlabel('Lambda')
    #f.subplots_adjust(hspace=0)
    plt.savefig(save_folder_adj + name + '.png', format='png')
    plt.clf()
    #plt.show()
    return 1
# # 

# # # # # # # # # # # # # # # # # # # # # #
#        NON-histogram plot               #
# # # # # # # # # # # # # # # # # # # # # #
# #
def vlos_plot_single(file1):
    ylimit = 100
    xlower = 180 
    xupper = -180
    w_overlap = 2.5
    w_adjacent = 1.5
    folder_hist = 'quick_plots/hists/'
    folder_outs = 'quick_plots/outputs/'
    save_folder_adj = 'quick_plots/comp_hist_plots/adj/'
    save_folder_ove = 'quick_plots/comp_hist_plots/overlap/'
    
    print "plot histogram 1: ", file1
    
    plot_hist1 = file1 + ".hist"
    
    output1 = file1 + ".out"
    
    label1 = '1'
    
    name = 'vlos_plots'
    print("plotting histograms\n")
    
    hist1 = nbody_histograms(folder_hist + plot_hist1)
     
    angle_cuttoffs = [-150.0, 150.0, 50, -15.0, 15.0, 1]
    
    out1 = nbody_outputs(folder_outs + output1)
    
    out1.binner_vlos(angle_cuttoffs)#bin the line of sight vels
    # plot overlapping #
    count_y_limit = 0.4
    rawcount_y_limit = 2000
    vel_disp_ylimit = 50
    
    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
    f.subplots_adjust(hspace=0)
    f.subplots_adjust(wspace=0)

    ax1 = plt.subplot(311)
    plt.bar(hist1.lbins, hist1.counts, width = w_adjacent, color='k', alpha=1)
    #plt.title(r'Line of Sight $\sigma_{line of sight}$ Distribution')
    plt.xlim((xlower, xupper))
    plt.ylim((0.0, count_y_limit))
    plt.ylabel('counts')
    plt.legend()
    
    #ax3 = plt.subplot(412)
    #plt.bar(hist1.lbins, hist1.count_err, width = w_adjacent, color='k', alpha=1)
    #plt.xlim((xlower, xupper))
    #plt.ylim((0.0, count_y_limit**0.5))
    #plt.ylabel('Count error')
    #plt.xlabel('Lambda')
    #plt.legend()
    
    ax2 = plt.subplot(312)
    plt.bar(hist1.lbins, hist1.vd, width = w_adjacent, color='k', alpha=1)
    plt.xlim((xlower, xupper))
    plt.ylim((0.0, vel_disp_ylimit))
    plt.ylabel(r'$\sigma$ (km/s)')
    #plt.legend()
    
    ax2 = plt.subplot(313)
    plt.scatter(out1.which_bin, out1.binned_vlos, s=2, marker= '.',  color='k', alpha=1, edgecolors='none')
    plt.xlim((xlower, xupper))
    #plt.ylim((0.0, vel_disp_ylimit))
    plt.ylabel(r'${v_{los}}$ (km/s)')
    #plt.legend()
    plt.xlabel(r'$\Lambda$')
    plt.savefig(save_folder_ove + name + '_overlapping_single.png', format='png', dpi=500)
    #plt.clf()
    #plt.show()
    
    return 1

def vlos_plot(file1, file2):
    ylimit = 100
    xlower = 180 
    xupper = -180
    w_overlap = 2.5
    w_adjacent = 1.5
    folder_hist = 'quick_plots/hists/'
    folder_outs = 'quick_plots/outputs/'
    save_folder_adj = 'quick_plots/comp_hist_plots/adj/'
    save_folder_ove = 'quick_plots/comp_hist_plots/overlap/'
    
    print "plot histogram 1: ", file1
    print "plot histogram 2: ", file2
    
    plot_hist1 = file1 + ".hist"
    plot_hist2 = file2 + ".hist"
    
    output1 = file1 + ".out"
    output2 = file2 + ".out"
    
    label1 = '1'
    label2 = '2'
    
    name = 'vlos_plots'
    print("plotting histograms\n")
    
    hist1 = nbody_histograms(folder_hist + plot_hist1)
    hist2 = nbody_histograms(folder_hist + plot_hist2)
     
    angle_cuttoffs = [-150.0, 150.0, 50, -15.0, 15.0, 1]
    
    out1 = nbody_outputs(folder_outs + output1)
    out2 = nbody_outputs(folder_outs + output2)
    
    out1.binner_vlos(angle_cuttoffs)#bin the line of sight vels
    out2.binner_vlos(angle_cuttoffs)
    
    # plot_adjacent #
    count_y_limit = 0.4
    rawcount_y_limit = 2000
    vel_disp_ylimit = 100
    
    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
    f.subplots_adjust(hspace=0)
    f.subplots_adjust(wspace=0)
    #plt.subplots(4, sharex = True, sharey = True)
    ax1 = plt.subplot(421)
    plt.bar(hist1.lbins, hist1.counts, width = w_adjacent, color='b')
    plt.title('Line of Sight Vel Disp Distribution')
    plt.xlim((xlower, xupper))
    plt.ylim((0.0, count_y_limit))
    plt.ylabel('counts')

    ax2 = plt.subplot(422)
    plt.bar(hist2.lbins, hist2.counts, width = w_adjacent, color='k')
    plt.xlim((xlower, xupper))
    plt.ylim((0.0, count_y_limit))
    plt.yticks([])

    ax5 = plt.subplot(423)
    plt.bar(hist1.lbins, hist1.count_err, width = w_adjacent, color='b')
    plt.xlim((xlower, xupper))
    plt.ylim((0.0, rawcount_y_limit))
    plt.ylabel('raw count')
    #plt.xlabel('Lambda')

    ax6 = plt.subplot(424)
    plt.bar(hist2.lbins, hist2.count_err, width = w_adjacent, color='k')
    plt.xlim((xlower, xupper))
    plt.ylim((0.0, rawcount_y_limit))
    plt.yticks([])
    #plt.xlabel('Lambda')
    
    ax3 = plt.subplot(425)
    #plt.subplots(2, sharex = True, sharey = False)
    plt.bar(hist1.lbins, hist1.vd, width = w_adjacent, color='b')
    plt.xlim((xlower, xupper))
    plt.ylim((0.0, vel_disp_ylimit))
    plt.ylabel('vel disp')

    ax4 = plt.subplot(426)
    plt.bar(hist2.lbins, hist2.vd, width = w_adjacent, color='k')
    plt.xlim((xlower, xupper))
    plt.ylim((0.0, vel_disp_ylimit))
    plt.yticks([])
    
    ax3 = plt.subplot(427)
    #plt.subplots(2, sharex = True, sharey = False)
    plt.scatter(out1.which_bin, out1.binned_vlos, color='b', s=.5, marker= 'o')
    plt.xlim((xlower, xupper))
    #plt.ylim((0.0, vel_disp_ylimit))
    plt.ylabel('vel disp')

    ax4 = plt.subplot(428)
    plt.scatter(out2.which_bin, out2.binned_vlos, color='k', s=.5, marker= 'o', )
    plt.xlim((xlower, xupper))
    #plt.ylim((0.0, vel_disp_ylimit))
    plt.yticks([])
    plt.savefig(save_folder_adj + name + '.png', format='png', dpi=1000)
    
     # plot overlapping #   
    count_y_limit = 0.4
    rawcount_y_limit = 2000
    vel_disp_ylimit = 100
    
    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
    f.subplots_adjust(hspace=0)
    f.subplots_adjust(wspace=0)

    ax1 = plt.subplot(411)
    plt.bar(hist1.lbins, hist1.counts, width = w_adjacent, color='k', alpha=1,    label= label1)
    plt.bar(hist2.lbins, hist2.counts, width = w_adjacent, color='r', alpha=0.75, label= label2)
    plt.title('Line of Sight Vel Disp Distribution')
    plt.xlim((xlower, xupper))
    plt.ylim((0.0, count_y_limit))
    plt.ylabel('counts')
    plt.legend()
    
    ax3 = plt.subplot(412)
    plt.bar(hist1.lbins, hist1.count_err, width = w_adjacent, color='k', alpha=1,    label= label1)
    plt.bar(hist2.lbins, hist2.count_err, width = w_adjacent, color='r', alpha=0.75, label= label2)
    plt.xlim((xlower, xupper))
    plt.ylim((0.0, rawcount_y_limit))
    plt.ylabel('raw count')
    plt.xlabel('Lambda')
    plt.legend()
    
    ax2 = plt.subplot(413)
    plt.bar(hist1.lbins, hist1.vd, width = w_adjacent, color='k', alpha=1,    label= label1)
    plt.bar(hist2.lbins, hist2.vd, width = w_adjacent, color='r', alpha=0.75, label= label2)
    plt.xlim((xlower, xupper))
    plt.ylim((0.0, vel_disp_ylimit))
    plt.ylabel('vel disp')
    plt.legend()
    
    ax2 = plt.subplot(414)
    plt.scatter(out1.which_bin, out1.binned_vlos, s=1, marker= '.',  color='red', alpha=1,label= label1, edgecolors='none')
    plt.scatter(out2.which_bin, out2.binned_vlos, s=1, marker= '.', color='blue', alpha=1, label= label2, edgecolors='none')
    plt.xlim((xlower, xupper))
    #plt.ylim((0.0, vel_disp_ylimit))
    plt.ylabel('vel disp')
    plt.legend()
    plt.savefig(save_folder_ove + name + '_overlapping.png', format='png', dpi=1000)
    #plt.clf()
    #plt.show()
    
    return 1

def lambda_beta_plot(file_name):
    path_charles = 'quick_plots/outputs/'
    path = 'quick_plots/'
    print file_name
    plot_lbr = y
    plot_light_and_dark = y
    plot_dm_alone = n
    
    out = nbody_outputs(path_charles + file_name + '.out')
    out.dark_light_split()
    out.convert_lambda_beta(True)
    
    fig = plt.figure()
    fig.subplots_adjust(hspace = 0.8, wspace = 0.8)
    # # # # # # # # # #
    if(plot_lbr):
        plt.figure(figsize=(10, 10))
        xlower = -180.0
        xupper = 180.0
        ylower = -80
        yupper = 80
        plt.xlim((xlower, xupper))
        plt.ylim((ylower, yupper))
        plt.xlabel(r'$\Lambda$')
        plt.ylabel(r'$\beta$')
        plt.title(r'$\Lambda$ vs $\beta$')
        #default to just plot lm
        plt.plot(out.light_lambda, out.light_beta, '.', markersize = .75, color = 'b', alpha=1.0, marker = '.',label = 'baryons')
        plt.legend()
        plt.savefig('/home/sidd/Desktop/research/quick_plots/' + file_name + '_lambdabeta_nodark', format='png')
        print "plotting:", len(out.light_l), " points"
        # # # # # # # # # #
        if(plot_light_and_dark):#plot lm and dm overlapping
            plt.xlim((xlower, xupper))
            plt.ylim((ylower, yupper))
            plt.xlabel(r'$\Lambda$')
            plt.ylabel(r'$\beta$')
            plt.title(r'$\Lambda$ vs $\beta$')
            plt.plot(out.dark_lambda, out.dark_beta, '.', markersize = .5, color = 'black', alpha=.75, marker = '.', label = 'dark matter')
            plt.legend()
            plt.savefig('/home/sidd/Desktop/research/quick_plots/' + file_name + '_lambdabeta', format='png')
            print "plotting:", len(out.light_l) + len(out.dark_l), " points"
        # # # # # # # # # #
        if(plot_dm_alone):#to plot just dm
            plt.clf()
            plt.figure(figsize=(20, 20))
            plt.xlim((xlower, xupper))
            plt.ylim((ylower, yupper))
            plt.xlabel(r'$\Lambda$')
            plt.ylabel(r'$\beta$')
            plt.title(r'$\Lambda$ vs $\beta$')
            plt.plot(out.dark_lambda, out.dark_beta, '.', markersize = 1, color = 'b', marker = '+')
            plt.legend()
            plt.savefig('/home/sidd/Desktop/research/quick_plots/tidal_stream_lambdabeta_dark', format='png')
            
    return 0 
 
def lb_plot(file_name):
    path_charles = 'quick_plots/outputs/'
    path = 'quick_plots/'
    print file_name
    plot_lbr = y
    plot_light_and_dark = y
    plot_dm = n
    plot_xyz = n
    plot_orbit = n
    
    out = nbody_outputs(path_charles + file_name + '.out')
    out.rescale_l()
    out.dark_light_split()
    
    fig = plt.figure()
    fig.subplots_adjust(hspace = 0.8, wspace = 0.8)
    # # # # # # # # # #
    if(plot_lbr):
        plt.figure(figsize=(10, 10))
        xlower = -180.0
        xupper = 180.0
        ylower = -80
        yupper = 80
        plt.xlim((xlower, xupper))
        plt.ylim((ylower, yupper))
        plt.xlabel('l')
        plt.ylabel('b')
        plt.title('l vs b')
        #default to just plot lm
        plt.plot(out.light_l, out.light_b, '.', markersize = 1., color = 'b', alpha=1.0, marker = '.', label = 'baryons')
        plt.legend()
        plt.savefig('/home/sidd/Desktop/research/quick_plots/' + file_name, format='png')
        print "plotting:", len(out.light_l), " points"
        # # # # # # # # # #
        if(plot_light_and_dark):#plot lm and dm overlapping
            plt.xlim((xlower, xupper))
            plt.ylim((ylower, yupper))
            plt.xlabel('l')
            plt.ylabel('b')
            plt.title('l vs b')
            plt.plot(out.dark_l, out.dark_b, '.', markersize = 1, color = 'red', alpha=.25, marker = '.', label = 'dark matter')
            plt.legend()
            plt.savefig('/home/sidd/Desktop/research/quick_plots/' + file_name, format='png')
            print "plotting:", len(out.light_l) + len(out.dark_l), " points"
        # # # # # # # # # #
        if(plot_orbit):
            orb = nbody_outputs(path + 'reverse_orbit.out')
            orb.rescale_l()
            orb.dark_light_split()
            
            plt.xlim((xlower, xupper))
            plt.ylim((ylower, yupper))
            plt.xlabel('l')
            plt.ylabel('b')
            plt.title('l vs b')
            plt.plot(orb.ls, orb.bs, '.', markersize = .15, color = 'g', alpha=1.0, marker = '.')
            #plt.savefig('/home/sidd/Desktop/research/quick_plots/tidal_stream_lbr_allmatter_orbit', format='png')
            
            
            orb = nbody_outputs(path + 'forward_orbit.out')
            orb.rescale_l()
            orb.dark_light_split()
                
            plt.xlim((xlower, xupper))
            plt.ylim((ylower, yupper))
            plt.xlabel('l')
            plt.ylabel('b')
            plt.title('l vs b')
            plt.plot(orb.ls, orb.bs, '.', markersize = .15, color = 'r', alpha=1.0, marker = '.')
            plt.savefig('/home/sidd/Desktop/research/quick_plots/tidal_stream_lbr_allmatter_orbit', format='png')
            plt.show()

        # # # # # # # # # #
        if(plot_dm):#to plot just dm
            plt.clf()
            plt.figure(figsize=(20, 20))
            plt.xlim((xlower, xupper))
            plt.ylim((ylower, yupper))
            plt.xlabel('l')
            plt.ylabel('b')
            plt.title('l vs b')
            plt.plot(out.dark_l, out.dark_b, '.', markersize = 1, color = 'b', marker = '+')
            plt.savefig('/home/sidd/Desktop/research/quick_plots/tidal_stream_lbr_dark', format='png')
            
    if(plot_xyz):
        xlower = 50
        xupper = -50
        fig.tight_layout()
        plt.axes().set_aspect('equal')
        plt.subplot(131, aspect='equal')
        plt.plot(out.light_x, out.light_y, '.', markersize = 1, color = 'r', marker = 'o')
        
        if(plot_light_and_dark == True):
            plt.plot(out.dark_x, out.dark_y, '.', markersize = 1, color = 'b', marker = '+')
        
        plt.xlim((xlower, xupper))
        plt.ylim((-80, 80))
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('x vs y')
        
        plt.subplot(132,aspect='equal')
        
        plt.plot(out.light_x, out.light_z, '.', markersize = 1, color = 'r', marker = 'o')
        
        if(plot_light_and_dark == True):
            plt.plot(out.dark_x, out.dark_z, '.', markersize = 1, color = 'b', marker = '+')
        
        plt.xlim((xlower, xupper))
        plt.ylim((-80, 80))
        plt.xlabel('x')
        plt.ylabel('z')
        plt.title('x vs z')
        
        plt.subplot(133, aspect='equal')
        
        plt.plot(out.light_z, out.light_y, '.', markersize = 1, color = 'r', marker = 'o')
        if(plot_light_and_dark == True):
            plt.plot(out.dark_z, out.dark_y, '.', markersize = 1, color = 'b', marker = '+')
        
        plt.xlim((xlower, xupper))
        plt.ylim((-80, 80))
        plt.xlabel('z')
        plt.ylabel('y')
        plt.title('z vs y')
        plt.savefig('/home/sidd/Desktop/research/quick_plots/tidal_stream_xyz', format='png')
    
    return 0

