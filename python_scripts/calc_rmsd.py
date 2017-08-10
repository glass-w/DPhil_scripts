import MDAnalysis as mda
import MDAnalysis.analysis.rms as rms
import numpy as np
import matplotlib.pyplot as plt

def rmsd_calc():
    ''' A function that calculates the RMSD of a user defined region
    and runs through a supplied trajectory
    '''

    #grofile = raw_input('Enter .gro filename with extension: ')
    #xtcfile = raw_input('Enter .xtc filname with extension: ')

    grofile = "repeat0_analysis_resid_90_190.gro"
    xtcfile = "repeat0_analysis_resid_90_190.xtc"

    u = mda.Universe(grofile, xtcfile)

    print ("Calculating RMSD...")

    selection = "resid 90:110 or resid 117:136 or resid 150:167 or resid 175:190 and name CA" # nav1.5 selection
    #selection = "resid 90:190 and name CA" # b3 selection	

    rmsd = []

    u.trajectory[1]
    ref = u.select_atoms(selection).positions # reference frame

    for i in range(len(u.trajectory)):

        u.trajectory[i]

        bb = u.select_atoms(selection).positions # updated frames

        rmsd.append((u.trajectory.time, rms.rmsd(ref, bb)))

        i += 1

        if i == len(u.trajectory)/4:
            print ("...25 %")
        elif i == len(u.trajectory)/2:
            print ("...50 %")
        elif i == len(u.trajectory)*(3/4):
            print ("...75 %")
        elif i == len(u.trajectory):
            print ("...100 %")

    rmsd = np.array(rmsd)

    np.savetxt("rmsd_data.csv", rmsd, delimiter=" ")
    
    return rmsd

def plot(data):

    ax = plt.subplot(111)
    ax.set_title("RMSD of the Mutant R3(0) Nav1.5 D1 S1-4")
    ax.plot((data[:,0]/1000), (data[:,1]/10), 'blue', lw=2, label=r"$R_G$", alpha = 0.3)
    ax.set_xlabel("Time / ns")
    ax.set_ylabel("RMSD from t = 0 / $\AA$")
    ax.figure.savefig("RMSD.png")
    plt.draw()

if __name__ == "__main__":

    data = rmsd_calc()
    fig = plot(data=data)





