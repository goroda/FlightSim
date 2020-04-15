import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fname = "../nrb_dyn.dat"
fname_lin = "../lrb_dyn.dat"
if __name__ == "__main__":

    data = np.loadtxt(fname, skiprows=1)

    data_lin = np.loadtxt(fname_lin, skiprows=1)
    nlin = int(data_lin.shape[0] / 2)
    time_lin = data_lin[:nlin,0]
    data_lin = data_lin[:nlin, :] + data_lin[nlin:, :] ## perturbation + steady_state
    # data_lin = data_lin[:nlin, :] # perturbation
    # data_lin = data_lin[nlin:, :] # steady_state    
    
    
    time = data[:, 0]

    ls_dyn = '-'
    c_dyn = 'r'
    c_lin = 'b'
    ls_lin = '--'

    
    fig, axs = plt.subplots(4,3,figsize=(10,8))
    axs[0,0].plot(time, data[:, 1], linestyle=ls_dyn, color=c_dyn, label='Nonlinear')
    axs[0,0].plot(time_lin, data_lin[:, 1], linestyle=ls_lin, color=c_lin, label='Linear')    
    axs[0,0].set_ylabel('x')
    axs[0,0].legend()
    
    axs[0,1].plot(time, data[:, 2], linestyle=ls_dyn, color=c_dyn)
    axs[0,1].plot(time_lin, data_lin[:, 2], linestyle=ls_lin, color=c_lin)
    axs[0,1].set_ylabel('y')

    axs[0,2].plot(time, data[:, 3], linestyle=ls_dyn, color=c_dyn)
    axs[0,2].plot(time_lin, data_lin[:, 3], linestyle=ls_lin, color=c_lin)    
    axs[0,2].set_ylabel('z')

    axs[1,0].plot(time, data[:, 4], linestyle=ls_dyn, color=c_dyn)
    axs[1,0].plot(time_lin, data_lin[:, 4], linestyle=ls_lin, color=c_lin)        
    axs[1,0].set_ylabel('U')

    axs[1,1].plot(time, data[:, 5], linestyle=ls_dyn, color=c_dyn)
    axs[1,1].plot(time_lin, data_lin[:, 5], linestyle=ls_lin, color=c_lin)            
    axs[1,1].set_ylabel('V')

    axs[1,2].plot(time, data[:, 6], linestyle=ls_dyn, color=c_dyn)
    axs[1,2].plot(time_lin, data_lin[:, 6], linestyle=ls_lin, color=c_lin)
    axs[1,2].set_ylabel('W')

    axs[2,0].plot(time, data[:, 7], linestyle=ls_dyn, color=c_dyn)
    axs[2,0].plot(time_lin, data_lin[:, 7], linestyle=ls_lin, color=c_lin)
    axs[2,0].set_ylabel('P')

    axs[2,1].plot(time, data[:, 8], linestyle=ls_dyn, color=c_dyn)
    axs[2,1].plot(time_lin, data_lin[:, 8], linestyle=ls_lin, color=c_lin)    
    axs[2,1].set_ylabel('Q')

    
    axs[2,2].plot(time, data[:, 9], linestyle=ls_dyn, color=c_dyn)
    axs[2,2].plot(time_lin, data_lin[:, 9], linestyle=ls_lin, color=c_lin)
    axs[2,2].set_ylabel('R')

    axs[3,0].plot(time, data[:, 10], linestyle=ls_dyn, color=c_dyn)
    axs[3,0].plot(time_lin, data_lin[:, 10], linestyle=ls_lin, color=c_lin)
    axs[3,0].set_ylabel('Roll')

    axs[3,1].plot(time, data[:, 11], linestyle=ls_dyn, color=c_dyn)
    axs[3,1].plot(time_lin, data_lin[:, 11], linestyle=ls_lin, color=c_lin)    
    axs[3,1].set_ylabel('Pitch')

    axs[3,2].plot(time, data[:, 12], linestyle=ls_dyn, color=c_dyn)
    axs[3,2].plot(time_lin, data_lin[:, 12], linestyle=ls_lin, color=c_lin)        
    axs[3,2].set_ylabel('Yaw')               


    # fig = plt.figure()
    # ax = fig.add_subplot(111,projection='3d')
    # ax.plot(data[:, 1], data[:, 2], data[:,3])
    # ax.set_xlabel('x')
    # ax.set_ylabel('y')
    # ax.set_zlabel('z')
    
    plt.show()
    
