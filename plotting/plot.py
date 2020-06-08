import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fname = "../examples/example_ic.json.nonlinrun"
fname_lin = "../examples/example_ic.json.run"
if __name__ == "__main__":

    data = np.loadtxt(fname, skiprows=1)
    time = data[:, 0]
    control = data[:, 13:]
    print("control.shape = ", control.shape)
    
    data_lin = np.loadtxt(fname_lin, skiprows=1)
    lin_ss = data_lin[0, :].reshape((1, 17))

    print("lin_ss = ", lin_ss)

    data_lin = data_lin[1:, :]
    print("data_lin[1:5]", data_lin[:5, :])
    time_lin = data_lin[:, 0]
    # exit(1)
    
    nlin = data_lin.shape[0]
    lin_ss = np.tile(lin_ss, ((nlin, 1)))
    print("lin_ss shape = ", lin_ss.shape)
    
    data_lin = data_lin + lin_ss  ## perturbation + steady_state
    control_lin = data_lin[:, 13:]

    print("\n after adding steady state")
    print("data_lin[1:5]", data_lin[-10:, :])
    # exit(1)

    ls_dyn = '-'
    c_dyn = 'r'
    c_lin = 'b'
    ls_lin = '--'

    labels = [["x", "y", "z"], ["U", "V", "W"], ["P", "Q", "R"], ["Roll", "Pitch", "Yaw"]] 
    

    fig, axs = plt.subplots(4, 4, figsize=(17,10))
    # States
    for ii in range(4):
        for jj in range(3):
            on_data = ii * 3 + jj + 1
            
            axs[ii, jj].plot(time, data[:, on_data], linestyle=ls_dyn, color=c_dyn, label='Nonlinear')
            axs[ii, jj].plot(time_lin, data_lin[:, on_data], linestyle=ls_lin, color=c_lin, label='Linear')    
            axs[ii, jj].set_ylabel(labels[ii][jj])
            if ii == 0 and jj == 0:
                axs[ii, jj].legend()

    # Controls
    labels = ["elevator", "aileron", "rudder", "thrust"]
    for jj in range(4):
        axs[jj, 3].plot(time, control[:, jj], linestyle=ls_dyn, color=c_dyn, label='Nonlinear')
        axs[jj, 3].plot(time_lin, control_lin[:, jj], linestyle=ls_lin, color=c_lin, label='Linear')    
        axs[jj,  3].set_ylabel(labels[jj])

    plt.show()
    # axs[0,1].plot(time, data[:, 2], linestyle=ls_dyn, color=c_dyn)
    # axs[0,1].plot(time_lin, data_lin[:, 2], linestyle=ls_lin, color=c_lin)
    # axs[0,1].set_ylabel('y')

    # axs[0,2].plot(time, data[:, 3], linestyle=ls_dyn, color=c_dyn)
    # axs[0,2].plot(time_lin, data_lin[:, 3], linestyle=ls_lin, color=c_lin)    
    # axs[0,2].set_ylabel('z')

    # axs[1,0].plot(time, data[:, 4], linestyle=ls_dyn, color=c_dyn)
    # axs[1,0].plot(time_lin, data_lin[:, 4], linestyle=ls_lin, color=c_lin)        
    # axs[1,0].set_ylabel('U')

    # axs[1,1].plot(time, data[:, 5], linestyle=ls_dyn, color=c_dyn)
    # axs[1,1].plot(time_lin, data_lin[:, 5], linestyle=ls_lin, color=c_lin)            
    # axs[1,1].set_ylabel('V')

    # axs[1,2].plot(time, data[:, 6], linestyle=ls_dyn, color=c_dyn)
    # axs[1,2].plot(time_lin, data_lin[:, 6], linestyle=ls_lin, color=c_lin)
    # axs[1,2].set_ylabel('W')

    # axs[2,0].plot(time, data[:, 7], linestyle=ls_dyn, color=c_dyn)
    # axs[2,0].plot(time_lin, data_lin[:, 7], linestyle=ls_lin, color=c_lin)
    # axs[2,0].set_ylabel('P')

    # axs[2,1].plot(time, data[:, 8], linestyle=ls_dyn, color=c_dyn)
    # axs[2,1].plot(time_lin, data_lin[:, 8], linestyle=ls_lin, color=c_lin)    
    # axs[2,1].set_ylabel('Q')

    
    # axs[2,2].plot(time, data[:, 9], linestyle=ls_dyn, color=c_dyn)
    # axs[2,2].plot(time_lin, data_lin[:, 9], linestyle=ls_lin, color=c_lin)
    # axs[2,2].set_ylabel('R')

    # axs[3,0].plot(time, data[:, 10], linestyle=ls_dyn, color=c_dyn)
    # axs[3,0].plot(time_lin, data_lin[:, 10], linestyle=ls_lin, color=c_lin)
    # axs[3,0].set_ylabel('Roll')

    # axs[3,1].plot(time, data[:, 11], linestyle=ls_dyn, color=c_dyn)
    # axs[3,1].plot(time_lin, data_lin[:, 11], linestyle=ls_lin, color=c_lin)    
    # axs[3,1].set_ylabel('Pitch')

    # axs[3,2].plot(time, data[:, 12], linestyle=ls_dyn, color=c_dyn)
    # axs[3,2].plot(time_lin, data_lin[:, 12], linestyle=ls_lin, color=c_lin)        
    # axs[3,2].set_ylabel('Yaw')               


    # fig = plt.figure()
    # ax = fig.add_subplot(111,projection='3d')
    # ax.plot(data[:, 1], data[:, 2], data[:,3])
    # ax.set_xlabel('x')
    # ax.set_ylabel('y')
    # ax.set_zlabel('z')
    
    # plt.show()
    
