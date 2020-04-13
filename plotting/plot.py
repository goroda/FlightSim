import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fname = "../nrb_dyn.dat"

if __name__ == "__main__":

    data = np.loadtxt(fname, skiprows=1)

    time = data[:, 0]

    fig, axs = plt.subplots(4,3,figsize=(10,8))
    axs[0,0].plot(time, data[:, 1])
    axs[0,0].set_ylabel('x')
    axs[0,1].plot(time, data[:, 2])
    axs[0,1].set_ylabel('y')
    axs[0,2].plot(time, data[:, 3])
    axs[0,2].set_ylabel('z')

    axs[1,0].plot(time, data[:, 4])
    axs[1,0].set_ylabel('U')
    axs[1,1].plot(time, data[:, 5])
    axs[1,1].set_ylabel('V')
    axs[1,2].plot(time, data[:, 6])
    axs[1,2].set_ylabel('W')

    axs[2,0].plot(time, data[:, 7])
    axs[2,0].set_ylabel('P')
    axs[2,1].plot(time, data[:, 8])
    axs[2,1].set_ylabel('Q')
    axs[2,2].plot(time, data[:, 9])
    axs[2,2].set_ylabel('R')

    axs[3,0].plot(time, data[:, 10])
    axs[3,0].set_ylabel('Roll')
    axs[3,1].plot(time, data[:, 11])
    axs[3,1].set_ylabel('Pitch')
    axs[3,2].plot(time, data[:, 12])
    axs[3,2].set_ylabel('Yaw')               


    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')
    ax.plot(data[:, 1], data[:, 2], data[:,3])
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    
    plt.show()
    
