#plotEQ.py
#Description:   Plots Equilibrium (2D) from gfile for pyqt5 application
#Engineer:      T Looby
#Date:          20190916
import matplotlib.pyplot as plt
import numpy as np
import MDSplus
import EFIT.equilParams_class as EP
from scipy import interpolate
from scipy.interpolate import interp1d

def EQ2Dplot(ep,shot,t,MachFlag,height=None):
    #Seperatrix
    rbdry = ep.g['lcfs'][:,0]
    zbdry = ep.g['lcfs'][:,1]
    if MachFlag == 'nstx':
        rlim, zlim = nstxu_wall(oldwall=False) #FOR NSTXU
    else:
        rlim = ep.g['wall'][:,0]
        zlim = ep.g['wall'][:,1]
    #Change the linspace to get higher psi_n (outside wall)
    psi = ep.g['psiRZn']
    #psi = ep.g['psiRZ']
    levels = sorted(np.append([0.0,0.05,0.1,0.25,0.5,0.75,1.0], np.linspace(1.01,psi.max(),15)))
    R, Z = np.meshgrid(ep.g['R'],ep.g['Z'])
    #psiMax = ep.psiFunc.ev(max(rlim),0.0)
    psiMax = psi.max()
    commonLev = np.linspace(1.0,psiMax,15)
    lcfs = [1.0]
    if height is None:
        plt.figure(figsize=(5,8))
    else:
        dpi = 80
        #w = width / dpi
        #h = 8.0/5.0 * w
        h = height / dpi
        w = 5.0/8.0 * h
        plt.figure(figsize=(w,h), dpi=dpi)
    #Color Contour Plot
    #CS = plt.contourf(R,Z,psi,levels,cmap=plt.cm.bone)
    CS = plt.contourf(R,Z,psi,levels,cmap=plt.cm.cividis)
    #Draw Flux Lines in Common Flux Region
    plt.contour(CS, levels = commonLev, colors=('white',),linestyles='dotted',linewidths=(1,))
    #Draw separatrix as red line
    plt.contour(CS, levels = lcfs, colors=('r',),linestyles=('-',),linewidths=(2,))
    plt.axes().set_aspect('equal')
    #ax.set_aspect('equal')
    plt.xlabel('R [m]', fontsize=14,color='w')
    plt.ylabel('Z [m]', fontsize=14,color='w')
    plt.tick_params(axis='both',colors='w')
    #plt.xlim(min(rlim)-0.02,1.6)
#    plt.xlim(0.0,1.6)
#    plt.ylim(-1.7,1.7)
    plt.title("{:06d} @ {:05d}ms".format(shot,t), fontsize=14, color='white')
    #plt.colorbar(CS, label=r'$\psi$')
    #Fill in missing limiter section
    rlim_patched = np.append(rlim[2:], rlim[2])
    zlim_patched = np.append(zlim[2:], zlim[2])
    #plt.plot(rlim_patched, zlim_patched,'k--')
    plt.plot(rlim, zlim, '--', color='lime', lw=2)

    return plt


def nstxu_wall(oldwall=False):
    """
    returns simplified wall.  Uses two different wall versions
    """
    if oldwall:
        R = np.array([0.1851, 0.1851, 0.2794, 0.2794, 0.2979, 0.5712,
                    1.0433, 1.3192, 1.3358,
                    1.4851, 1.4791, 1.5174, 1.5313, 1.5464, 1.5608,
                    1.567, 1.5657, 1.5543, 1.5341, 1.5181, 1.4818,
                    1.4851, 1.3358, 1.3192, 1.0433,
                    0.5712, 0.2979, 0.2794, 0.2794, 0.1851, 0.1851])
        Z = np.array([0.0, 1.0081, 1.1714, 1.578, 1.6034, 1.6034,
                    1.43, 1.0397, 0.9976,
                    0.545, 0.4995, 0.306, 0.2355, 0.1586, 0.0801,
                    0.0, -0.0177, -0.1123, -0.221, -0.3026, -0.486,
                    -0.545, -0.9976, -1.0397, -1.43,
                    -1.6034, -1.6034, -1.578, -1.1714, -1.0081, 0])
    else:
      R = np.array([ 0.3147568,  0.3147568,  0.4441952,  0.4441952,  0.443484 ,
           0.443484 ,  0.6000496,  0.7672832,  0.8499856,  1.203452,  1.3192,  1.3358,  1.4851,  1.489 ,
           1.5638,  1.57  ,  1.5737,  1.575 ,  1.5737,  1.57  ,  1.5638,
           1.489 ,  1.4851,  1.3358,  1.3192,  1.203452 ,  0.8499856,  0.7672832,  0.6000496,  0.443484 ,
           0.443484 ,  0.4441952,  0.4441952,  0.3147568,  0.3147568 ])
      Z = np.array([ 0.       ,  1.0499344,  1.2899136,  1.5104872,  1.5104872,
            1.6028416,  1.6028416,  1.5367   ,  1.5367   ,  1.397508,  1.0397,  0.9976,  0.545 ,  0.49  ,
            0.1141,  0.0764,  0.0383,  0.    , -0.0383, -0.0764, -0.1141,
            -0.49  , -0.545 , -0.9976, -1.0397, -1.397508 , -1.5367   , -1.5367   , -1.6028416, -1.6028416,
            -1.5104872, -1.5104872, -1.2899136, -1.0499344,  0.])
    return R,Z
