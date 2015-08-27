import numpy as np
import math
import matplotlib.pyplot as plt
import os, sys
sys.path.append("/home/prateek/Desktop/Public_Codes/PLUTO41/Tools/pyPLUTO/pyPLUTO")
import pyPLUTO as pp 

imax = 200; kmax = 200; kpc=3.086e21;
for i in range(0,100): # 100 frames
        print i
        D=pp.pload(i);
        plt.clf()
        plt.subplot(111,aspect='equal');
        #title(r'$\log_{10}|B|$, $\Omega t/2\pi$=' + str(0.05*i), fontsize=16); 
        plt.xlabel('x (kpc)');
        plt.ylabel('y (kpc)');
        plt.contourf(D.x1/kpc,D.x2/kpc,np.log10(D.rho));
	#plt.clim([-25.6, -22.8]);
        plt.colorbar();
        plt.hold(True);
        nx=np.size(D.x1); ny=np.size(D.x2);
        sk=5;
        xp=D.x1[::sk]; yp=D.x2[::sk];
	#bx=D.bx1.reshape(nx,ny); by=D.bx2.reshape(nx,ny);
        #bxp=D.bx1[::4][::4]; byp=D.bx2[::4][::4];
	bxp=np.zeros( ( (nx-1)/sk+1, (nx-1)/sk+1 ) );
        byp=np.zeros( ( (nx-1)/sk+1, (nx-1)/sk+1 ) );
	for k in range( (nx-1)/sk + 1 ):
		bxp[k][:]=D.bx1[sk*k][::sk];
		byp[k][:]=D.bx2[sk*k][::sk]; 
        plt.quiver( xp/kpc, yp/kpc, np.divide(bxp, np.sqrt(bxp*bxp+byp*byp)), 
	np.divide(byp, np.sqrt(bxp*bxp+byp*byp)), scale=40. );
        #plt.contourf(x,z,np.log10(e));
        #plt.colorbar(); #plt.clim(0,2.e-18);
        #plt.axis([0,800,-800,800]);
        plt.show();
        fname = 'log10den%03d.png'%i
        print 'Saving frame', fname
        plt.savefig(fname)
        #plt.waitforbuttonpress(); 
        #plt.clf()
