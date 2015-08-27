import os
import os, sys
sys.path.append("/home/deovrat/PLUTO/Tools/pyPLUTO/pyPLUTO")
import numpy as np
import math
import pylab
import matplotlib.pyplot as plt
import pyPLUTO as pp

mu =0.62
mp=1.67e-24
kb=1.38e-16
kpc=3.086e21
for i in range(0,1):
	D = pp.pload(i) 
	den = D.rho 
	pres = D.prs
        T= pres*mu*mp*1.16e-7/(den*kb)
        x = D.x1
        y = D.x2
#        bx=D.bx1
#        by=D.bx2
	plt.clf() 
	plt.subplot(111, aspect = 'equal') 
        plt.xlabel('x (kpc)');
        plt.ylabel('y (kpc)');
	plt.contourf(x/kpc,y/kpc,np.log10(T)) 
	plt.axis([0, 40, 0, 40]) 
	plt.colorbar(orientation='vertical') 
        plt.hold(True);
#	plt.clim([-27, -24])
        nx=np.size(D.x1); ny=np.size(D.x2);
        sk=5;
        xp=D.x1[::sk]; yp=D.x2[::sk]; 
        bxp=np.zeros( ( (nx-1)/sk+1, (nx-1)/sk+1 ) );
        byp=np.zeros( ( (nx-1)/sk+1, (nx-1)/sk+1 ) );
        for k in range( (nx-1)/sk + 1 ):
                bxp[k][:]=D.bx1[sk*k][::sk];
                byp[k][:]=D.bx2[sk*k][::sk];
        plt.quiver(xp/kpc, yp/kpc, bxp/np.sqrt(bxp*bxp + byp*byp),
        byp/np.sqrt(bxp*bxp+byp*byp), scale=40.0)
	plt.show() 
	fname = 'time_step_%03d.png'%i
        print 'Saving frame', fname
        plt.savefig(fname)
        #plt.waitforbuttonpress()  
#        plt.clf()
