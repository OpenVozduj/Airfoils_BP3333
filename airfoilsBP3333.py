import numpy as np
import BPbs as bp
import matplotlib.pyplot as plt

#%% Variables BP3333
Cj = np.array([-0.0011,   #r_le
               0.4486,    #x_t
               0.0401,    #y_t
               -0.2847,   #k_t
               0.0380,    #beta_te
               0.3145,    #gamma_le
               0.5141,    #x_c
               0.0140,    #y_c
               -0.0934,   #k_c
               0.4446,    #alpha_te
               0.001,     #dz_te
               0.0])      #z_te

nc = 24 #Points per construction curve


#%% Building airfoil
u = np.zeros(nc)
w = np.zeros(nc)
for i in range(nc):
    u[i] = 1-np.cos((i*np.pi)/(2*(nc-1)))
    w[i] = np.sin((i*np.pi)/(2*(nc-1)))
b9 = bp.BPb9(Cj[0], Cj[1], Cj[2], Cj[3])
b1 = bp.BPb1(Cj[5], Cj[7], Cj[8], Cj[11], Cj[9])
xlt, ylt, xtt, ytt = bp.Pgrosor(b9,Cj[1],Cj[2],Cj[3],Cj[4],Cj[10])
xlc, ylc, xtc, ytc = bp.Pconvadura(b1,Cj[6],Cj[7],Cj[8],Cj[9],Cj[5],Cj[11])
xul = bp.Bezier3(xlt,u)
yul = bp.Bezier3(ylt,u)
xut = bp.Bezier3(xtt,w[1:])
yut = bp.Bezier3(ytt,w[1:])
Xg = np.concatenate((xul,xut), axis=0)
Yg = np.concatenate((yul,yut), axis=0)
xcl = bp.Bezier3(xlc,u)
ycl = bp.Bezier3(ylc,u)
xct = bp.Bezier3(xtc,w[1:])
yct = bp.Bezier3(ytc,w[1:])
Xc = np.concatenate((xcl,xct), axis=0)
Yc = np.concatenate((ycl,yct), axis=0)
theta = bp.funtheta(Xc,Yc)
XU, YU, XL, YL = bp.coorper(Xg,Yg,Xc,Yc,theta) #Airfoil points

#%% Creation File .txt
nl = np.size(XU)
Upper = ['Upper \n']
for i in range(nl):
    Upper.append(str(XU[i])+'\t'+str(YU[i])+'\n')
Lower = ['Lower \n']
for i in range(nl):
    Lower.append(str(XL[i])+'\t'+str(YL[i])+'\n')

Coord = Upper + Lower
gmt = open('coordAirfoil.txt','w')
gmt.writelines(Coord)
gmt.close()

#%% Display airfoil
fig = plt.figure(3)
ax = fig.add_subplot(111)
ax.plot(XL, YL, '-r')
ax.plot(XU, YU, '-g')
ax.set_xlabel('x/b', fontsize=14)
ax.set_ylabel('y/b', fontsize=14)
ax.grid()
for label in (ax.get_xticklabels() + ax.get_yticklabels()):
	label.set_fontsize(12)
ax.set_aspect('equal')
