def BPb9(rle, xt, yt, kt):
    import numpy as np
    def funb9(x):
        f = (27/4)*kt**2*x**4-27*kt**2*xt*x**3+(9*kt*yt+(81/2)*kt**2*xt**2)*x**2+(2*rle-18*kt*xt*yt-27*kt**2*xt**3)*x+(3*yt**2+9*kt*xt**2*yt+(27/4)*kt**2*xt**4)
        return f
    vcheck = xt-np.sqrt(-2*yt/(3*kt))
    limi = max([0, vcheck])
    lims = xt
    xi = (lims-limi)/2
    dx = 0.01
    es = 0.01
    for i in range(50):
        xi1 = xi-funb9(xi)/((funb9(xi+dx)-funb9(xi-dx))/(2*dx))
        ea = np.abs((xi1-xi)/xi1)
        if ea <= es:
            break
        else:
            xi = xi1;
    b9 = xi1
    return b9

def Bezier3(z, u):
    zu = z[0]*(1-u)**3 + 3*z[1]*u*(1-u)**2 + 3*z[2]*u**2*(1-u)+z[3]*u**3
    return zu

def BPb1(gammale, yc, kc, zte, alphate):
    import numpy as np
    def cotg(alpha):
        c = 1/np.tan(alpha)
        return c
    t1 = 3*kc*(cotg(gammale)+cotg(alphate))**2
    t2 = 16 + 3*kc*(cotg(gammale)+cotg(alphate))*(1+zte*cotg(alphate))
    t3 = 4*np.sqrt(16 + 6*kc*(cotg(gammale)+cotg(alphate))*(1-yc*(cotg(gammale)+cotg(alphate))+zte*cotg(alphate)))
    b = (t2 + t3)/t1
    if b<yc and b>0.0:
        b1 = b
    else:
        b1 = (t2 - t3)/t1
    return b1

def funtheta(X,Y):
    import numpy as np
    gr = np.gradient(Y,X)
    theta = np.arctan(gr)
    return theta 

def checkpossx(a):
    v = a.size
    s = 0
    for i in range(v-1):
        if a[i]<a[i+1]:
            s = s + 0
        else:
            s = s + 1
    return s

def Pgrosor(b9, xt, yt, kt, betate, dzte):
    import numpy as np
    xlt = np.array([0,
                    0,
                    b9,
                    xt])
    ylt = np.array([0,
                    1.5*kt*(xt-b9)**2+yt,
                    yt,
                    yt])
    xtt = np.array([xt,
                    2*xt-b9,
                    1+(dzte-(1.5*kt*(xt-b9)**2+yt))/np.tan(betate),
                    1.0])
    ytt = np.array([yt,
                    yt,
                    1.5*kt*(xt-b9)**2+yt,
                    dzte])
    return xlt, ylt, xtt, ytt

def Pconvadura(b1, xc, yc, kc, alphate, gammale, zte):
    import numpy as np
    xlc = np.array([0.0,
                    b1/np.tan(gammale),
                    xc-np.sqrt(2*(b1-yc)/(3*kc)),
                    xc])
    ylc = np.array([0.0,
                    b1,
                    yc,
                    yc])
    xtc = np.array([xc,
                    xc+np.sqrt(2*(b1-yc)/(3*kc)),
                    1+(zte-b1)/np.tan(alphate),
                    1.0])
    ytc = np.array([yc,
                    yc,
                    b1,
                    zte])
    return xlc, ylc, xtc, ytc

def coorper(Xg, Yg, Xc, Yc, theta):
    import numpy as np
    from scipy.interpolate import InterpolatedUnivariateSpline
    t = np.size(Xg)
    Yin = np.zeros(t)
    for i in range(t):
        if Xc[i]==Xg[i]:
            Yin[i] = Yg[i]
        else:
            P = InterpolatedUnivariateSpline(Xg,Yg,k=3)
            Yin[i] = P(Xc[i])
    XU = Xc-Yin*np.sin(theta)
    YU = Yc+Yin*np.cos(theta)
    XL = Xc+Yin*np.sin(theta)
    YL = Yc-Yin*np.cos(theta)
    return XU, YU, XL, YL
