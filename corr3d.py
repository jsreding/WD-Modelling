def corr3d(teff, logg):
    import numpy as np
    Tx = (teff - 10000.)/1000.
    gx = logg - 8
    a = [-1.0461690e-3, -2.6846737e-1, 3.0654611e-1, 1.8025848, 1.5006909e-1, 1.0125295e-1, -5.2933335e-2, -1.3414353e-1]
    b = [1.1922481e-3, -2.7230889e-1, -6.7437328e-2, -8.7753624e-1, 1.4936511e-1, -1.9749393e-1, 4.1687626e-1, 3.8195432e-1, -1.4141054e-1, -2.9439950e-2, 1.1908339e-1]
    delT = (a[0] + (a[1] + a[4]*Tx + (a[5] + a[6]*Tx + a[7]*gx)*gx)*np.exp(-a[2]*(Tx-a[3])**2))*1000
    dellg = b[0] + b[1]*np.exp(-(b[2] + (b[4] + b[6]*np.exp(-b[7]*(Tx - b[8])**2))*Tx + (b[5] + b[9]*Tx + b[10]*gx)*gx)**2*(Tx - b[3])**2)
    tfin = teff + delT
    lgfin = logg + dellg
    return tfin, lgfin
