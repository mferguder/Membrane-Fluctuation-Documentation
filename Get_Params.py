
#----------------------------------------------------------------------
def hq_func(x,a,b,c): return (a+b*np.power(x,2))/(1-np.power(x,2)*c)
def nqll_func(x,a,b): return a/(1-b*np.power(x,2))
def nqL_func(x,a,b): return 1/(1/a+b*np.power(x,2))
def nllL_func(x,a,b,c,d):
    for i in range(1,len(x)):
        if x[i]==x[0]:
            x2=x[:i]
            x3=x[i:]   
    nll=nqll_func(x2,a,b/a*np.power(c,2)/4)
    nL=nqL_func(x3,b,d)
    return list(nll) + list(nL)
#----------------------------------------------------------------------
def Bootstrp(qx,zx,dzx,xq_func,px):
    badness=0.999#9
    Nboot=2000
    div_chi=1

    a,b = curve_fit(xq_func,   qx, zx, px, sigma=dzx)

    params_xq=np.zeros((Nboot,len(a)))
    dparams_xq=np.zeros((Nboot,len(a),len(a)))

    chi_s=np.power(np.divide(np.subtract(zx,xq_func(qx,*a)),dzx),2)
    c=sc.gammainc((len(qx)-len(a))/2.,sum(chi_s)/4./div_chi)   
    #     c=stats.chi2.cdf(sum(chi_s) , len(qx)-len(a)-1)

    while c>badness:
        qx=qx[:-1]
        zx=zx[:-1]
        dzx=dzx[:-1]
        a,b = curve_fit(xq_func,   qx, zx, px, sigma=dzx)
        chi_s=np.power(np.divide(np.subtract(zx,xq_func(qx,*a)),dzx),2)
        c=sc.gammainc((len(qx)-len(a))/2.,sum(chi_s)/4./div_chi)
    #         c=stats.chi2.cdf(sum(chi_s) , len(qx)-len(a)-1)
    print(len(qx),c,np.mean(abs(np.divide(np.subtract(zx,xq_func(qx,*a)),dzx))) )

    for i in range(Nboot):
        zxk=np.random.normal(zx,dzx,len(qx))
        try: params_xq[i], dparams_xq[i] = curve_fit(xq_func,   qx, zxk, px, sigma=dzx)
        except RuntimeError: print("Error - curve_fit failed")
    params_xq =np.transpose(params_xq)
    qzdzx=np.concatenate((np.array([qx]),np.array([zx]),np.array([dzx])))
#     print params_xq,qzdzx
    return params_xq,qzdzx
#----------------------------------------------------------------------
def Bootstrp_2D(qx,qy,zx,dzx,xq_func,_xq_func,px,qmax):
    badness=0.999#9
    Nboot=2000
    div_chi=1

    xdatar=np.vstack((qx, qy))
    a,b = curve_fit(_xq_func,   xdatar, zx, px, sigma=dzx)
    chi_s=np.power(np.divide(np.subtract(zx,xq_func(qx,qy,*a)),dzx),2)
    c=sc.gammainc((len(qx)-len(a))/2.,sum(chi_s)/4./div_chi)   

    qzdzx=np.concatenate((np.array([qx]),np.array([qy]),np.array([zx]),np.array([dzx])))
    params_xq=np.zeros((Nboot,len(a)))
    dparams_xq=np.zeros((Nboot,len(a),len(a)))

    while c>badness or qmax**2<qx[-1]**2+qy[-1]**2:
    #     print('c='+str(c))
    #     if c<1: print( len(qx),c,np.mean(abs(np.divide(np.subtract(zx,xq_func(qx,qy,*a)),dzx))))
        qx=qx[:-1]
        qy=qy[:-1]
        zx=zx[:-1]
        dzx=dzx[:-1]
        xdatar=np.vstack((qx, qy))
        a,b = curve_fit(_xq_func,   xdatar, zx, px, sigma=dzx)
        chi_s=np.power(np.divide(np.subtract(zx,xq_func(qx,qy,*a)),dzx),2)
        c=sc.gammainc((len(qx)-len(a))/2.,sum(chi_s)/4./div_chi)
    #         c=stats.chi2.cdf(sum(chi_s) , len(qx)-len(a)-1)
    print( len(qx),c,np.mean(abs(np.divide(np.subtract(zx,xq_func(qx,qy,*a)),dzx))))
    q=np.sqrt(qx**2+qy**2)
    zx=zx[q.argsort()]
    dzx=dzx[q.argsort()]
    qx=qx[q.argsort()]
    qy=qy[q.argsort()]
    q=q[q.argsort()]
    qzdzx=np.concatenate((np.array([qx]), np.array([qy]), np.array([q]), np.array([zx]), np.array([dzx])))

    i=0
    while i<Nboot:
        zxk=np.random.normal(zx,dzx,len(qx))
        try:
            params_xq[i], dparams_xq[i] = curve_fit(_xq_func,   xdatar, zxk, px, sigma=dzx)
            i+=1
        except RuntimeError: print("Error - curve_fit failed")
    params_xq =np.transpose(params_xq)
    return params_xq,qzdzx
#----------------------------------------------------------------------
def Bootstrp2(qx,qy,zx,dzx,zy,dzy,xq_func,px,kx,ky,badness=0.9999):
    Nboot=2000
    qxn=qx*1.
    qyn=qy*1.
    zxn=zx*1.
    zyn=zy*1.
    dzxn=dzx
    dzyn=dzy
    qqn= np.concatenate((qxn,qyn))
    zzn = np.concatenate((zxn,zyn))
    dzzn= np.concatenate((dzxn,dzyn))
#     print len(qqn),len(qxn),len(qyn), zzn

    a,b = curve_fit(xq_func,   qqn, zzn, px, sigma=dzzn, maxfev=1100)
    chi_s=np.power(np.divide(np.subtract(zzn,xq_func(qqn,*a)),dzzn),2)
    c=max(sc.gammainc((len(qxn)-kx)/2.,sum(chi_s[:len(qxn)])/4.),
          sc.gammainc((len(qyn)-ky)/2.,sum(chi_s[len(qxn):])/4.))
#     c=max(stats.chi2.cdf(sum(chi_s[:len(qxn)]) , len(qxn)-2-1),
#           stats.chi2.cdf(sum(chi_s[len(qxn):]) , len(qyn)-2-1))
    
    params_xq=np.zeros((Nboot,len(a)))
    dparams_xq=np.zeros((Nboot,len(a),len(a)))
    while c>badness:
        #print('nll&nL sim. fit is worse than seperate')
        which=cmp(sc.gammainc((len(qxn)-kx)/2.,sum(chi_s[:len(qxn)])/4.),
                  sc.gammainc((len(qyn)-ky)/2.,sum(chi_s[len(qxn):])/4.))
#         which=cmp(stats.chi2.cdf(sum(chi_s[:len(qxn)]) , len(qxn)-2-1),
#                   stats.chi2.cdf(sum(chi_s[len(qxn):]) , len(qyn)-2-1))
        if which>=0:
            qxn=qxn[:-1]
            zxn=zxn[:-1]
            dzxn=dzxn[:-1]
        if which<=0:
            qyn=qyn[:-1]
            zyn=zyn[:-1]
            dzyn=dzyn[:-1]
        qqn= np.concatenate((qxn,qyn)) 
        zzn = np.concatenate((zxn,zyn))
        dzzn= np.concatenate((dzxn,dzyn))
        a,b = curve_fit(xq_func,   qqn, zzn, px, sigma=dzzn)
        chi_s=np.power(np.divide(np.subtract(zzn,xq_func(qqn,*a)),dzzn),2)
        c=max(sc.gammainc((len(qxn)-kx)/2.,sum(chi_s[:len(qxn)])/4.),
              sc.gammainc((len(qyn)-ky)/2.,sum(chi_s[len(qxn):])/4.))
#         c=max(stats.chi2.cdf(sum(chi_s[:len(qxn)]) , len(qxn)-2-1),
#               stats.chi2.cdf(sum(chi_s[len(qxn):]) , len(qyn)-2-1))
    print(len(qxn),len(qyn),c,np.mean(abs(np.divide(np.subtract(zzn,xq_func(qqn,*a)),dzzn))))

    for i in range(Nboot):
        zznk=np.random.normal(zzn,dzzn,len(qqn))
        try: params_xq[i], dparams_xq[i] = curve_fit(xq_func,   qqn, zznk, px, sigma=dzzn)
        except RuntimeError: print("Error - curve_fit failed")
    params_xq =np.transpose(params_xq)
    qzdzx=np.concatenate((np.array([qqn]),np.array([zzn]),np.array([dzzn]))) 
    return params_xq,qzdzx
#----------------------------------------------------------------------
   
   
params_hq,qzdz1=Bootstrp(q,z1,dz1,hq_func,ph)
params_nqll,qzdz2=Bootstrp(q,z2,dz2,nqll_func,pll)
params_nqL,qzdz3=Bootstrp(q,z3,dz3,nqL_func,pL);
params_nllL,qzdz4=Bootstrp2(q,z2,dz2,z3,dz3,nllL_func,pllL,2,2);
#params_hll,qzdz5=Bootstrp2(q,z2,dz2,z1,dz1,hll_func,phll,3,2);
