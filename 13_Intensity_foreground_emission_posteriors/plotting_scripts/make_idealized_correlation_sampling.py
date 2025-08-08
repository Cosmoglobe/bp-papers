import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as scconst
from matplotlib.cm import get_cmap

plttype='.pdf'    
msize=7
N=10000
fs=16
mean_x=0.0
mean_y=0.0
std_x=1.0
std_y=1.0
corr=0.99
CLs=np.array([0.95,0.68])
fh=2 #height value of any conditional

samp_frac=0.1
xmin_plt=-4
xmax_plt=4
ymin_plt=-3
ymax_plt=4
p_init=np.array([1.5,1.3])

#plotting order for scatter
scat='cond first' #{'marg first','cond first'}

#pf max
f_min=fh/200.0


def asymp_line(x,corr,mx=0.0,my=0.0,sx=1.0,sy=1.0):
    if (abs(corr) > 0.0):
        y = my + np.sign(corr)*(sy/sx) * (x-mx)
    else:
        y = my
    return y


def make_circle(r,N=1000):
    theta=np.linspace(0,2*scconst.pi,N)
    theta[-1]=0.0 #make sure we loop all the way
    x=r*np.cos(theta)
    y=r*np.sin(theta)
    return x,y

def f_bivar(x,y,corr,mx=0.0,my=0.0,sx=0.0,sy=0.0):
    t1 = 1.0/(2.0*scconst.pi*sx*sy*np.sqrt(1-corr**2))
    t2 = -1/(2*(1.0-corr**2))
    t3 = ( (x-mx)/sx )**2 - 2*corr*((x-mx)/sx)*((y-my)/sy) + ( (y-my)/sy )**2

    return t1*np.exp(t2*t3)

def f_bivar_marg(par,mp=0.0,sp=1.0):
    t1 = 1.0/(np.sqrt(2.0*scconst.pi)*sp)
    t2 = -1/(2.0*sp**2)
    t3 = (par-mp)**2

    return t1*np.exp(t2*t3)

def get_conditional_likelihood(x,y,mu_x=0.0,mu_y=0.0,sig_x=1.0,sig_y=1.0,corr=0.0):
    # given inputs: x, y, axis, mu, sigma, correlation, 
    M=np.zeros([2,1])
    M[0,0]=mu_x*1.0
    M[1,0]=mu_y*1.0
    S=np.zeros([2,2],dtype='double')
    S[0,0]=sig_x**2
    S[1,1]=sig_y**2
    S[0,1]=sig_x*sig_y*corr
    S[1,0]=sig_x*sig_y*corr

    iS = np.linalg.inv(S)
    denom=2*scconst.pi*sig_x*sig_y*np.sqrt(1.0-corr**2)
    if (denom > 0.0):
        A=1.0/denom
    else:
        return 0.0

    X=np.zeros([2,1])
    X[0,0]=x
    X[1,0]=y

    Y = X-M
    
    ex = -1.0/2.0 * Y.T@iS@Y
    temp=A*np.exp(ex)
    return temp[0,0]


def get_CL_line_rot(cl,corr,mx=0.0,my=0.0,sx=1.0,sy=1.0,N=1000):
    if (abs(corr)>= 1.0):
        print('No level possible. The level is a straight line')
        exit()

    if (cl >= 1.0):
        print('No contour level possible for CL:',cl)
        return np.array([-1e30]),np.array([-1e30])
    elif (cl <= 0.0):
        print('No contour level possible for CL:',cl)
        return np.array([-1e30]),np.array([-1e30])
    
    r = np.sqrt(-2*np.log(1.0-cl))
    
    S=np.array([[sx**2, sx*sy*corr],
               [corr*sx*sy, sy**2]])

    #iS = np.linalg.inv(S)
    #print(S)
    #print(iS)
    #print(S@iS)
    
    w,v=np.linalg.eig(S)

    D = np.zeros([2,2])
    D[0,0]=w[0]
    D[1,1]=w[1]

    Dsq=np.sqrt(D)
    
    Ssqrt=v@Dsq@v.T
    vec=np.zeros([2,1])
    xarr,yarr=make_circle(1.0,N=N)
    mvec=np.array([[mx],[my]])
    for i in range(N):
        vec[0,0]=xarr[i]*1.0
        vec[1,0]=yarr[i]*1.0

        Y = r*Ssqrt@vec + mvec
        xarr[i]=Y[0,0].copy()
        yarr[i]=Y[1,0].copy()

    return xarr,yarr


t20=get_cmap('tab20c')
col20=[]
for i in range(20):
    col20.append(t20(i/20.0))


############## Get points ###################

#get new conditional sample
xf=np.linspace(xmin_plt,xmax_plt,N)
yf=xf*0.0 + p_init[1]
yf=f_bivar(xf,yf,corr,mx=mean_x,my=mean_y,sx=std_x,sy=std_y)

Pf=yf/np.sum(yf) #normalize
j=-1
for i in range(1,N):
    Pf[i]=Pf[i]+Pf[i-1]
    if (j<0):
        if Pf[i] > samp_frac:
            j=i
p_cond=np.array([0.0,0.0])
p_cond[0]=xf[j]

yfc=np.linspace(ymin_plt,ymax_plt,N)
xfc=yfc*0.0 + p_cond[0]

xfc=f_bivar(xfc,yfc,corr,mx=mean_x,my=mean_y,sx=std_x,sy=std_y)
Pf=xfc/np.sum(xfc) #normalize
j=-1
for i in range(1,N):
    Pf[i]=Pf[i]+Pf[i-1]
    if (j<0):
        if Pf[i] > samp_frac:
            j=i
p_cond[1]=yfc[j]

#get new marginal sample
xf=np.linspace(xmin_plt,xmax_plt,N)
yf=xf*0.0 + p_init[1]
yf=f_bivar_marg(xf,mp=mean_x,sp=std_x)

Pf=yf/np.sum(yf) #normalize
j=-1
for i in range(1,N):
    Pf[i]=Pf[i]+Pf[i-1]
    if (j<0):
        if Pf[i] > samp_frac:
            j=i
p_marg=np.array([0.0,0.0])
p_marg[0]=xf[j]

yfm=np.linspace(ymin_plt,ymax_plt,N)
xfm=yfm*0.0 + p_marg[0]

xfm=f_bivar(xfm,yfm,corr,mx=mean_x,my=mean_y,sx=std_x,sy=std_y)
Pf=xfm/np.sum(xfm) #normalize
j=-1
for i in range(1,N):
    Pf[i]=Pf[i]+Pf[i-1]
    if (j<0):
        if Pf[i] > samp_frac:
            j=i
p_marg[1]=yfc[j]


    
##################### v1 #######################
# Plot conditional vs marginal graphs version 1

#plot 68 and 95% CL
for i in range(len(CLs)):
    xplt,yplt = get_CL_line_rot(CLs[i],corr,mx=mean_x,my=mean_y,sx=std_x,sy=std_y)
    plt.plot(xplt,yplt,'k')

mu=np.array([mean_x,mean_y])
sigma=np.array([std_x,std_y])
    

#plot conditional and marginal (init)
xf=np.linspace(xmin_plt,xmax_plt,N)
yf=xf*0.0 + p_init[1]
yf=f_bivar(xf,yf,corr,mx=mean_x,my=mean_y,sx=std_x,sy=std_y)
yf=yf*fh/np.amax(yf)
plt.plot(xf,yf+p_init[1],color=col20[4],label='Conditional')

yf=xf*0.0 + p_init[1]
yf=f_bivar_marg(xf,mp=mean_x,sp=std_x)
yf=yf*fh/np.amax(yf)
plt.plot(xf,yf+p_init[1],color=col20[0],label='Marginal')

#plot start line
xline=np.array([xmin_plt,xmax_plt])
yline=np.array([p_init[1],p_init[1]])
plt.plot(xline,yline,'tab:grey')

#plot cond marg sampling
xfc=xfc*fh/np.amax(xfc)
plt.plot(xfc+p_cond[0],yfc,color=col20[5])

#plot cond marg sampling
xfm=xfm*fh/np.amax(xfm)
plt.plot(p_marg[0]-xfm,yfm,color=col20[1])

#plot new sampling lines
#marg
xline=np.array([p_marg[0],p_marg[0]])
yline=np.array([ymin_plt,ymax_plt])
plt.plot(xline,yline,'tab:grey')
#plot new sampling line
#cond
xline=np.array([p_cond[0],p_cond[0]])
yline=np.array([ymin_plt,ymax_plt])
plt.plot(xline,yline,'tab:grey')

#plot points
plt.plot(p_init[0],p_init[1],'o',color='tab:purple',ms=msize)
plt.plot(p_cond[0],p_cond[1],'o',color=col20[5],ms=msize)
plt.plot(p_marg[0],p_marg[1],'o',color=col20[1],ms=msize)
#plt.legend()
plt.xlim(xmin_plt, xmax_plt)
plt.ylim(ymin_plt, ymax_plt)
plt.gca().set_aspect('equal', adjustable='box')

plt.xticks([])
plt.yticks([])
plt.xlabel(r'$\beta$',fontsize=fs)
plt.ylabel(r'$a$',fontsize=fs)
plt.text(xmin_plt+0.2,ymax_plt-0.4, 'Conditional',fontsize=fs-1,color=col20[4])
plt.text(xmin_plt+0.2,ymax_plt-0.8, 'Marginal',fontsize=fs-1,color=col20[0])
plt.savefig('cond_vs_marg_graph_plot'+plttype,bbox_inches='tight',  pad_inches=0.02)
plt.clf()

##################### v2 #######################
# Plot conditional vs marginal graphs version 2
    
#plot 68 and 95% CL
for i in range(len(CLs)):
    xplt,yplt = get_CL_line_rot(CLs[i],corr,mx=mean_x,my=mean_y,sx=std_x,sy=std_y)
    plt.plot(xplt,yplt,'k')

mu=np.array([mean_x,mean_y])
sigma=np.array([std_x,std_y])

#plot conditional and marginal (init)
xf=np.linspace(xmin_plt,xmax_plt,N)
yf=xf*0.0 + p_init[1]
yf=f_bivar(xf,yf,corr,mx=mean_x,my=mean_y,sx=std_x,sy=std_y)
yf=yf*fh/np.amax(yf)
inds=np.where(yf > f_min)
plt.plot(xf[inds],yf[inds]+p_init[1],color=col20[4],label='Conditional')

yf=xf*0.0 + p_init[1]
yf=f_bivar_marg(xf,mp=mean_x,sp=std_x)
yf=yf*fh/np.amax(yf)
inds=np.where(yf > f_min)
plt.plot(xf[inds],yf[inds]+p_init[1],color=col20[0],label='Marginal')

#plot start line
xline=np.array([xmin_plt,xmax_plt])
yline=np.array([p_init[1],p_init[1]])
plt.plot(xline,yline,'--',color='tab:grey')

#plot cond marg sampling
xfc=xfc*fh/np.amax(xfc)
inds=np.where(xfc > f_min)
plt.plot(xfc[inds]+p_cond[0],yfc[inds],color=col20[5])

#plot cond marg sampling
xfm=xfm*fh/np.amax(xfm)
inds=np.where(xfm > f_min)
plt.plot(p_marg[0]-xfm[inds],yfm[inds],color=col20[1])

#plot new sampling lines
#marg
xline=np.array([p_marg[0],p_marg[0]])
yline=np.array([ymin_plt,ymax_plt])
plt.plot(xline,yline,'--',color='tab:grey')
#plot new sampling line
#cond
xline=np.array([p_cond[0],p_cond[0]])
yline=np.array([ymin_plt,ymax_plt])
plt.plot(xline,yline,'--',color='tab:grey')

#plot points
plt.plot(p_init[0],p_init[1],'o',color='tab:purple',ms=msize)
plt.plot(p_cond[0],p_cond[1],'o',color=col20[5],ms=msize)
plt.plot(p_marg[0],p_marg[1],'o',color=col20[1],ms=msize)
#plt.legend()
plt.xlim(xmin_plt, xmax_plt)
plt.ylim(ymin_plt, ymax_plt)
plt.gca().set_aspect('equal', adjustable='box')

plt.xticks([])
plt.yticks([])
plt.xlabel(r'$\beta$',fontsize=fs)
plt.ylabel(r'$a$',fontsize=fs)
plt.text(xmin_plt+0.2,ymax_plt-0.4, 'Conditional',fontsize=fs-1,color=col20[4])
plt.text(xmin_plt+0.2,ymax_plt-0.8, 'Marginal',fontsize=fs-1,color=col20[0])
plt.savefig('cond_vs_marg_graph_plot_v2'+plttype,bbox_inches='tight',  pad_inches=0.02)
plt.clf()

##################### v3 #######################
# Plot conditional vs marginal graphs version 3
    
#plot 68 and 95% CL
for i in range(len(CLs)):
    xplt,yplt = get_CL_line_rot(CLs[i],corr,mx=mean_x,my=mean_y,sx=std_x,sy=std_y)
    plt.plot(xplt,yplt,'k')

mu=np.array([mean_x,mean_y])
sigma=np.array([std_x,std_y])


#plot conditional and marginal (init)
xf=np.linspace(xmin_plt,xmax_plt,N)
yf=xf*0.0 + p_init[1]
yf=f_bivar(xf,yf,corr,mx=mean_x,my=mean_y,sx=std_x,sy=std_y)
yf=yf*fh/np.amax(yf)
inds=np.where(yf > f_min)
plt.plot(xf[inds],yf[inds]+p_init[1],color=col20[4],label='Conditional')

yf=xf*0.0 + p_init[1]
yf=f_bivar_marg(xf,mp=mean_x,sp=std_x)
yf=yf*fh/np.amax(yf)
inds=np.where(yf > f_min)
plt.plot(xf[inds],yf[inds]+p_init[1],color=col20[0],label='Marginal')

#plot start line
xline=np.array([xmin_plt,xmax_plt])
yline=np.array([p_init[1],p_init[1]])
plt.plot(xline,yline,'--',color='tab:grey')

#plot cond sampling
xfc=xfc*fh/np.amax(xfc)
inds=np.where(xfc > f_min)
plt.plot(p_cond[0]-xfc[inds],yfc[inds],color=col20[5])

#plot cond sampling from marg sample
xfm=xfm*fh/np.amax(xfm)
inds=np.where(xfm > f_min)
plt.plot(p_marg[0]-xfm[inds],yfm[inds],color=col20[1])

#plot new sampling lines
#marg
xline=np.array([p_marg[0],p_marg[0]])
yline=np.array([ymin_plt,ymax_plt])
plt.plot(xline,yline,'--',color='tab:grey')
#plot new sampling line
#cond
xline=np.array([p_cond[0],p_cond[0]])
yline=np.array([ymin_plt,ymax_plt])
plt.plot(xline,yline,'--',color='tab:grey')

#plot points
plt.plot(p_init[0],p_init[1],'o',color='tab:purple',ms=msize)
plt.plot(p_cond[0],p_cond[1],'o',color=col20[5],ms=msize)
plt.plot(p_marg[0],p_marg[1],'o',color=col20[1],ms=msize)
#plt.legend()
plt.xlim(xmin_plt, xmax_plt)
plt.ylim(ymin_plt, ymax_plt)
plt.gca().set_aspect('equal', adjustable='box')

plt.xticks([])
plt.yticks([])
plt.xlabel(r'$\beta$',fontsize=fs)
plt.ylabel(r'$a$',fontsize=fs)

plt.text(p_cond[0]+0.2,ymax_plt-0.25, r'$P(\beta \mid a_{\mathrm{init}})$',ha='left',va='top', fontsize=fs-1,color=col20[4])
plt.text(p_marg[0]+0.2,ymax_plt-0.2, r'$P(\beta)$',ha='left',va='top',fontsize=fs-1,color=col20[0])
plt.text(p_cond[0]+0.2,p_init[1]-0.8, r'$P(a \mid \beta_{\mathrm{cond}})$', fontsize=fs-1,color=col20[5])
plt.text(p_marg[0]-0.2,p_marg[1]+0.8, r'$P(a \mid \beta_{\mathrm{marg}})$', ha='right',fontsize=fs-1,color=col20[1])
plt.savefig('cond_vs_marg_graph_plot_v3'+plttype,bbox_inches='tight',  pad_inches=0.02)
plt.clf()


#create and plot Gibbs chains scatter plots

seeds=[12350]
Npoints=100

pc=np.zeros([2,Npoints+1])
pm=np.zeros([2,Npoints+1])
pc[:,0]=p_init.copy()
pm[:,0]=p_init.copy()

for seed in seeds:
    print('seed',seed)
    #first draw marginal samples
    np.random.seed(seed)
    #get new marginal samples
    for s in range(1,Npoints+1):
        xf=np.linspace(xmin_plt,xmax_plt,N)
        yf=xf*0.0 + pm[1,s-1]
        yf=f_bivar_marg(xf,mp=mean_x,sp=std_x)

        Pf=yf/np.sum(yf) #normalize
        j=-1
        ran=np.random.uniform(0,1.0)
        for i in range(1,N):
            Pf[i]=Pf[i]+Pf[i-1]
            if (j<0):
                if Pf[i] >= ran:
                    j=i
                    break
        pm[0,s]=xf[j]

        yf=np.linspace(ymin_plt,ymax_plt,N)
        xf=yf*0.0 + pm[0,s]

        xf=f_bivar(xf,yf,corr,mx=mean_x,my=mean_y,sx=std_x,sy=std_y)
        Pf=xf/np.sum(xf) #normalize
        j=-1
        ran=np.random.uniform(0,1.0)
        for i in range(1,N):
            Pf[i]=Pf[i]+Pf[i-1]
            if (j<0):
                if Pf[i] >= ran:
                    j=i
                    break
        pm[1,s]=yf[j]

    #secondly draw conditional samples
    np.random.seed(seed) #reset seed (draw same random variables 
    #get new conditional samples
    for s in range(1,Npoints+1):
        xf=np.linspace(xmin_plt,xmax_plt,N)
        yf=xf*0.0 + pc[1,s-1]
        yf=f_bivar(xf,yf,corr,mx=mean_x,my=mean_y,sx=std_x,sy=std_y)

        Pf=yf/np.sum(yf) #normalize
        j=-1
        ran=np.random.uniform(0,1.0)
        for i in range(1,N):
            Pf[i]=Pf[i]+Pf[i-1]
            if (j<0):
                if Pf[i] >= ran:
                    j=i
                    break
        pc[0,s]=xf[j]

        yf=np.linspace(ymin_plt,ymax_plt,N)
        xf=yf*0.0 + pc[0,s]

        xf=f_bivar(xf,yf,corr,mx=mean_x,my=mean_y,sx=std_x,sy=std_y)
        Pf=xf/np.sum(xf) #normalize
        j=-1
        ran=np.random.uniform(0,1.0)
        for i in range(1,N):
            Pf[i]=Pf[i]+Pf[i-1]
            if (j<0):
                if Pf[i] >= ran:
                    j=i
                    break
        pc[1,s]=yf[j]


    #plot the data
    #plot 68 and 95% CL
    for i in range(len(CLs)):
        xplt,yplt = get_CL_line_rot(CLs[i],corr,mx=mean_x,my=mean_y,sx=std_x,sy=std_y)
        plt.plot(xplt,yplt,'k')

    #plot points
    if (scat=='marg first'):
        plt.scatter(pm[0,1:],pm[1,1:],color=col20[0],s=msize-2)
        plt.scatter(pc[0,1:],pc[1,1:],color=col20[4],s=msize-2)
    else:
        plt.scatter(pc[0,1:],pc[1,1:],color=col20[4],s=msize-2)
        plt.scatter(pm[0,1:],pm[1,1:],color=col20[0],s=msize-2)
    plt.plot(p_init[0],p_init[1],'o',color='tab:purple',ms=msize)
    #plt.legend()
    plt.xlim(xmin_plt, xmax_plt)
    plt.ylim(ymin_plt, ymax_plt)
    plt.gca().set_aspect('equal', adjustable='box')

    plt.xticks([])
    plt.yticks([])
    plt.xlabel(r'$\beta$',fontsize=fs)
    plt.ylabel(r'$a$',fontsize=fs)
    plt.text(xmin_plt+0.2,ymax_plt-0.4, 'Conditional',fontsize=fs-1,color=col20[4])
    plt.text(xmin_plt+0.2,ymax_plt-0.8, 'Marginal',fontsize=fs-1,color=col20[0])
    plt.savefig('cond_vs_marg_scatt_plot_seed%i'%(seed)+plttype,bbox_inches='tight',  pad_inches=0.02)
    plt.clf()



exit()
Nseeds=10
seed=12345
Npoints=100

pc=np.zeros([2,Npoints+1])
pm=np.zeros([2,Npoints+1])
pc[:,0]=p_init.copy()
pm[:,0]=p_init.copy()

for k in range(Nseeds):
    print('seed',seed)
    #first draw marginal samples
    np.random.seed(seed)
    #get new marginal samples
    for s in range(1,Npoints+1):
        xf=np.linspace(xmin_plt,xmax_plt,N)
        yf=xf*0.0 + pm[1,s-1]
        yf=f_bivar_marg(xf,mp=mean_x,sp=std_x)

        Pf=yf/np.sum(yf) #normalize
        j=-1
        ran=np.random.uniform(0,1.0)
        for i in range(1,N):
            Pf[i]=Pf[i]+Pf[i-1]
            if (j<0):
                if Pf[i] >= ran:
                    j=i
                    break
        pm[0,s]=xf[j]

        yf=np.linspace(ymin_plt,ymax_plt,N)
        xf=yf*0.0 + pm[0,s]

        xf=f_bivar(xf,yf,corr,mx=mean_x,my=mean_y,sx=std_x,sy=std_y)
        Pf=xf/np.sum(xf) #normalize
        j=-1
        ran=np.random.uniform(0,1.0)
        for i in range(1,N):
            Pf[i]=Pf[i]+Pf[i-1]
            if (j<0):
                if Pf[i] >= ran:
                    j=i
                    break
        pm[1,s]=yf[j]

    #secondly draw conditional samples
    np.random.seed(seed) #reset seed (draw same random variables 
    #get new conditional samples
    for s in range(1,Npoints+1):
        xf=np.linspace(xmin_plt,xmax_plt,N)
        yf=xf*0.0 + pc[1,s-1]
        yf=f_bivar(xf,yf,corr,mx=mean_x,my=mean_y,sx=std_x,sy=std_y)

        Pf=yf/np.sum(yf) #normalize
        j=-1
        ran=np.random.uniform(0,1.0)
        for i in range(1,N):
            Pf[i]=Pf[i]+Pf[i-1]
            if (j<0):
                if Pf[i] >= ran:
                    j=i
                    break
        pc[0,s]=xf[j]

        yf=np.linspace(ymin_plt,ymax_plt,N)
        xf=yf*0.0 + pc[0,s]

        xf=f_bivar(xf,yf,corr,mx=mean_x,my=mean_y,sx=std_x,sy=std_y)
        Pf=xf/np.sum(xf) #normalize
        j=-1
        ran=np.random.uniform(0,1.0)
        for i in range(1,N):
            Pf[i]=Pf[i]+Pf[i-1]
            if (j<0):
                if Pf[i] >= ran:
                    j=i
                    break
        pc[1,s]=yf[j]


    #plot the data
    #plot 68 and 95% CL
    for i in range(len(CLs)):
        xplt,yplt = get_CL_line_rot(CLs[i],corr,mx=mean_x,my=mean_y,sx=std_x,sy=std_y)
        plt.plot(xplt,yplt,'k')

    #plot points
    if (scat=='marg first'):
        plt.scatter(pm[0,1:],pm[1,1:],color=col20[0],s=msize-2)
        plt.scatter(pc[0,1:],pc[1,1:],color=col20[4],s=msize-2)
    else:
        plt.scatter(pc[0,1:],pc[1,1:],color=col20[4],s=msize-2)
        plt.scatter(pm[0,1:],pm[1,1:],color=col20[0],s=msize-2)
    plt.plot(p_init[0],p_init[1],'o',color='tab:purple',ms=msize)
    #plt.legend()
    plt.xlim(xmin_plt, xmax_plt)
    plt.ylim(ymin_plt, ymax_plt)
    plt.gca().set_aspect('equal', adjustable='box')

    plt.xticks([])
    plt.yticks([])
    plt.xlabel(r'$\beta$',fontsize=fs)
    plt.ylabel(r'$a$',fontsize=fs)
    plt.text(xmin_plt+0.2,ymax_plt-0.4, 'Conditional',fontsize=fs-1,color=col20[4])
    plt.text(xmin_plt+0.2,ymax_plt-0.8, 'Marginal',fontsize=fs-1,color=col20[0])
    plt.savefig('cond_vs_marg_scatt_plot_seed%i'%(seed)+plttype,bbox_inches='tight',  pad_inches=0.02)
    plt.clf()

    seed+=1
