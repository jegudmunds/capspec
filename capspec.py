import os
import numpy as np
import datetime as dt
import viscosity as vs
from matplotlib.pyplot import *
from scipy.optimize import minimize
import matplotlib.cm as cm
import matplotlib.colors as colors

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap


def discrete_cmap(N, base_cmap=None):
    """Create an N-bin discrete colormap from the specified input map"""

    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string    , None, or a colormap instance:

    base = cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    return base.from_list(cmap_name, color_list, N)

class FlowMeas():
    '''
    The flow rate class
    '''

    def __init__(self,dp,dt,DeltaP,Tmt,Tcap,Tsft,
                 name=None,date=None,gas=None,eta=None,
                 Vsft=18.0,Vman=2.0,Troom=300.):

        self.dp = dp+1.0 # Pressure increase [torr]
        self.dt = dt # Integration time [min]
        self.DeltaP = DeltaP # Driving pressure [psia]
        self.Tmt = Tmt # Main tank [K]
        self.Tcap = Tcap # Temperature of the capillaries [K]
        self.Tsft = Tsft # Average temperature of the SFT [K]
        self.eta = eta # Viscosity in Poise
        self.gas = gas
        self.Q = None

        # Volume of the SFT and manifold in liters
        self.Vsft = Vsft #[L]
        self.Vman = Vman #[L]
        self.Troom = Troom
        self.R = 8.31 # J/K/mole

    def set_Q(self, a=0.5, b=1.0):
        '''
        Calculates the quality factor derived from this flow measurement
        '''
        self.get_eta(b)
        T = (a*self.Tsft+self.Troom)/(1+a)
        self.Q = 1e9*(self.dp/self.dt)*((self.Vsft/T+ \
          self.Vman/T)/self.R)*(self.eta/self.DeltaP)
        #self.Q = 1e9*(self.dp/self.dt)*((self.Vsft/self.Tsft+ \
        #  self.Vman/self.Troom)/self.R)*(self.eta/self.DeltaP)

        # OLDER VERSIONS
        #self.Q = 1e9*(self.dp/self.dt)*((self.Vsft/((a*self.Tsft+self.Troom)
        #/(1+a))+self.Vman/self.Troom)/self.R)*(self.eta/self.DeltaP)
        #self.Q = 1e9*(self.dp/self.dt)*((self.Vsft/self.Tsft+self.Vman/self.Troom)
        #/self.R)*(self.eta/self.DeltaP)

    def get_volumetric(self):
        atm1 = 760.0
        return 1e3*(self.Vsft+self.Vman)*self.dp/(60*self.dt)/atm1

    def get_Q(self, a=0.5, b=1.0):
        self.set_Q(a,b); return self.Q

    def get_eta(self, b=1.0):
        '''
        Calculates the viscosity given some gas type
        '''
        Tcap_mean = ((b*self.Tcap+self.Tmt)/(1.0+b))
        #Tcap_mean = self.Tcap #Tcap_mean = self.Tmt

        if self.gas == "He":
            self.eta = vs.eta_he(Tcap_mean)
        elif self.gas == "N2":
            self.eta = vs.eta_n2(Tcap_mean)
        elif self.gas==None and Tcap_mean > 250:
            self.eta = vs.eta_n2(Tcap_mean)
        elif self.gas==None:
            self.eta = vs.eta_he(Tcap_mean)
        else:
            raise ValueError('Need to define something')
        self.eta = float(self.eta)

def genlist():
    '''
    This function contains all flow measurements that we have performed.
    It's a silly way to keep track of things, but haven't got a more
    elegant solution. The function generates a list of FlowMeas objects that
    containt measurements from various runs in the past.

    Return:
    ftlist - the list of FlowMeas objects
    warmlist - the list of FlowMeas objects for wamr flow rates
    coldflow - the superfluid flow rate at for that corresponding run
    '''

    ## Run 10
    # This was the first run with current capillary assembly design
    # We used

    ## Run 11
    # Warm measurement
    r11_warm = FlowMeas(26.0,87,15.7,300,300,300,date=dt.datetime(2012,10,20))
    # 77K measurement
    # Verify temperatures
    r11_cold0 = FlowMeas(17.0,36,15.2,85,90,90,
        date=dt.datetime(2012,11,16),gas='He')

    ## Run 12
    # Warm measurement
    # Made up values for dp/dt. Only have that flow = 13.1e-2 cm^3/s
    r12_warm = FlowMeas(27.5,87,15.7,300,300,300,date=dt.datetime(2012,10,20))
    # 77K measurement
    r12_cold0 = FlowMeas(5.5,16,15.7,85,85,90,
        date=dt.datetime(2013,2,28),gas='He')


    # Why did we never get a flow rate measurement for

    ## Run 13
    # Warm measurement
    r13_warm = FlowMeas(46.0,169,15.7,300,300,300,date=dt.datetime(2013,6,6))
    # 77K measurement

    # We never really got cold enough to characterize capillary perforance
    # from 6/25 flow test
    #r13_cold0 = FlowMeas(14.5,63,14.7,45,35,date=dt.datetime(2013,6,25,12,42))

    # Once we got rid of the MT vent constriction we performed a series of
    # flow tests
    r13_cold1 = FlowMeas(14.5,63,14.7,11,42,35,
                         date=dt.datetime(2013,6,29,8,55),gas='He')
    r13_cold2 = FlowMeas(8.5,55,14.7,11,104,36,
                         date=dt.datetime(2013,6,29,9,56),gas='He')
    r13_cold3 = FlowMeas(11.0,120,14.7,12,156,37,
                         date=dt.datetime(2013,6,29,11,18),gas='He')
    r13_cold4 = FlowMeas(11.0,120,14.7,13.5,148,40,
                         date=dt.datetime(2013,6,29,13,41),gas='He')

    # We decided to warm up to see if increased flow rate through capillaries
    # could help with the burp events, we were successful in getting the flow
    # to go up but it didn't help with the burps
    r13_cold5 = FlowMeas(12.8,78,15.7,95,262,92,
                         date=dt.datetime(2013,7,14,18,25),gas='He')

    # Measured on Feb 27, based on logbook entry
    r14_warm1 =  FlowMeas(53.2,238,15.7,290,290,290,
                         date=dt.datetime(2014,2,27,21,59),gas='He')

    # Measured on Mar 27 2014
    r14_cold1 = FlowMeas(7.75,49,15.7,88,270,95,
                         date=dt.datetime(2014,3,27,17,13),gas='He')

    # Measured on Jul 2 2014
    #dp = 0.13 psi = 6.7 torr, dt = 34 min,
    r15_warm1 = FlowMeas(13.4,68,14.8,300,300,300,date=dt.datetime(2014,7,2),
                         gas='N2')

    r15_cold1 = FlowMeas(8.4,48,15.7,89,220,94,
                         date=dt.datetime(2014,7,6,13,12),gas='He')
    r15_cold2 = FlowMeas(5.51,27,15.86,90,221,100,
                         date=dt.datetime(2014,7,6,13,12),gas='He')

    # Measured on Nov 4 2014
    r16_warm1 = FlowMeas(19.5,101.5,15.1,295,295,295,
        date=dt.datetime(2014,11,4), gas='N2')
    # Measured on Nov 5 2014 (on the ice)
    r16_warm2 = FlowMeas(17.9,96,15.1,295,295,295,date=dt.datetime(2014,11,5),
                         gas='N2')

    r16_cold1 = FlowMeas(8.4,48,15.7,97,274,96,
                         date=dt.datetime(2014,11,25,14,45),gas='He')

    # Measured in the lab prior to installing in Lloro
    l2_warm1 = FlowMeas(103.43,11.88,16.69,295,295,295,
        date=dt.datetime(2017,4,1), gas='N2',Vsft=0.,Vman=0.580,Troom=295.)

    # Measurements in Lloro after closeup (warm)
    l2_warm2 = FlowMeas(116.87, 175, 14.7, 295, 295, 295, Vsft=6.464,
        Vman=0.5, gas='N2')

    l2_cold1 = FlowMeas(22.2, 50, 15.7, 290, 290, 290, Vsft=6.464,
        Vman=1.0, gas='He')

    # Warm (~300 K) flow measurements
    warm_run = np.array([11, 12, 13, 14, 15, 16, 2])
    warmlist = [r11_warm, r12_warm, r13_warm, r14_warm1,
        r15_warm1, r16_warm1, l2_warm1]

    # LN2 flow measurements
    ftlist = [r11_cold0, r13_cold1,
              r13_cold2, r13_cold3, r13_cold4, r13_cold5, r14_cold1]

    # Realized superfluid flow rates
    run = np.array([11.0, 13, 13, 13, 13, 13, 14])
    coldflows = np.array([2.45, 0.8, 0.8, 0.8, 0.8, 2.2, 1.7])

    volumetric = [ft.get_volumetric() for ft in ftlist]
    warm_volumetric = [ft.get_volumetric() for ft in warmlist]

    if len(ftlist) != len(coldflows):
        print "Make sure flow rate values are as many as flow measurements"

    return run, ftlist, volumetric, coldflows, warm_run, \
        warmlist, warm_volumetric

def f2min(x):

    run, f2list, volumetric, coldflows, warm_run, \
        warmlist, warm_volumetric = genlist()

    Q = [ft.get_Q(x[0], x[1]) for ft in ftlist]
    Qf = Q/coldflows

    return np.sum(np.abs(Qf-1))/float(len(Qf))


### ----------------------------------------------------------------------------
###
### Making some plots for analysis
###
### ----------------------------------------------------------------------------


print('\nQuality factor is a unitless parameter that is derived from LN2 flow rates ')
print('The parameters used to calculate the quality factor have been tweaked so that ')
print('the quality factor roughly maps to cold flow rate in SLPM. Note that Q/f from')
print('runs 11, 13, and 14 is 1.2 on average.\n')
print('\n-----\n')

# Obtaining data from previous cryo runs
run, ftlist, volumetric, coldflows, warm_run, warmlist,\
    warm_volumetric = genlist()

# Calculating warm quality factor for kicks
Q_warm = [ft.get_Q(a=0.51,b=10.0) for ft in warmlist];

#Q = [ft.get_Q(a=0.51,b=10.0) for ft in ftlist]; Q = np.array(Q)
Q = [ft.get_Q(a=0.51,b=5.0) for ft in ftlist]; Q = np.array(Q)
Qf = Q/coldflows

print('Run:\t |' + '\t|'.join(['{:d}'.format(int(r)) for r in run]))
print('Q:\t |' + '\t|'.join(['{:5.2f}'.format(q) for q in Q]))
print('Q/f:\t |' + '\t|'.join(['{:5.2f}'.format(qf) for qf in Qf]))
print('Vdot:\t |' + '\t|'.join(['{:5.3f}'.format(v) for v in volumetric]) +
    ' [cm^3/s]')
print('Flow:\t |' + '\t|'.join(['{:5.2f}'.format(v) for v in coldflows]) +
    ' [SLPM]')

print('\nRun:\t\t |' + '\t|'.join(['{:d}'.format(int(r)) for r in warm_run]))
print('Vdot (warm):\t |' + '\t|'.join(['{:5.3f}'.format(v) \
    for v in warm_volumetric])+' [cm^3/s]')

Qf_mean = np.mean(Qf); Qf_std = np.std(Qf)
print('\nMean quality factor: {:5.2f} '.format(Qf_mean))
print('STD of quality factor: {:5.2f} \n'.format(Qf_std))

#r14_cold1 = FlowMeas(6.72,50,15.7,88,270,95,
#                     date=dt.datetime(2014,3,27,17,13),gas='He')
r14_cold1 = FlowMeas(7.75,49,15.7,88,270,95,
                     date=dt.datetime(2014,3,27,17,13),gas='He')

Q14 = r14_cold1.get_Q(a=0.47,b=5.0)
r15_cold1 = FlowMeas(8.4,48,15.7,89,220,94,
                     date=dt.datetime(2014,7,6,13,12),gas='He')
r15_cold2 = FlowMeas(5.51,27,15.86,90,221,100,
                     date=dt.datetime(2014,7,6,13,12),gas='He')

r16_cold1 = FlowMeas(8.4,48,15.7,97,274,96,
                     date=dt.datetime(2014,11,25,14,45),gas='He')

r16_cold1 = FlowMeas(7.75,50,15.86,97,280,96,
                     date=dt.datetime(2014,11,25,14,45),gas='He')

# Measured in the lab prior to installing in Lloro
l2_warm1 = FlowMeas(103.43,11.88,16.69,295,295,295,date=dt.datetime(2017,4,1),
                     gas='N2',Vsft=0.,Vman=0.580,Troom=295.)

# Measurements in Lloro after closeup (warm)
l2_warm2 = FlowMeas(116.87, 175, 14.7, 295, 295, 295, Vsft=6.464,
    Vman=0.5, gas='N2')

l2_cold1 = FlowMeas(22.2, 50, 15.7, 290, 290, 290, Vsft=6.464,
    Vman=1.0, gas='He')

Q15_1 = r15_cold1.get_Q(a=0.47,b=5.0)
Q15_2 = r15_cold2.get_Q(a=0.47,b=5.0)
Q16_1 = r16_cold1.get_Q(a=0.47,b=5.0)
Q2_1 = l2_cold1.get_Q(a=0.47,b=5.0)

#Q = np.append(Q,Q14)
#Qf = np.append(Qf,Q14/1.67)

print('Quality factor from last few runs:')
print('  Q14:  \t{:5.2f}'.format(Q14))
print('  Q15_1:  \t{:5.2f}'.format(Q15_1))
print('  Q15_2:  \t{:5.2f}'.format(Q15_2))
print('  Q16:  \t{:5.2f}'.format(Q16_1))
print('  Q2 (cold):  \t{:5.2f}\n\n'.format(Q2_1))

###
### Plotting the Qualify factor
###
plot(1+np.arange(len(Q)),Q,marker='*',lw=0,ms=12,markeredgecolor='black',
     markerfacecolor='None')
locs, labels = xticks()
xlim([0, len(Q)+1])
labels2use = ['R11','R13a','R13a','R13a','R13a','R13b','R14']
ax = gca()
ax.set_xticks(locs)
ax.set_xticklabels(labels2use)
xlabel('Run')

# The quality factor is essentially the predicted flow
ylabel('Predicted Flow [SLPM]')
ylim([0.5, 3.5])
savefig('img/predicted_flow.png',figsize=(8, 6), dpi=200)
close()

###
### Plotting the ratio of predicted to realized flow
###
plot(1+np.arange(len(Qf)),Qf,marker='*',lw=0,ms=12,markeredgecolor='black',
     markerfacecolor='None')
locs, labels = xticks()
xlim([0, len(Q)+1])
labels2use = ['R11','R13a','R13a','R13a','R13a','R13b','R14']
ax = gca()
ax.set_xticks(locs)
ax.set_xticklabels(labels2use)
xlabel('Run'); ylabel('$f_\mathrm{pred}/f_\mathrm{rel}$')
savefig('img/flow_ratio.png',figsize=(8, 6), dpi=200)
close()

###
### Plotting the predicted versus realized flow
###
plot(2.5/3.1*Q, coldflows, marker='*',
    lw=0, ms=12,markeredgecolor='black',
     markerfacecolor='None')
locs, labels = xticks()
ax = gca()
xlim([0.0, 3.0]); ylim([0.0, 3.0])
xlabel('Predicted flow [SLPM]');
ylabel('Realized flow [SLPM]')
savefig('img/flow_vs_pred.png',figsize=(8, 6), dpi=200)
close()

run_norm = (run-np.min(run))/(np.max(run)-np.min(run))
scatter(2.5/3.1*Q, coldflows, c=run_norm, s=100, cmap='summer',
    edgecolor='black')
colorbar()

xlim([0.0, 3.0]); ylim([0.0, 3.0])
axis('equal')
xlabel('Predicted flow [SLPM]');
ylabel('Realized flow [SLPM]')
savefig('img/flow_vs_pred2.png',figsize=(8, 6), dpi=200)
close()

###
### Plotting the Qf variation as a function of alpha
###
N = 1000
alpha = np.linspace(0.2,10.0,N)
val = np.zeros(np.shape(alpha))
for idx, a in enumerate(alpha):
    val[idx] = f2min([a, 10.0])

plot(alpha, val,lw=0,marker='.')
xlabel('a')
ylabel('Residual')
savefig('img/average_variation.png')
close()

nlev = 5
levels = np.linspace(0.01,1.0,nlev)
cmap_tmp = cm.get_cmap('summer',len(levels)-1)
cmap2use = truncate_colormap(cmap_tmp,0.0,0.50)

run = np.array([9, 10, 11, 12, 13, 14, 15])
coldflows = np.array([0.45, 3.1, 2.5, 0.65, 2.2, 1.7, 3.14])
impedance = np.array([240, 36, 47, 41, 49, 59, 55.5])

run = np.array([10, 11, 13, 14, 15])
coldflows = np.array([3.1, 2.5, 2.2, 1.7, 3.14])
impedance = np.array([36, 47, 49, 59, 55.5])

r0 = int(np.min(run))
rd = int(np.max(run)-np.min(run))+1

figure(figsize=(6,3))


z = np.polyfit(impedance[:-1], coldflows[:-1], 1)
f = np.poly1d(z)
impi = np.linspace(0.9*np.min(impedance), 1.1*np.max(impedance), 50)
coldi = f(impi)

print 'Best fit linear model:'
print f

scatter(impedance, coldflows, c=run, s=50,
    cmap=discrete_cmap(rd, 'RdYlBu'), edgecolor='black')

plot(impi, coldi, color='black', label='Fit')
cbar = colorbar(ticks=range(r0,r0+rd))
clim(r0-0.5, r0+rd-0.5)

cbar.set_ticklabels(['Run 10','Run 11','Run 12','Run 13','Run 14','Run 2'])

xlabel('Flow impedance [cm$^{-3}$]')
ylabel('Capillary cold flow rate [SLPM]')
savefig('img/impedance_vs_cold.png',dpi=200,
    bbox_inches='tight')
close()

###
### Plotting the sfd
###

labels2use = ['R11','R13a','R13a','R13a','R13a','R13b','R14']

run, ftlist, volumetric, coldflows, warm_run, warmlist, \
    warm_volumetric = genlist()

# Calculating warm quality factor for kicks
#Q = [ft.get_Q(a=0.51,b=10.0) for ft in ftlist]; Q = np.array(Q)
Q = [ft.get_Q(a=1., b=37.) for ft in ftlist]; Q = np.array(Q)
Qf = Q/coldflows

run = np.array([11, 13, 13, 13, 13, 13, 14])

r0 = int(np.min(run))
rd = int(np.max(run)-np.min(run))+1

z = np.polyfit(Q, coldflows, 1)
f = np.poly1d(z)
Qi = np.linspace(0.9*np.min(Q), 1.1*np.max(Q), 50)
coldi = f(Qi)

print 'Best fit linear model:'
print f

figure(figsize=(6,3))
scatter(Q, coldflows, c=run, s=50,
    cmap=discrete_cmap(rd, 'RdYlBu'), edgecolor='black')

plot(Qi, coldi, color='black', label='Fit')

cbar = colorbar(ticks=range(r0,r0+rd))
clim(r0-0.5, r0+rd-0.5)

cbar.set_ticklabels(['Run 11','Run 12','Run 13','Run 14'])

xlabel('Quality factor')
ylabel('Capillary cold flow rate [SLPM]')
savefig('img/q_vs_cold.png', dpi=200, bbox_inches='tight')
close()

###
### Some crude minimization routine
###

a0 = 0.3; b0 = 1.0
bnds = ((0.2,5.0), (0.2,10.0))
#print len(bnds)
x0 = [a0, b0]
#res = minimize(f2min,(a0,b0),method='L-BFGS-B',bounds=bnds)
res = minimize(f2min, x0, method='L-BFGS-B', bounds=bnds)
print res
print res.x
print f2min(res.x)/len(Q)

## Mapping the 2D parameter space
N = 80
a1, a2 = 0.1, 10.
#b1, b2 = 39.0, 40.
b1, b2 = 1., 10.

alpha = np.linspace(a1, a2, N)
beta = np.linspace(b1, b2, N)
resid = np.zeros((N,N))

for i in xrange(N):
    for j in xrange(N):
        resid[i,j] = f2min([alpha[i], beta[j]])

#imshow(resid, extent=[a1, a2, b1, b2], cmap='RdYlBu',
#    aspect='auto',vmin=resid.min(), vmax=0.2)
imshow(resid, extent=[a1, a2, b1, b2], cmap='RdYlBu',
    aspect='auto')

colorbar()
xlabel('alpha')
ylabel('beta')
savefig('img/brute_force.png')
close()



