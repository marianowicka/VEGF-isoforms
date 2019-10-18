import numpy as np
from scipy.integrate import odeint
from scipy import arange
import PyDDE.pydde as p



# fixed rates for 165
ap165,am165,bp165,bm165,kint165,krec165,kdegR165,kintC165,bmf165,amf165,kdegC165,ksyn165,krecC165, NRT = np.genfromtxt('fixed_pams_165.txt')

# fixed rates for 121
ap121,am121,bp121,bm121,kint121,krec121,kdegR121,kintC121,bmf121,amf121,kdegC121,ksyn121,krecC121, NRT = np.genfromtxt('fixed_pams_121.txt')
#initial conditions
nrs = NRT*0.6
nre = NRT*0.2

Rs0 = nrs
Ms = 0
Ps0 = 0
Re0 = nre
Me0 = 0
Pe0 = 0
S0 = 0

conc = np.array([0.025,0.25,1.25]) #nM
lig = list(conc*(10**5)*6.022) #molecules

print nrs, nre, lig

def odegradErk(s, c, t):
	ap = c[0]; am=c[1]; bp= c[2]; bm=c[3]; kint=c[4]; krec=c[5]; kdegR=c[6]; kintC=c[7]; bmf=c[8];
	amf = c[9]; kdegC= c[10]; ksyn=c[11]; krecC= c[12]; mus=c[13]; lambdas=c[14]; kappas=c[15]; tau = c[16]

	L = s[0]
	R = s[1]
	M = s[2]
	P = s[3]
	Re = s[4]
	Me = s[5]
	Pe = s[6]
	S = s[7]

	dLdt = -2*ap*L*R+am*M
	dRsdt = -2*ap*L*R+am*M-bp*M*R+2*bm*P-kint*R+krec*Re+ksyn
	dMsdt = 2*ap*L*R-am*M-bp*M*R+2*bm*P-kint*M
	dPsdt = bp*M*R-2*bm*P-kintC*P
	dRedt = kint*R-krec*Re-kdegR*Re+2*bmf*Pe+amf*Me
	dMedt = kint*M+2*bmf*Pe-amf*Me
	dPedt = kintC*P-2*bmf*Pe
	dSdt = 	-mus*(S)+lambdas*(Pe)/(Pe+kappas)
	return np.array([dLdt, dRsdt, dMsdt, dPsdt, dRedt, dMedt, dPedt, dSdt])	



def ddegradErk(s, c, t):
	ap = c[0]; am=c[1]; bp= c[2]; bm=c[3]; kint=c[4]; krec=c[5]; kdegR=c[6]; kintC=c[7]; bmf=c[8];
	amf = c[9]; kdegC= c[10]; ksyn=c[11]; krecC= c[12]; mus=c[13]; lambdas=c[14]; kappas=c[15];
	tau = c[16]

	L = s[0]
	R = s[1]
	M = s[2]
	P = s[3]
	Re = s[4]
	Me = s[5]
	Pe = s[6]
	S = s[7]

	Pelag = 0.0
    	if (t>tau):
        	Pelag = p.pastvalue(6,t-tau,0)

	dLdt = -2*ap*L*R+am*M
	dRsdt = -2*ap*L*R+am*M-bp*M*R+2*bm*P-kint*R+krec*Re+ksyn
	dMsdt = 2*ap*L*R-am*M-bp*M*R+2*bm*P-kint*M
	dPsdt = bp*M*R-2*bm*P-kintC*P
	dRedt = kint*R-krec*Re-kdegR*Re+2*bmf*Pe+amf*Me
	dMedt = kint*M+2*bmf*Pe-amf*Me
	dPedt = kintC*P-2*bmf*Pe
	dSdt = 	-mus*(S)+lambdas*(Pelag)/(Pelag+kappas)

	return np.array([dLdt, dRsdt, dMsdt, dPsdt, dRedt, dMedt, dPedt, dSdt])		

    
def ddesthistErk(g, s, c, t):
    return (s, g)	

def find_sol(ap,am,bp,bm,kint,krec,kdegR,kintC,bmf,amf,kdegC,ksyn,krecC,mus, lambdas, kappas,taus,L0,Rs0,Ms,Ps0,Re0,Me0,Pe0,S0):
	#  Duration of the simulation
	endTimeErk = 3601 
	# The coefficients (constants) in the equations 
	odeconsErk = np.array([ap,am,bp,bm,kint,krec,kdegR,kintC,bmf,amf,kdegC,ksyn,krecC,mus, lambdas, kappas,taus])
	#  Initial conditions
	odeistErk = np.array([L0,Rs0,Ms,Ps0,Re0,Me0,Pe0,S0])
	#  Create an ode/dde object
	ode_egErk = p.dde()

	ode_egErk.initproblem(no_vars=8, no_cons=17, nlag=0, nsw=0,t0=0.0, t1=endTimeErk,initstate=odeistErk, c=odeconsErk, otimes=arange(0.0,endTimeErk,1),grad=odegradErk)

	odestscErk = np.array([L0,Rs0,0,0,Re0,0,0,0])

	ode_egErk.initsolver(tol=1*10**(-8), hbsize=10**4,dt=0.1,statescale=odestscErk)

	ode_egErk.solve()

	dde_egErk = p.dde()

	ddeistErk = odeistErk
	ddeconsErk = odeconsErk
	ddestscErk = odestscErk

	dde_egErk.dde(y=ddeistErk, times=arange(0.0, endTimeErk, 1), func=ddegradErk, parms=ddeconsErk, tol=0.000005, dt=0.1, hbsize=10**4, nlag=1, ssc=ddestscErk)

	return(dde_egErk)


conc = np.array([0.025,0.25,1.25]) #nM
lig = list(conc*(10**5)*6.022) #molecules
time = [5*60,15*60,30*60,60*60] #sec
#tode5 = np.linspace(0,300,300)
#tode60 = np.linspace(0,3600,3600)

rates = np.genfromtxt('endosome/post_end.txt')

#arange(a,b,c)  from a to b every c
tpoints = arange(0,3601,20)

for i in xrange(0,len(rates)):
	mus, lambdas, kappas, taus = rates[i][1:5]
	sim0025iso165 = find_sol(ap165,am165,bp165,bm165,kint165,krec165,kdegR165,kintC165,bmf165,amf165,kdegC165,ksyn165,krecC165,mus, lambdas, kappas,taus,lig[0],Rs0,Ms,Ps0,Re0,Me0,Pe0,S0)
	sim025iso165 = find_sol(ap165,am165,bp165,bm165,kint165,krec165,kdegR165,kintC165,bmf165,amf165,kdegC165,ksyn165,krecC165,mus, lambdas, kappas,taus,lig[1],Rs0,Ms,Ps0,Re0,Me0,Pe0,S0)
	sim125iso165 = find_sol(ap165,am165,bp165,bm165,kint165,krec165,kdegR165,kintC165,bmf165,amf165,kdegC165,ksyn165,krecC165,mus, lambdas, kappas,taus,lig[2],Rs0,Ms,Ps0,Re0,Me0,Pe0,S0)

	sim0025iso121 = find_sol(ap121,am121,bp121,bm121,kint121,krec121,kdegR121,kintC121,bmf121,amf121,kdegC121,ksyn121,krecC121,mus, lambdas, kappas,taus,lig[0],Rs0,Ms,Ps0,Re0,Me0,Pe0,S0)
	sim025iso121 = find_sol(ap121,am121,bp121,bm121,kint121,krec121,kdegR121,kintC121,bmf121,amf121,kdegC121,ksyn121,krecC121,mus, lambdas, kappas,taus,lig[1],Rs0,Ms,Ps0,Re0,Me0,Pe0,S0)
	sim125iso121 = find_sol(ap121,am121,bp121,bm121,kint121,krec121,kdegR121,kintC121,bmf121,amf121,kdegC121,ksyn121,krecC121,mus, lambdas, kappas,taus,lig[2],Rs0,Ms,Ps0,Re0,Me0,Pe0,S0)

	control = sim125iso165.data[:,8][299]

	a_0025iso165 = sim0025iso165.data[:,8]/control
	a_025iso165 = sim025iso165.data[:,8]/control
	a_125iso165 = sim125iso165.data[:,8]/control

	a_0025iso121 = sim0025iso121.data[:,8]/control
	a_025iso121 = sim025iso121.data[:,8]/control
	a_125iso121 = sim125iso121.data[:,8]/control

	for j in xrange(0,len(tpoints)):
		print 	a_0025iso165[tpoints[j]], a_025iso165[tpoints[j]], a_125iso165[tpoints[j]], a_0025iso121[tpoints[j]],a_025iso121[tpoints[j]],a_125iso121[tpoints[j]]


