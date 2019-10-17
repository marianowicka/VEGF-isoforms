import numpy as np
from scipy import arange
import matplotlib.pyplot as plt


import seaborn as sns#; sns.set(color_codes=True)

conc = np.array([0.025,0.25,1.25]) #nM
lig = list(conc*(10**5)*6.022) #molecules

#erk data

mean165_0025 = np.array([0.3738466957,0.7398699765,0.3635385957,0.2146770428])
mean165_025  = np.array([0.9624716313,1.5622983355,1.3475041738,0.4158901363])
mean165_125  = np.array([1,1.565011909,1.1390334686,0.2538133392])
Mean165 = np.array([mean165_0025,mean165_025,mean165_125[1:len(mean165_125)]])

var165_0025 = np.array([0.0140604898,0.0260077877,0.0952533245,0.0169104176])
var165_025  = np.array([0.0284406991,1.5919126561,0.5981518517,0.0082594898])
var165_125  = np.array([0,1.0557334733,0.211088122,0.0002276298])
Var165 = np.array([var165_0025,var165_025,var165_125[1:len(var165_125)]])

mean121_0025 = np.array([0.1048327379,0.2048746446,0.1988403698,0.0972506629])
mean121_025  = np.array([0.1398482568,0.860565077,0.3075486139,0.1623867909])
mean121_125  = np.array([0.4747238464,1.4566145324,0.6686451279,0.1468588427])
Mean121 = np.array([mean121_0025,mean121_025,mean121_125])

var121_0025 = np.array([0.0015428588,0.0093699883,0.0158658862,0.0010963897])
var121_025  = np.array([0.0286968007,0.4767832349,0.0265586378,0.0046858272])
var121_125  = np.array([0.0645179203,1.4083926413,0.2027224829,0.0036057395])
Var121 = np.array([var121_0025,var121_025,var121_125])

std165_0025 = np.array([0.1185769361,0.1612693019,0.3086313732,0.1300400613])
std165_025 = np.array([0.1686437046,1.2617102108,0.7734027746,0.0908817353])
std165_125 = np.array([0,1.0274889164,0.4594432748,0.0150874056])
Std165 = np.array([std165_0025, std165_025, std165_125[1:len(std165_125)]])

std121_0025 = np.array([0.0392792416,0.0967986999,0.1259598595,0.0331117761])
std121_025 = np.array([0.0535693949,0.6904949202,0.1629682109,0.068453102])
std121_125 = np.array([0.2540037801,1.1867571956,0.4502471354,0.0600478104])
Std121 = np.array([std121_0025, std121_025, std121_125])



conc = np.array([0.025,0.25,1.25]) #nM
lig = list(conc*(10**5)*6.022) #molecules
time = [5*60,15*60,30*60,60*60] #sec
tick1 = ['  0',10,20,' 30',40,50,60]


tpoints = arange(0,3601,20)

per5_end = np.genfromtxt('per5_end.txt')
per50_end = np.genfromtxt('per50_end.txt')
per95_end = np.genfromtxt('per95_end.txt')

#per5_end = np.genfromtxt('per5_sur_no_delay.txt')
#per50_end = np.genfromtxt('per50_sur_no_delay.txt')
#per95_end = np.genfromtxt('per95_sur_no_delay.txt')

#per5_sur = np.genfromtxt('per5_end_no_delay.txt')
#per50_sur = np.genfromtxt('per50_end_no_delay.txt')
#per95_sur = np.genfromtxt('per95_end_no_delay.txt')

per5_sur = np.genfromtxt('per5_sur.txt')
per50_sur = np.genfromtxt('per50_sur.txt')
per95_sur = np.genfromtxt('per95_sur.txt')

fs = 20
fs2 = 15

plt.grid()
ax=plt.subplot(2,3,1)
plt.ylim(0,7)
plt.plot(tpoints,per50_sur.T[0],'r-', label = r'$\overline{sim}(t,c_l,iso)$ with $H_1$')
plt.plot(tpoints,per50_end.T[0],'b-', label = r'$\overline{sim}(t,c_l,iso)$ with $H_2$')
#plt.plot(tpoints,per50_end.T[0],'b-', label = 'hypothesis H1 \nfrom surface \n'r'$\tau = 0$')
plt.fill_between(tpoints,per5_end.T[0],per95_end.T[0],color = 'b', alpha = 0.2)
#plt.plot(tpoints,per50_sur.T[0],'r-', label = 'hypothesis H2 \nfrom endosome \n'r'$\tau = 0$')
plt.fill_between(tpoints,per5_sur.T[0],per95_sur.T[0],color = 'r', alpha = 0.2)
plt.errorbar([5*60, 15*60, 30*60, 60*60], mean165_0025, yerr=std165_0025, barsabove=True, ls="none", marker="o", mfc="k", color="k")
plt.title(r'VEGF-A$_{165}$, $c_l = 0.025$ $nM$', fontsize = fs)
#plt.ylabel(r'$\frac{S(t)}{S(t=5min)| c_l = 1.25 nM \ of \ VEGF-A_{165}}$', fontsize = 14)
plt.xticks(np.arange(0, 3601, 600))
a=ax.get_xticks().tolist()
a = tick1	
ax.set_xticklabels(a, fontsize = fs2)
plt.yticks(fontsize = fs2)
plt.legend(loc='center left', bbox_to_anchor=(0, -1.7),fancybox=True, shadow=True, ncol=2,prop={'size':30})


ax=plt.subplot(2,3,2)
plt.ylim(0,7)
plt.plot(tpoints,per50_end.T[1],'b-')
plt.fill_between(tpoints,per5_end.T[1],per95_end.T[1],color = 'b', alpha = 0.2)
plt.plot(tpoints,per50_sur.T[1],'r-')
plt.fill_between(tpoints,per5_sur.T[1],per95_sur.T[1],color = 'r', alpha = 0.2)
plt.errorbar([5*60, 15*60, 30*60, 60*60], mean165_025, yerr=std165_025, barsabove=True, ls="none", marker="o", mfc="k", color="k")
plt.title(r'VEGF-A$_{165}$, $c_l=0.25$ $nM$', fontsize = fs)
plt.xticks(np.arange(0, 3601, 600))
a=ax.get_xticks().tolist()
a = tick1	
ax.set_xticklabels(a, fontsize = fs2)
plt.yticks(fontsize = fs2)

ax=plt.subplot(2,3,3)
plt.ylim(0,7)
plt.plot(tpoints,per50_end.T[2],'b-')
plt.fill_between(tpoints,per5_end.T[2],per95_end.T[2],color = 'b', alpha = 0.2)
plt.plot(tpoints,per50_sur.T[2],'r-')
plt.fill_between(tpoints,per5_sur.T[2],per95_sur.T[2],color = 'r', alpha = 0.2)
plt.errorbar([5*60, 15*60, 30*60, 60*60], mean165_125, yerr=std165_125, barsabove=True, ls="none", marker="o", mfc="k", color="k")
plt.title(r'VEGF-A$_{165}$, $c_l=1.25$ $nM$', fontsize = fs)
plt.xticks(np.arange(0, 3601, 600))
a=ax.get_xticks().tolist()
a = tick1	
ax.set_xticklabels(a, fontsize = fs2)
plt.yticks(fontsize = fs2)

ax=plt.subplot(2,3,4)
plt.ylim(0,7)
plt.plot(tpoints,per50_end.T[3],'b-')
plt.fill_between(tpoints,per5_end.T[3],per95_end.T[3],color = 'b', alpha = 0.2)
plt.plot(tpoints,per50_sur.T[3],'r-')
plt.fill_between(tpoints,per5_sur.T[3],per95_sur.T[3],color = 'r', alpha = 0.2)
plt.errorbar([5*60, 15*60, 30*60, 60*60], mean121_0025, yerr=std121_0025, barsabove=True, ls="none", marker="o", mfc="k", color="k")
plt.title(r'VEGF-A$_{121}$, $c_l=0.025$ $nM$', fontsize = fs)
#plt.ylabel(r'$\frac{S(t)}{S(t=5min)| c_l = 1.25 nM \ of \ VEGF-A_{165}}$', fontsize = 14)
plt.xticks(np.arange(0, 3601, 600))
a=ax.get_xticks().tolist()
a = tick1	
ax.set_xticklabels(a, fontsize = fs2)
plt.yticks(fontsize = fs2)

ax=plt.subplot(2,3,5)
plt.ylim(0,7)
plt.plot(tpoints,per50_end.T[4],'b-')
plt.fill_between(tpoints,per5_end.T[4],per95_end.T[4],color = 'b', alpha = 0.2)
plt.plot(tpoints,per50_sur.T[4],'r-')
plt.fill_between(tpoints,per5_sur.T[4],per95_sur.T[4],color = 'r', alpha = 0.2)
plt.errorbar([5*60, 15*60, 30*60, 60*60], mean121_025, yerr=std121_025, barsabove=True, ls="none", marker="o", mfc="k", color="k")
plt.title(r'VEGF-A$_{121}$, $c_l=0.25$ $nM$', fontsize = fs)
plt.xticks(np.arange(0, 3601, 600))
a=ax.get_xticks().tolist()
a = tick1	
ax.set_xticklabels(a, fontsize = fs2)
plt.yticks(fontsize = fs2)

ax=plt.subplot(2,3,6)
plt.ylim(0,7)
plt.plot(tpoints,per50_end.T[5],'b-')
plt.fill_between(tpoints,per5_end.T[5],per95_end.T[5],color = 'b', alpha = 0.2)
plt.plot(tpoints,per50_sur.T[5],'r-')
plt.fill_between(tpoints,per5_sur.T[5],per95_sur.T[5],color = 'r', alpha = 0.2)
plt.errorbar([5*60, 15*60, 30*60, 60*60], mean121_125, yerr=std121_125, barsabove=True, ls="none", marker="o", mfc="k", color="k")
plt.title(r'VEGF-A$_{121}$, $c_l=1.25$ $nM$', fontsize = fs)
plt.xticks(np.arange(0, 3601, 600))
a=ax.get_xticks().tolist()
a = tick1	
ax.set_xticklabels(a, fontsize = fs2)
plt.yticks(fontsize = fs2)
text = r'$Time$ $(min)$'
plt.text(1800,-1.5,text, fontsize = fs)
plt.subplots_adjust(left=0.05, bottom=0.20, right=0.7, top=0.95, wspace=0.2, hspace=0.35)
plt.show()
