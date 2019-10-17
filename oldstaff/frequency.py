import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns#; sns.set(color_codes=True)
import matplotlib.patches as mpatches
import matplotlib.lines as mlines


delta = np.linspace(1.8,2.4,21)

h1 = np.genfromtxt('surface/post_sur.txt').T[0]
fh1=[]

h2 = np.genfromtxt('endosome/post_end.txt').T[0]
fh2=[]

for i in range(0,len(delta)):
	fh1.append(len(np.where(h1<delta[i])[0]))
	fh2.append(len(np.where(h2<delta[i])[0]))

fh1 = np.array(fh1)*1.
fh2 = np.array(fh2)*1.

pfh1 = fh1/(fh1+fh2)
pfh2 = fh2/(fh1+fh2)


plt.subplot(2,2,1)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
plt.plot(delta, pfh1, 'r', label=r'$p(H_1 | \delta = \delta^*)$')
plt.plot(delta, pfh2,'b', label=r'$p(H_2 | \delta = \delta^*)$')
plt.xlabel(r'$\delta^*$', fontsize = 20)
plt.title('Relative probability', fontsize = 20)
#plt.title(r'$p(H_i| \delta = \delta^*) = \frac{f(H_i|\delta = \delta^*)}{f(H_1| \delta = \delta^*)+f(H_2| \delta = \delta^*)}$' "\n", fontsize = 20)
plt.legend(loc = 1,fontsize = 20)

plt.subplot(2,2,2)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
plt.plot(delta, fh1, 'r', label = r'$f(H_1 | \delta = \delta^*)$')
plt.plot(delta, fh2, 'b', label = r'$f(H_2 | \delta = \delta^*)$')
plt.legend(loc = 1,fontsize = 20)
plt.title('Number of accepted parameters \n for the hypotheses $H_1$ and $H_2$ for threshold $\delta = \delta^*$', fontsize = 20)
#plt.title(r'$f(H_i | \delta = \delta^*)$ = number of accepted samples given $\delta^*$' "\n", fontsize = 20)
plt.xlabel(r'$\delta^*$', fontsize = 20)

'''
plt.subplot(2,2,3)

plt.plot(delta,fh2/fh1, color='k')
plt.fill_between(delta, [150]*len(delta),[200]*len(delta), color='r', alpha = 0.2)
plt.fill_between(delta, [150]*len(delta),[20]*len(delta), color='g', alpha = 0.2)
plt.fill_between(delta, [20]*len(delta),[3]*len(delta), color='b', alpha = 0.2)
plt.fill_between(delta, [3]*len(delta),[1]*len(delta), color='k', alpha = 0.5)

red_patch = mpatches.Patch(color='r', label='very strong',alpha=0.2)
green_patch = mpatches.Patch(color='g', label='strong',alpha=0.2)
blue_patch = mpatches.Patch(color='b', label='positive',alpha=0.2)
black_patch = mpatches.Patch(color='k', label='weak',alpha=0.5)
black_line = mlines.Line2D([], [], color='k',label='Bayes factor')

plt.legend(handles=[black_line, red_patch, green_patch, blue_patch, black_patch],loc='center left', bbox_to_anchor=(1, 0.6),fancybox=True, shadow=True, ncol=1,prop={'size':20})
'''
plt.show()

'''
 \\
	}

'''
