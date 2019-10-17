import numpy as np

def units(value, unit, order):
	# value - give the value of your parameter
	# unit - give the 'type' of the value's unit: "pico","nano", "micro", "milli", "centi", "deci", "zero"
	# ordr - give the order of the unit
	# example: 16 cm2 -> value = 16, unit = "centi", order = 2
	r = 0
	u = ["pico","nano", "micro", "milli", "centi", "deci", "zero"]
	v = [-12,       -9,      -6,      -3,      -2,     -1,      0]
	for i in range(0,len(u)):
		if unit == u[i]:
			r = value*10.**(v[i]*order)
	return(r)


def kplus(par):
	#provide parameters in the following order : aphat, dl, r, rt, sc
	# return value in units volume/time
	kon = par[0]; dl = par[1]; nr = par[2]; nrt = par[3]; sc = par[4];
	r = np.sqrt(nr*sc/(nrt*np.pi))
	kdl = 4*np.pi*dl*r
	return(kon*kdl/(kdl-kon*nr))


def kminus(par):
	#provide parameters in the following order : aphat, dl, r, rt, sc, am
	# return value in units 1/time
	kon = par[0]; dl = par[1]; nr = par[2]; nrt = par[3]; sc = par[4]; koff = par[5]
	r = np.sqrt(nr*sc/(nrt*np.pi))
	kdl = 4*np.pi*dl*r
	return(koff*(kdl+nr*kplus(par))/kdl)


#beta plus

def bplus(par):
	#provide parameters in the following order : aphat, dl, r, rt, sc, am, dr, b
	# return value in units area/time
	nrt = par[3]; sc = par[4]; dr = 2*par[6]; b = par[7]; hm=par[8];
	kp2D = kplus(par)/hm
	w = np.sqrt(sc/(nrt*np.pi))
	kdr = 2*np.pi*dr/(np.log(w/b))
	return(kp2D*kdr/(kdr+kp2D))	

#beta minus

def bminus(par):
	#provide parameters in the following order : aphat, dl, r, rt, sc, am, dr, b,
	# return value in units 1/time
	nrt = par[3]; sc = par[4]; dr = 2*par[6]; b = par[7]; hm = par[8];
	km = kminus(par)
	kp2D = kplus(par)/hm
	w = np.sqrt(sc/(nrt*np.pi))
	kdr = 2*np.pi*dr/(np.log(w/b))
	return(km*kdr/(kdr+kp2D))
	
	

