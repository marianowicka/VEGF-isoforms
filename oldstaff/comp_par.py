import numpy as np
import def_par as dp

def par_comp(R2tot,iso):
	na = 6.022*10**(23)			#mol-1
	am = 10.**(-3) 	#s-1  Interactions of VEGF isoforms with VEGFR-1, VEGFR-2, andneuropilin in vivo: a computational model of human skeletal muscle
	kd = 100.				#pM Differential binding of VEGF isoforms to VEGF receptor 2 in the presence of neuropilin-1: a computational model, 
							# ref Interactions of VEGF isoforms with VEGFR-1, VEGFR-2, and neuropilin in vivo: a computational model of human skeletal muscle
	kd = dp.units(kd, "pico",1)
	rt= R2tot
	per = 1 		#percent of the cell 
	r = rt*per
	ap = am/kd				#M-1 s-1
	ap = ap/na					#dm3 s-1
	ap = dp.units(ap,"deci",3)
	if iso == 165:
		dl = 5.2*10**(-7)				#cm2 s-1 Differential binding of VEGF isoforms to VEGF receptor 2 in the presence of neuropilin-1: a computational model
	if iso == 121:
		dl = 5.8*10**(-7)
	dl = dp.units(dl, "centi", 2)
	sc = 1000.					#um2  --------------Popel 2006
	sc = dp.units(sc, "micro",2)  

	dr = 10.**(-10)				#cm2 s-1
	dr = dp.units(dr, "centi",2)
	 
	b = 0.5						#nm
	b = dp.units(b, "nano",1)

	#divide by the membrane thickness
	hm = 0.5					#um Interactions of VEGF isoforms with VEGFR-1, VEGFR-2, and neuropilin in vivo: a computational model of human skeletal muscle
								#ref Capillary changes in skeletal muscle of patients with essential hypertension
	hm = dp.units(hm, "micro", 1)
	
	h = 1.						#mm
	h = dp.units(h, "milli", 1)

	par = [ap, dl, r, rt, sc,am, dr, b, hm]
	kp = dp.kplus(par)
	km = dp.kminus(par)
	bp = dp.bplus(par)
	bm = dp.bminus(par)

	scf = r*sc/rt

	AP = ap/(scf*h)
	AM = am
	BP = dp.bplus(par)/scf
	BM = dp.bminus(par)
	return(np.array([AP,AM,BP,BM]))


