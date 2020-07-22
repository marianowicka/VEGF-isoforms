import numpy as np
import pandas as pd
from scipy.integrate import odeint
import statsmodels.api as sm
from scipy.special import logit
from scipy.special import expit
import scipy.stats
import math
from joblib import Parallel, delayed
import multiprocessing

# define data
df = pd.DataFrame()

df['time'] = [5.,5.,5.,5.,5.,5.,5.,15.,15.,15.,15.,15.,15.,15.,30.,30.,30.,30.,30.,30.,30.,60.,60.,60.,60.,60.,60.,60.]
df['conc'] = [0.,0.025,0.25,1.25,0.025,0.25,1.25,0.,0.025,0.25,1.25,0.025,0.25,1.25,0.,0.025,0.25,1.25,0.025,0.25,
              1.25,0.,0.025,0.25,1.25,0.025,0.25,1.25]
df['iso'] = ['none',165,165,165,121,121,121,'none',165,165,165,121,121,121,'none',165,165,165,121,121,121,'none',165,165,
             165,121,121,121]
df['set1'] = [0.0035575714,0.11132365,1.2740077358,1.,0.0765509492,0.1668638667,0.5170657593,0.1236339806,0.1926835614,
              0.7642575023,0.7775098059,0.0309980504,0.0966172999,0.5099644114,0.0162035263,0.0361912303,0.168570485,
              0.1705809873,0.0908178704,0.3711193794,0.3720809751,0.0513725437,0.024138434,0.1034273654,0.0766612266,
              0.0567229579,0.0780333965,0.064149561]
df['set2'] = [0.0290623029,0.2016807021,2.2291898407,1.,0.0505096377,0.1705836867,0.2320744114,0.0515953935,0.1353737104,
              0.7975269282,0.4734967272,0.0657081701,0.0878916565,0.1907239925,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,
              np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]
df['set3'] = [0.030616119,0.3291722246,0.8762754833,1.,0.0468036704,0.0738414136,0.2977818851,0.0167982909,0.2485932472,
              0.9202194614,0.8342452534,0.0625457161,0.1511175166,0.8811723915,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,
              np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]
df['set4'] = [0.0729483283,0.7092198582,1.4204660588,1.,0.3100303951,0.6828774063,0.5075987842,0.0729483283,0.2401215805,
              0.2077001013,0.2998986829,0.1438703141,0.3505572442,0.7375886525,0.0729483283,0.2279635258,0.3242147923,
              0.2907801418,0.1529888551,0.4498480243,0.7700101317,0.0729483283,0.1063829787,0.1610942249,0.3698074975,
              0.0759878419,0.1458966565,0.2887537994]
df['set5'] = [0.1094470046,0.465437788,1.9965437788,1.,0.1555299539,0.1785714286,0.418202765,0.1647465438,0.3006912442,
              1.1198156682,1.1140552995,0.1762672811,0.3675115207,0.8283410138,0.1347926267,0.2315668203,0.2580645161,
              0.2281105991,0.1900921659,0.2718894009,0.3306451613,0.1647465438,0.2096774194,0.2822580645,0.2511520737,
              0.1543778802,0.2085253456,0.2085253456]

df['set6'] = [0.0433121019,0.4968152866,0.7949044586,1.,0.0624203822,0.098089172,0.272611465,0.0840764331,0.778343949,
              0.6280254777,0.8535031847,0.1694267516,0.5006369427,0.7452229299,0.2038216561,0.7656050955,0.8828025478,
              0.8700636943,0.3630573248,0.5057324841,0.6777070064,0.2687898089,0.3834394904,0.5515923567,0.2420382166,
              0.0993630573,0.2203821656,0.1732484076]

def sort_data(DF, cl, iso):
    ind_con = DF['conc'] == cl
    ind_iso = DF['iso'] == iso
    index = ind_con & ind_iso
    newdf = DF.copy()
    newdf = newdf.loc[index, :]
    newdf = newdf.reset_index(drop=True)
    return newdf

def ode_model1(xini, t, par):
    alpp, alpm, betp, betm, krec, kdeg, kintM, kintP, f, nrt = par
    ksyn =kdeg*nrt/5.
    kintR =(krec+kdeg)/3.
    ls, rs, ms, ds, re, me, de = xini
    dls = -2*alpp*ls*rs+alpm*ms
    drs = ksyn -2*alpp*ls*rs+alpm*ms-betp*ms*rs+2*betm*ds-kintR*rs+krec*re
    dms = 2*alpp*ls*rs-alpm*ms-betp*ms*rs+2*betm*ds-kintM*ms
    dds = betp*ms*rs-2*betm*ds-kintP*ds
    dre = kintR*rs-krec*re-kdeg*re+2*f*betm*de+f*alpm*me
    dme = kintM*ms+2*f*betm*de-f*alpm*me
    dde = kintP*ds-2*f*betm*de
    
    dres = np.array([dls, drs, dms, dds, dre, dme, dde])
    return dres


def sim_data(xini, t, par, parNorm):
    simData = odeint(ode_model1, xini, t, args=(par,))
    xiniNorm = xini
    xiniNorm[0] = 1.25*10**5*6.022
    simDataNorm = odeint(ode_model1, xiniNorm, t, args=(parNorm,))
    tind = np.where(t==5*60)[0][0]
    res = (simData.T[3]+simData.T[6])/(simDataNorm.T[3][tind]+simDataNorm.T[6][tind])
    return res

def get_data_to_compare(DF, iso, cl):
    newDF = sort_data(DF, cl, iso)
    t5 = newDF.T[0][3:10].values
    t15 = newDF.T[1][3:10].values
    tint = newDF.T[2][3:10].values
    t30 = np.array([value for value in tint if not math.isnan(value)])
    tint = newDF.T[3][3:10].values
    t60 = np.array([value for value in tint if not math.isnan(value)])
    return t5, t15, t30, t60
 
def samp_prior():
    alpp165 =  10.**np.random.normal(-6,0.5)
    betp165 =  10.**np.random.normal(-2,0.5)
    alpm165 =  10.**np.random.normal(-6,2)
    betm165 =  10.**np.random.normal(-2,0.5)
    kdeg = 10.**np.random.normal(-3,1)
    krec = 10.**np.random.normal(-3,1)
    kintM165 = 10.**np.random.normal(-3,1)
    kintP165 = 10.**np.random.normal(-1.5, 0.25)
    f165 = 1.+10.**np.random.normal(2,1)
    nrt = 10.**np.random.normal(5,0.5)
    alpp121 =  10.**np.random.normal(-6,0.5)
    betp121 =  10.**np.random.normal(-2,0.5)
    alpm121 =  10.**np.random.normal(-6,2)
    betm121 =  10.**np.random.normal(-2,0.5)
    kintM121 = 10.**np.random.normal(-3,1)
    kintP121 = 10.**np.random.normal(-1.5,0.25)
    f121 = 1.+10.**np.random.normal(2,1)
    return alpp121,alpm121, betp121, betm121, kintM121, kintP121, f121, alpp165,alpm165, betp165, betm165,kdeg, krec, kintM165, kintP165, f165, nrt

conc = np.array([0.025,0.25,1.25]) #nM
lig = list(conc*(10**5)*6.022) #molecules
t_vec = np.arange(0,3601,1)

def dist(iteration,DF):
    print(iteration)
    while True:
        try:
            dis = []
            np.random.seed()
            alpp121,alpm121, betp121, betm121, kintM121, kintP121, f121, alpp165,alpm165, betp165, betm165,kdeg, krec, kintM165, kintP165, f165, nrt = samp_prior()
            par165 = [alpp165, alpm165, betp165, betm165, krec, kdeg, kintM165,kintP165, f165, nrt]
            par121 = [alpp121, alpm121, betp121, betm121, krec, kdeg, kintM121, kintP121, f121, nrt]
            parNorm = par165 
            t5 = np.where(t_vec == 5*60)[0][0]
            t15 = np.where(t_vec == 15*60)[0][0]
            t30 = np.where(t_vec == 30*60)[0][0]
            t60 = np.where(t_vec == 60*60)[0][0]  
            for i in range(0,3):
                cini = [lig[i], 0.6*nrt, 0., 0., 0.2*nrt, 0., 0.]
                res165 = sim_data(cini, t_vec, par165, parNorm)
                r1655 = res165[t5]
                r16515 = res165[t15]
                r16530 = res165[t30]
                r16560 = res165[t60]
                d1655, d16515, d16530, d16560 = get_data_to_compare(DF, 165, conc[i])
                if i < 2:
                    int_sum1655 = (r1655 - d1655)**2/(np.std(d1655)**2)
                elif i == 2:
                    int_sum1655 = 0.
                int_sum16515 = (r16515 - d16515)**2/(np.std(d16515)**2)
                int_sum16530 =  (r16530 - d16530)**2/(np.std(d16530)**2)
                int_sum16560 = (r16560 - d16560)**2/(np.std(d16560)**2)
                res121 = sim_data(cini, t_vec, par121, parNorm)
                r1215 = res121[t5]
                r12115 = res121[t15]
                r12130 = res121[t30]
                r12160 = res121[t60]
                d1215, d12115, d12130, d12160 = get_data_to_compare(DF, 121, conc[i])
                int_sum1215 = (r1215 - d1215)**2/(np.std(d1215)**2)
                int_sum12115 = (r12115 - d12115)**2/(np.std(d12115)**2)
                int_sum12130 = (r12130 - d12130)**2/(np.std(d12130)**2)
                int_sum12160 = (r12160 - d12160)**2/(np.std(d12160)**2)
                dis.append(np.sum(int_sum1215)+np.sum(int_sum12115)+np.sum(int_sum12130)+np.sum(int_sum12160)+np.sum(int_sum1655)+np.sum(int_sum16515)+np.sum(int_sum16530)+np.sum(int_sum16560))
            dis = np.array(dis)
            dis = np.sqrt(np.sum(dis))
        except:
            continue
        else:
            return np.array([dis, alpp165,alpm165,  alpp121, alpm121, betp165, betm165, betp121, betm121,  kdeg, krec, kintM165, kintM121,kintP165,kintP121, f165, f121,nrt])
            break


def run_par(DF, numIterations, nCores):
    Res = Parallel(n_jobs=nCores)(delayed(dist)(i,DF) for i in range(0,numIterations)) 
    return Res


Res = run_par(df, 1000000, 70)

 
Res = np.array(Res)
Rdf = pd.DataFrame()
Rdf['dis'] = Res.T[0]
Rdf['alpp165'] = Res.T[1]
Rdf['alpm165'] = Res.T[2]
Rdf['alpp121'] = Res.T[3]
Rdf['alpm121'] = Res.T[4]
Rdf['betp165'] = Res.T[5]
Rdf['betm165'] = Res.T[6]
Rdf['betp121'] = Res.T[7]
Rdf['betm121'] = Res.T[8]

Rdf['kdeg'] = Res.T[9]
Rdf['krec'] = Res.T[10]
Rdf['kintM165'] = Res.T[11]
Rdf['kintM121'] = Res.T[12]
Rdf['kintP165'] = Res.T[13]
Rdf['kintP121'] = Res.T[14]
Rdf['f165'] = Res.T[15]
Rdf['f121'] = Res.T[16]
Rdf['nrt'] = Res.T[17]

Rdf.to_csv('model1lastrun.csv')
