import sys
import os
import numpy as np
from numpy import mean
np.seterr(divide='ignore', invalid='ignore')
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

import julian
import datetime

from statistics import stdev
from statistics import mean

from matplotlib import rc
rc('text', usetex=True)
rc('font', family='serif',size=14)

DP_Leo_Schwope_2002 = open("Schwope_2002_TDB.dat",'r').readlines()
N_dpleo_Schwope_2002 = len(DP_Leo_Schwope_2002)

#Read datat
Cycle = []
T_obs = []
T_obs_err = []
#Please change the input file
for line in open("Schwope_2002_TDB.dat"):
    li=line.strip()
    if not li.startswith("#"):
        Cycle.append(float(li.split(" ")[0]))
        T_obs.append(float(li.split(" ")[1]))
        T_obs_err.append(float(li.split(" ")[2]))


#Arrays
LT_a = [i for i in range(len(DP_Leo_Schwope_2002))]

Local_Time = []
Flux = []
for i in range(len(DP_Leo_Schwope_2002)):
    LT = julian.from_jd(T_obs[i], fmt='mjd')
    LT_a[i] = LT
    Local_Time.append(LT_a[i])
    print(i, T_obs[i], LT_a[i])


#New ephemeris
#Linear phemeris equation(From equation 1) convert from TT to TDB
#T0_bjd = 2456717.36832561 - 2400000
#T0_bjd_err = 0.000000001
#P0_day = 0.0623628426
#P0_day_err = 0.0000000006

#Ephemeris of Schwope(2002)
#2448773.215071(18) + 0.06236283691(70)xE
#T0_bjd = 2448773.215071 - 2400000
#T0_bjd_err = 0.000018
#P0_day = 0.06236283691
#P0_day_err = 0.00000000070

#Ephemeris of Beuermann(2011)
#2454914.8322920(20) + 0.06236285648(90)xE
#T0_bjd = 2454914.8322920 - 2400000
#T0_bjd_err = 0.0000020
#P0_day = 0.06236285648
#P0_day_err = 0.00000000090


#New ephemeris
T0_bjd = 2454914.8322920 - 2400000
T0_bjd_err = 0.000000001
P0_day = 0.0623628426
P0_day_err = 0.0000000006


#Arrays
BJD_time_a = [i for i in range(N_dpleo_Schwope_2002)]
Delta_aT = [i for i in range(N_dpleo_Schwope_2002)]
Delta_aT_err = [i for i in range(N_dpleo_Schwope_2002)]
E_af = [i for i in range(N_dpleo_Schwope_2002)] #float number
E_ak = [i for i in range(N_dpleo_Schwope_2002)] #integer number
E_aj = [i for i in range(N_dpleo_Schwope_2002)] #integer number
P_aE = [i for i in range(N_dpleo_Schwope_2002)]
P_err_aE = [i for i in range(N_dpleo_Schwope_2002)]
T_aC_linear = [i for i in range(N_dpleo_Schwope_2002)]
T_aO_linear = [i for i in range(N_dpleo_Schwope_2002)]
P_aver_a = [i for i in range(N_dpleo_Schwope_2002)]
P_aver_std_a = [i for i in range(N_dpleo_Schwope_2002)]

delta_tdb_tt = 0.0013/(24*60*60)
OC_cal = []
#print ('-----------------------------------------------------------------------------')
#print ('Cycle \t\t T_O \t   T_C \t\t BJD - 2450000 \t OC_lin OC_err_Lin OC_occ')
print('No. \t BJD_time \t Cycle \t T_O_linear \t T_C_linear \t OC_s \t\t OC_s_err')
#print ('-----------------------------------------------------------------------------')
for i in range (0,N_dpleo_Schwope_2002):
    BJD_time = np.array(T_obs)+ delta_tdb_tt + 2400000
    BJD_time_a[i] = BJD_time
    Delta_T = np.array(T_obs) - np.array(T0_bjd)
    Delta_aT[i] = Delta_T #arrays
    Delta_T_err = np.sqrt((np.array(T_obs_err)/np.array(T_obs))**2 + (np.array(T0_bjd_err)/np.array(T0_bjd))**2)
    E_k = Cycle
    E_ak[i] = E_k #arrays
    #    print (Delta_T_err[i])
    E_f = Delta_T / P0_day                      #Calculate cycle with float number
    ##    print (E_f)                                 #print cycle with float number
    E_af[i] = E_f #arrays
    E_j = np.round(Delta_T / P0_day)           #Calculate cycle with integer number
##print (Delta_T)
    if  E_j[i] != 0:
        P_E_day = Delta_T[i] / E_j[i]
        P_aE[i] = P_E_day
        P_E_err_day = np.abs((np.array(T_obs_err[i]) - np.array(T0_bjd_err)) / E_j[i])
        P_err_aE[i] = P_E_err_day
        T_O_linear = T0_bjd + P_aE[i]*E_j[i]               #Linear
        T_aO_linear[i] = T_O_linear #arrays
    else:
        E_k[i] = 1
        P_E_day = Delta_T[i] / E_k[i]
#        print (P_E_day)
        P_aE[i] = P_E_day
        P_E_err_day = np.abs((np.array(T_obs_err[i]) - np.array(T0_bjd_err)) / E_k[i])
        P_err_aE[i] = P_E_err_day
        T_O_linear = T0_bjd + P_aE[i]*E_k[i]               #Linear
        T_aO_linear[i] = T_O_linear #arrays
#    print ('%0.6f' %(T_O))
#print (E_j)                                #print cycle with integer number
    if  E_j[i] != 0:
        P_E_day = Delta_T[i] / E_j[i]
        P_aE[i] = P_E_day
        P_E_err_day = np.abs((np.array(T_obs_err[i]) - np.array(T0_bjd_err)) / E_j[i])
#        print (P_E_err_day)
        P_err_aE[i] = P_E_err_day
        T_C_linear = T0_bjd + P0_day*E_j[i]              #Linear
        T_aC_linear[i] = T_C_linear #arrays
#    print (T_O, T_C)
        OC = np.array(T_O_linear) - np.array(T_C_linear)
        OC_s = (np.array(T_O_linear) - np.array(T_C_linear))*24*60*60
        OC_err = np.abs(np.sqrt((np.array(P_err_aE[i])**2 + (np.array(P0_day_err)**2))) * np.array(E_j[i]))
#        print (OC_err)
        OC_s_err = OC_err*24*60*60
    else:
        P_E_day = Delta_T[i] / E_k[i]
        P_aE[i] = P_E_day
        P_E_err_day = np.abs((np.array(T_obs_err[i]) - np.array(T0_bjd_err)) / E_k[i])
        P_err_aE[i] = P_E_err_day
        T_C_linear = T0_bjd + P0_day*E_j[i]              #Linear
        T_aC_linear[i] = T_C_linear #arrays
#    print (T_O, T_C)
        OC = np.array(T_O_linear) - np.array(T_C_linear)
        OC_s = (np.array(T_O_linear) - np.array(T_C_linear))*24*60*60
        OC_err = np.abs(np.sqrt((np.array(P_err_aE[i])**2)) *np.array(E_k[i]))
#        print (OC_err)
        OC_s_err = OC_err*24*60*60
    print ('%0.0f\t%0.6f\t%0.0f\t%0.6f\t%0.6f\t%0.6f\t%0.6f' %(i, BJD_time[i], Cycle[i], T_O_linear, T_C_linear, OC_s, OC_s_err))
    OC_cal.append('%0.6f\t%0.0f\t%0.6f\t%0.6f\t%0.6f\t%0.6f' %(BJD_time[i], Cycle[i], T_O_linear, T_C_linear, OC_s, OC_s_err))

#P_aver = mean(P_aE[i])
#P_aver_a[i] = P_aver
#P_aver_std = np.std(P_aE[i])
#P_aver_std = mean(P_err_aE[i])
#P_aver_std_a[i] = P_aver_std
#print('%0.11f %0.11f' %(P_aver, P_aver_std))

rerults = OC_cal
f = open("Schwope_2002_TDB.out", 'w')
for i in range(len(rerults)):
    f.write(str(rerults[i])+ '\n')
f.close()

#Plot O-C vs BJD
Input_Schwope_2002  = 'Schwope_2002_TDB.out'
Data_Schwope_2002   = np.genfromtxt(Input_Schwope_2002)

N = 2450000
BJD_time_Schwope_2002 = Data_Schwope_2002[:,0] - N
Cycle_Schwope_2002 = Data_Schwope_2002[:,1]
T_O_linear_Schwope_2002 = Data_Schwope_2002[:,2]
T_C_linear_Schwope_2002 = Data_Schwope_2002[:,3]
OC_s_Schwope_2002 = Data_Schwope_2002[:,4]
OC_s_err_Schwope_2002 = Data_Schwope_2002[:,5]

 ##Plotgraph
fig=plt.figure(figsize=(15, 5), tight_layout=True)
plt.tick_params(direction='in', which='both', bottom='on',top='on', right = 'on')

x1 = min(BJD_time_Schwope_2002)
x2 = max(BJD_time_Schwope_2002)
#plt.errorbar(BJD_time, OC_s, yerr=OC_s_err, fmt='o', color='limegreen')
plt.errorbar(BJD_time_Schwope_2002, OC_s_Schwope_2002, yerr= OC_s_err_Schwope_2002, fmt='o', markersize=8, color='green',
                    ecolor='lightgray', label = 'Schwope et al. (2002)')

#Schwope_2002
plt.text(BJD_time_Schwope_2002[0], 250, '1979')
plt.text(BJD_time_Schwope_2002[3], 250, '1981')
plt.text(BJD_time_Schwope_2002[13], 250, '1984')
plt.text(BJD_time_Schwope_2002[17], 250, '1985')
plt.text(BJD_time_Schwope_2002[19], 250, '1991')
plt.text(BJD_time_Schwope_2002[20]+100, 250, '1992')
plt.text(BJD_time_Schwope_2002[22], 250, '1993')
plt.text(BJD_time_Schwope_2002[30], 250, '2000')
plt.text(BJD_time_Schwope_2002[32], 250, '2002')

#plt.hlines(y= 0, xmin= x1, xmax= x2, colors='k', linestyles='dotted')
#plt.xlim(x1,x2)
#plt.ylim(-20,20)
plt.xlabel('BJD - 2450000')
plt.ylabel('O-C (sec)')
plt.legend()
plt.grid(linestyle='dotted')
######plt.title('O-C diagram: DP Leo')
#output_filename = os.path.splitext(__file__)[0] + '.png'
#plt.savefig(output_filename, dpi=1000)
plt.savefig("OC_Schwope_2002_rev1.png", dpi=1000)
plt.show()


#Input file: Beuermann_2011
DP_Leo_Beuermann_2011 = open("Beuermann_2011_TDB.dat",'r').readlines()
N_dpleo_Beuermann_2011 = len(DP_Leo_Beuermann_2011)

#Read datat
Cycle = []
T_obs = []
T_obs_err = []
#Please change the input file
for line in open("Beuermann_2011_TDB.dat"):
    li=line.strip()
    if not li.startswith("#"):
        Cycle.append(float(li.split(" ")[0]))
        T_obs.append(float(li.split(" ")[1]))
        T_obs_err.append(float(li.split(" ")[2]))
        
#Arrays
LT_a = [i for i in range(len(DP_Leo_Beuermann_2011))]

Local_Time = []
Flux = []
for i in range(len(DP_Leo_Beuermann_2011)):
    LT = julian.from_jd(T_obs[i], fmt='mjd')
    LT_a[i] = LT
    Local_Time.append(LT_a[i])
    print(i, T_obs[i], LT_a[i])

#Arrays
BJD_time_a = [i for i in range(N_dpleo_Beuermann_2011)]
Delta_aT = [i for i in range(N_dpleo_Beuermann_2011)]
Delta_aT_err = [i for i in range(N_dpleo_Beuermann_2011)]
E_af = [i for i in range(N_dpleo_Beuermann_2011)] #float number
E_ak = [i for i in range(N_dpleo_Beuermann_2011)] #integer number
E_aj = [i for i in range(N_dpleo_Beuermann_2011)] #integer number
P_aE = [i for i in range(N_dpleo_Beuermann_2011)]
P_err_aE = [i for i in range(N_dpleo_Beuermann_2011)]
T_aC_linear = [i for i in range(N_dpleo_Beuermann_2011)]
T_aO_linear = [i for i in range(N_dpleo_Beuermann_2011)]
P_aver_a = [i for i in range(N_dpleo_Beuermann_2011)]
P_aver_std_a = [i for i in range(N_dpleo_Beuermann_2011)]

delta_tdb_tt = 0.0013/(24*60*60)
OC_cal = []
#print ('-----------------------------------------------------------------------------')
#print ('Cycle \t\t T_O \t   T_C \t\t BJD - 2450000 \t OC_lin OC_err_Lin OC_occ')
print('No. \t BJD_time \t Cycle \t T_O_linear \t T_C_linear \t OC_s \t\t OC_s_err')
#print ('-----------------------------------------------------------------------------')
for i in range (0,N_dpleo_Beuermann_2011):
    BJD_time = np.array(T_obs)+ delta_tdb_tt + 2400000
    BJD_time_a[i] = BJD_time
    Delta_T = np.array(T_obs) - np.array(T0_bjd)
    Delta_aT[i] = Delta_T #arrays
    Delta_T_err = np.sqrt((np.array(T_obs_err)/np.array(T_obs))**2 + (np.array(T0_bjd_err)/np.array(T0_bjd))**2)
    E_k = Cycle
    E_ak[i] = E_k #arrays
    #    print (Delta_T_err[i])
    E_f = Delta_T / P0_day                      #Calculate cycle with float number
    ##    print (E_f)                                 #print cycle with float number
    E_af[i] = E_f #arrays
    E_j = np.round(Delta_T / P0_day)           #Calculate cycle with integer number
##print (Delta_T)
    if  E_j[i] != 0:
        P_E_day = Delta_T[i] / E_j[i]
        P_aE[i] = P_E_day
        P_E_err_day = np.abs((np.array(T_obs_err[i]) - np.array(T0_bjd_err)) / E_j[i])
        P_err_aE[i] = P_E_err_day
        T_O_linear = T0_bjd + P_aE[i]*E_j[i]               #Linear
        T_aO_linear[i] = T_O_linear #arrays
    else:
        E_k[i] = 2
        P_E_day = Delta_T[i] / E_k[i]
#        print (P_E_day)
        P_aE[i] = P_E_day
        P_E_err_day = np.abs((np.array(T_obs_err[i]) - np.array(T0_bjd_err)) / E_k[i])
        P_err_aE[i] = P_E_err_day
        T_O_linear = T0_bjd + P_aE[i]*E_k[i]               #Linear
        T_aO_linear[i] = T_O_linear #arrays
#    print ('%0.6f' %(T_O))
#print (E_j)                                #print cycle with integer number
    if  E_j[i] != 0:
        P_E_day = Delta_T[i] / E_j[i]
        P_aE[i] = P_E_day
        P_E_err_day = np.abs((np.array(T_obs_err[i]) - np.array(T0_bjd_err)) / E_j[i])
#        print (P_E_err_day)
        P_err_aE[i] = P_E_err_day
        T_C_linear = T0_bjd + P0_day*E_j[i]              #Linear
        T_aC_linear[i] = T_C_linear #arrays
#    print (T_O, T_C)
        OC = np.array(T_O_linear) - np.array(T_C_linear)
        OC_s = (np.array(T_O_linear) - np.array(T_C_linear))*24*60*60
        OC_err = np.abs(np.sqrt((np.array(P_err_aE[i])**2 + (np.array(P0_day_err)**2))) * np.array(E_j[i]))
#        print (OC_err)
        OC_s_err = OC_err*24*60*60
    else:
        P_E_day = Delta_T[i] / E_k[i]
        P_aE[i] = P_E_day
        P_E_err_day = np.abs((np.array(T_obs_err[i]) - np.array(T0_bjd_err)) / E_k[i])
        P_err_aE[i] = P_E_err_day
        T_C_linear = T0_bjd + P0_day*E_j[i]              #Linear
        T_aC_linear[i] = T_C_linear #arrays
#    print (T_O, T_C)
        OC = np.array(T_O_linear) - np.array(T_C_linear)
        OC_s = (np.array(T_O_linear) - np.array(T_C_linear))*24*60*60
        OC_err = np.abs(np.sqrt((np.array(P_err_aE[i])**2)) *np.array(E_k[i]))
#        print (OC_err)
        OC_s_err = OC_err*24*60*60
    print ('%0.0f\t%0.6f\t%0.0f\t%0.6f\t%0.6f\t%0.6f\t%0.6f' %(i, BJD_time[i], Cycle[i], T_O_linear, T_C_linear, OC_s, OC_s_err))
    OC_cal.append('%0.6f\t%0.0f\t%0.6f\t%0.6f\t%0.6f\t%0.6f' %(BJD_time[i], Cycle[i], T_O_linear, T_C_linear, OC_s, OC_s_err))

#P_aver = mean(P_aE[i])
#P_aver_a[i] = P_aver
#P_aver_std = np.std(P_aE[i])
#P_aver_std = mean(P_err_aE[i])
#P_aver_std_a[i] = P_aver_std
#print('%0.11f %0.11f' %(P_aver, P_aver_std))

rerults = OC_cal
f = open("Beuermann_2011_TDB.out", 'w')
for i in range(len(rerults)):
    f.write(str(rerults[i])+ '\n')
f.close()


#Plot O-C vs BJD
Input_Beuermann_2011  = 'Beuermann_2011_TDB.out'
Data_Beuermann_2011   = np.genfromtxt(Input_Beuermann_2011)

N = 2450000
BJD_time_Beuermann_2011 = Data_Beuermann_2011[:,0] - N
Cycle_Beuermann_2011 = Data_Beuermann_2011[:,1]
T_O_linear_Beuermann_2011 = Data_Beuermann_2011[:,2]
T_C_linear_Beuermann_2011 = Data_Beuermann_2011[:,3]
OC_s_Beuermann_2011= Data_Beuermann_2011[:,4]
OC_s_err_Beuermann_2011 = Data_Beuermann_2011[:,5]

##Plotgraph
fig=plt.figure(figsize=(15, 5), tight_layout=True)
plt.tick_params(direction='in', which='both', bottom='on',top='on', right = 'on')

x1 = min(BJD_time_Beuermann_2011)
x2 = max(BJD_time_Beuermann_2011)
plt.errorbar(BJD_time_Beuermann_2011, OC_s_Beuermann_2011, yerr= OC_s_err_Beuermann_2011, fmt='o', markersize=8, color='blue',
                    ecolor='lightgray', label = 'Beuermann et al. (2011)')

#Beuermann_2011
plt.text(BJD_time_Beuermann_2011[0], -25, '2009')
plt.text(BJD_time_Beuermann_2011[42], -25, '2010')


#plt.hlines(y= 0, xmin= x1, xmax= x2, colors='k', linestyles='dotted')
#plt.xlim(x1,x2)
#plt.ylim(-20,20)
plt.xlabel('BJD - 2450000')
plt.ylabel('O-C (sec)')
plt.legend()
plt.grid(linestyle='dotted')
######plt.title('O-C diagram: DP Leo')
#output_filename = os.path.splitext(__file__)[0] + '.png'
#plt.savefig(output_filename, dpi=1000)
plt.savefig("OC_Beuermann_2011_rev1.png", dpi=1000)
plt.show()

#Input file: Beuermann_2011
DP_Leo_Kittipong_2020 = open("Kittipong_2020.dat",'r').readlines()
N_dpleo_Kittipong_2020 = len(DP_Leo_Kittipong_2020)

#Read datat
Cycle = []
T_obs = []
T_obs_err = []
#Please change the input file
for line in open("Kittipong_2020.dat"):
    li=line.strip()
    if not li.startswith("#"):
        Cycle.append(float(li.split(" ")[0]))
        T_obs.append(float(li.split(" ")[1]))
        T_obs_err.append(float(li.split(" ")[2]))
        

#Arrays
LT_a = [i for i in range(len(DP_Leo_Kittipong_2020))]

Local_Time = []
Flux = []
for i in range(len(DP_Leo_Kittipong_2020)):
    LT = julian.from_jd(T_obs[i], fmt='mjd')
    LT_a[i] = LT
    Local_Time.append(LT_a[i])
    print(i, T_obs[i], LT_a[i])
    
    
#Arrays
BJD_time_a = [i for i in range(N_dpleo_Kittipong_2020)]
Delta_aT = [i for i in range(N_dpleo_Kittipong_2020)]
Delta_aT_err = [i for i in range(N_dpleo_Kittipong_2020)]
E_af = [i for i in range(N_dpleo_Kittipong_2020)] #float number
E_ak = [i for i in range(N_dpleo_Kittipong_2020)] #integer number
E_aj = [i for i in range(N_dpleo_Kittipong_2020)] #integer number
P_aE = [i for i in range(N_dpleo_Kittipong_2020)]
P_err_aE = [i for i in range(N_dpleo_Kittipong_2020)]
T_aC_linear = [i for i in range(N_dpleo_Kittipong_2020)]
T_aO_linear = [i for i in range(N_dpleo_Kittipong_2020)]
P_aver_a = [i for i in range(N_dpleo_Kittipong_2020)]
P_aver_std_a = [i for i in range(N_dpleo_Kittipong_2020)]


OC_cal = []
#print ('-----------------------------------------------------------------------------')
#print ('Cycle \t\t T_O \t   T_C \t\t BJD - 2450000 \t OC_lin OC_err_Lin OC_occ')
print('No. \t BJD_time \t Cycle \t T_O_linear \t T_C_linear \t OC_s \t\t OC_s_err')
#print ('-----------------------------------------------------------------------------')
for i in range (0,N_dpleo_Kittipong_2020):
    BJD_time = np.array(T_obs) + 2400000
    BJD_time_a[i] = BJD_time
    Delta_T = np.array(T_obs) - np.array(T0_bjd)
    Delta_aT[i] = Delta_T #arrays
    Delta_T_err = np.sqrt((np.array(T_obs_err)/np.array(T_obs))**2 + (np.array(T0_bjd_err)/np.array(T0_bjd))**2)
    E_k = Cycle
    E_ak[i] = E_k #arrays
    #    print (Delta_T_err[i])
    E_f = Delta_T / P0_day                      #Calculate cycle with float number
    ##    print (E_f)                                 #print cycle with float number
    E_af[i] = E_f #arrays
    E_j = np.round(Delta_T / P0_day)           #Calculate cycle with integer number
##print (Delta_T)
    if  E_j[i] != 0:
        P_E_day = Delta_T[i] / E_j[i]
        P_aE[i] = P_E_day
        P_E_err_day = np.abs((np.array(T_obs_err[i]) - np.array(T0_bjd_err)) / E_j[i])
        P_err_aE[i] = P_E_err_day
        T_O_linear = T0_bjd + P_aE[i]*E_j[i]               #Linear
        T_aO_linear[i] = T_O_linear #arrays
    else:
        E_k[i] = 1
        P_E_day = Delta_T[i] / E_k[i]
#        print (P_E_day)
        P_aE[i] = P_E_day
        P_E_err_day = np.abs((np.array(T_obs_err[i]) - np.array(T0_bjd_err)) / E_k[i])
        P_err_aE[i] = P_E_err_day
        T_O_linear = T0_bjd + P_aE[i]*E_k[i]               #Linear
        T_aO_linear[i] = T_O_linear #arrays
#    print ('%0.6f' %(T_O))
#print (E_j)                                #print cycle with integer number
    if  E_j[i] != 0:
        P_E_day = Delta_T[i] / E_j[i]
        P_aE[i] = P_E_day
        P_E_err_day = np.abs((np.array(T_obs_err[i]) - np.array(T0_bjd_err)) / E_j[i])
#        print (P_E_err_day)
        P_err_aE[i] = P_E_err_day
        T_C_linear = T0_bjd + P0_day*E_j[i]              #Linear
        T_aC_linear[i] = T_C_linear #arrays
#    print (T_O, T_C)
        OC = np.array(T_O_linear) - np.array(T_C_linear)
        OC_s = (np.array(T_O_linear) - np.array(T_C_linear))*24*60*60
        OC_err = np.abs(np.sqrt((np.array(P_err_aE[i])**2 + (np.array(P0_day_err)**2))) * np.array(E_j[i]))
#        print (OC_err)
        OC_s_err = OC_err*24*60*60
    else:
        P_E_day = Delta_T[i] / E_k[i]
        P_aE[i] = P_E_day
        P_E_err_day = np.abs((np.array(T_obs_err[i]) - np.array(T0_bjd_err)) / E_k[i])
        P_err_aE[i] = P_E_err_day
        T_C_linear = T0_bjd + P0_day*E_j[i]              #Linear
        T_aC_linear[i] = T_C_linear #arrays
#    print (T_O, T_C)
        OC = np.array(T_O_linear) - np.array(T_C_linear)
        OC_s = (np.array(T_O_linear) - np.array(T_C_linear))*24*60*60
        OC_err = np.abs(np.sqrt((np.array(P_err_aE[i])**2)) *np.array(E_k[i]))
#        print (OC_err)
        OC_s_err = OC_err*24*60*60
    print ('%0.0f\t%0.6f\t%0.0f\t%0.6f\t%0.6f\t%0.6f\t%0.6f' %(i, BJD_time[i], Cycle[i], T_O_linear, T_C_linear, OC_s, OC_s_err))
    OC_cal.append('%0.6f\t%0.0f\t%0.6f\t%0.6f\t%0.6f\t%0.6f' %(BJD_time[i], Cycle[i], T_O_linear, T_C_linear, OC_s, OC_s_err))

#P_aver = mean(P_aE[i])
#P_aver_a[i] = P_aver
#P_aver_std = np.std(P_aE[i])
#P_aver_std = mean(P_err_aE[i])
#P_aver_std_a[i] = P_aver_std
#print('%0.11f %0.11f' %(P_aver, P_aver_std))

rerults = OC_cal
f = open("Kittipong_2020.out", 'w')
for i in range(len(rerults)):
    f.write(str(rerults[i])+ '\n')
f.close()

#Plot O-C vs BJD
Input_Kittipong_2020  = 'Kittipong_2020.out'
Data_Kittipong_2020   = np.genfromtxt(Input_Kittipong_2020)

N = 2450000
BJD_time_Kittipong_2020 = Data_Kittipong_2020[:,0] - N
Cycle_Kittipong_2020 = Data_Kittipong_2020[:,1]
T_O_linear_Kittipong_2020 = Data_Kittipong_2020[:,2]
T_C_linear_Kittipong_2020 = Data_Kittipong_2020[:,3]
OC_s_Kittipong_2020 = Data_Kittipong_2020[:,4]
OC_s_err_Kittipong_2020 = Data_Kittipong_2020[:,5]

##Plotgraph
fig=plt.figure(figsize=(15, 5), tight_layout=True)
plt.tick_params(direction='in', which='both', bottom='on',top='on', right = 'on')

x1 = min(BJD_time_Kittipong_2020)
x2 = max(BJD_time_Kittipong_2020)
plt.errorbar(BJD_time_Kittipong_2020, OC_s_Kittipong_2020, yerr= OC_s_err_Kittipong_2020, fmt='o', markersize=8, color='red',
                    ecolor='lightgray', label = 'Kittipong et al. (2020)')

#Kittipong_2020
plt.text(BJD_time_Kittipong_2020[0], 56, '2014')
plt.text(BJD_time_Kittipong_2020[3], 56, '2015')
plt.text(BJD_time_Kittipong_2020[5]-50, 56, '2016')
plt.text(BJD_time_Kittipong_2020[7], 56, '2017')
plt.text(BJD_time_Kittipong_2020[10], 56, '2018')
plt.text(BJD_time_Kittipong_2020[12], 56, '2019')
plt.text(BJD_time_Kittipong_2020[14], 56, '2020')



#plt.hlines(y= 0, xmin= x1, xmax= x2, colors='k', linestyles='dotted')
#plt.xlim(x1,x2)
#plt.ylim(-20,20)
plt.xlabel('BJD - 2450000')
plt.ylabel('O-C (sec)')
plt.legend()
plt.grid(linestyle='dotted')
######plt.title('O-C diagram: DP Leo')
#output_filename = os.path.splitext(__file__)[0] + '.png'
#plt.savefig(output_filename, dpi=1000)
plt.savefig("OC_Kittipong_2020_rev1.png", dpi=1000)
plt.show()

##Plotgraph
fig=plt.figure(figsize=(15, 5), tight_layout=True)
plt.tick_params(direction='in', which='both', bottom='on',top='on', right = 'on')

#x1 = min(BJD_time_Schwope_2002)
#x2 = max(BJD_time_Schwope_2002)
#plt.errorbar(BJD_time, OC_s, yerr=OC_s_err, fmt='o', color='limegreen')
plt.errorbar(BJD_time_Schwope_2002, OC_s_Schwope_2002, yerr= OC_s_err_Schwope_2002, fmt='o', markersize=8, color='green',
                    ecolor='lightgray', label = 'Schwope et al. (2002)')

#Schwope_2002
plt.text(BJD_time_Schwope_2002[0], 155, '1979')
plt.text(BJD_time_Schwope_2002[3]-100, 155, '1981')
plt.text(BJD_time_Schwope_2002[13]-200, 155, '1984')
plt.text(BJD_time_Schwope_2002[17], 155, '1985')
plt.text(BJD_time_Schwope_2002[19]-200, 155, '1991')
plt.text(BJD_time_Schwope_2002[20]+100, 155, '1992')
plt.text(BJD_time_Schwope_2002[22]+200, 155, '1993')
plt.text(BJD_time_Schwope_2002[30]-50, 155, '2000')
plt.text(BJD_time_Schwope_2002[32]+50, 155, '2002')


plt.errorbar(BJD_time_Beuermann_2011, OC_s_Beuermann_2011, yerr= OC_s_err_Beuermann_2011, fmt='o', markersize=8, color='blue',
                    ecolor='lightgray', label = 'Beuermann et al. (2011)')

#Beuermann_2011
#plt.text(BJD_time_Beuermann_2011[0], 300, '2002')
plt.text(BJD_time_Beuermann_2011[1]-100, 155, '2009')
plt.text(BJD_time_Beuermann_2011[42]+100, 155, '2010')


plt.errorbar(BJD_time_Kittipong_2020, OC_s_Kittipong_2020, yerr= OC_s_err_Kittipong_2020, fmt='o', markersize=8, color='red',
                    ecolor='lightgray', label = 'Kittipong et al. (2020)')

#Kittipong_2020
#plt.text(BJD_time_Kittipong_2020[0]-200, 155, '2014')
#plt.text(BJD_time_Kittipong_2020[1]+50, 155, '2015')
#plt.text(BJD_time_Kittipong_2020[2]+100, 155, '2020')
#plt.text(BJD_time_Kittipong_2020[11]-300, 155, '2018')
#plt.text(BJD_time_Kittipong_2020[12]-100, 155, '2019')

plt.text(BJD_time_Kittipong_2020[0]-200, 155, '2014')
plt.text(BJD_time_Kittipong_2020[3], 155, '2015')
plt.text(BJD_time_Kittipong_2020[5]-100, 155, '2016')
#plt.text(BJD_time_Kittipong_2020[7], 155, '2017')
plt.text(BJD_time_Kittipong_2020[10], 155, '2018')
#plt.text(BJD_time_Kittipong_2020[12], 155, '2019')
plt.text(BJD_time_Kittipong_2020[14], 155, '2020')



#plt.hlines(y= 0, xmin= x1, xmax= x2, colors='k', linestyles='dotted')
#plt.xlim(x1,x2)
plt.ylim(-150,150)
plt.xlabel('BJD - 2450000')
plt.ylabel('O-C (sec)')
plt.legend()
plt.grid(linestyle='dotted')
######plt.title('O-C diagram: DP Leo')
#output_filename = os.path.splitext(__file__)[0] + '.png'
#plt.savefig(output_filename, dpi=1000)
plt.savefig("Schwope_Beuermann_Kittipong_rev1.png", dpi=1000)
plt.show()

##Plotgraph
fig=plt.figure(figsize=(15, 5), tight_layout=True)
plt.tick_params(direction='in', which='both', bottom='on',top='on', right = 'on')

#x1 = min(BJD_time_Schwope_2002)
#x2 = max(BJD_time_Schwope_2002)
#plt.errorbar(BJD_time, OC_s, yerr=OC_s_err, fmt='o', color='limegreen')

plt.errorbar(BJD_time_Schwope_2002[0:2], OC_s_Schwope_2002[0:2], yerr= OC_s_err_Schwope_2002[0:2], fmt='o', markersize=8, color='magenta',
                    ecolor='lightgray', label = 'X-ray')
plt.errorbar(BJD_time_Schwope_2002[15:16], OC_s_Schwope_2002[15:16], yerr= OC_s_err_Schwope_2002[15:16], fmt='o', markersize=8, color='magenta',
                    ecolor='lightgray')
plt.errorbar(BJD_time_Schwope_2002[20:31], OC_s_Schwope_2002[20:31], yerr= OC_s_err_Schwope_2002[20:31], fmt='o', markersize=8, color='magenta',
                    ecolor='lightgray')

plt.errorbar(BJD_time_Schwope_2002[19], OC_s_Schwope_2002[19], yerr= OC_s_err_Schwope_2002[19], fmt='o', markersize=8, color='darkorange',
                    ecolor='lightgray', label = 'UV')

plt.errorbar(BJD_time_Schwope_2002[3:14], OC_s_Schwope_2002[3:14], yerr= OC_s_err_Schwope_2002[3:14], fmt='o', markersize=8, color='teal',
                    ecolor='lightgray', label = 'Optical')
plt.errorbar(BJD_time_Schwope_2002[17:18], OC_s_Schwope_2002[17:18], yerr= OC_s_err_Schwope_2002[17:18], fmt='o', markersize=8, color='teal',
                    ecolor='lightgray')
plt.errorbar(BJD_time_Schwope_2002[32], OC_s_Schwope_2002[32], yerr= OC_s_err_Schwope_2002[32], fmt='o', markersize=8, color='teal',
                    ecolor='lightgray')

#Schwope_2002
plt.text(BJD_time_Schwope_2002[0], 155, '1979')
plt.text(BJD_time_Schwope_2002[3]-100, 155, '1981')
plt.text(BJD_time_Schwope_2002[13]-200, 155, '1984')
plt.text(BJD_time_Schwope_2002[17], 155, '1985')
plt.text(BJD_time_Schwope_2002[19]-200, 155, '1991')
plt.text(BJD_time_Schwope_2002[20]+100, 155, '1992')
plt.text(BJD_time_Schwope_2002[22]+200, 155, '1993')
plt.text(BJD_time_Schwope_2002[30]-50, 155, '2000')
plt.text(BJD_time_Schwope_2002[32]+50, 155, '2002')

plt.errorbar(BJD_time_Beuermann_2011, OC_s_Beuermann_2011, yerr= OC_s_err_Beuermann_2011, fmt='o', markersize=8, color='teal',
                    ecolor='lightgray')

#Beuermann_2011
plt.text(BJD_time_Beuermann_2011[1]-100, 155, '2009')
plt.text(BJD_time_Beuermann_2011[42]+100, 155, '2010')


plt.errorbar(BJD_time_Kittipong_2020, OC_s_Kittipong_2020, yerr= OC_s_err_Kittipong_2020, fmt='o', markersize=8, color='teal',
                    ecolor='lightgray')

#Kittipong_2020
#plt.text(BJD_time_Kittipong_2020[0]-200, 155, '2014')
#plt.text(BJD_time_Kittipong_2020[1]+50, 155, '2015')
#plt.text(BJD_time_Kittipong_2020[2]+100, 155, '2020')
#plt.text(BJD_time_Kittipong_2020[11]-300, 155, '2018')
#plt.text(BJD_time_Kittipong_2020[12]-100, 155, '2019')

plt.text(BJD_time_Kittipong_2020[0]-200, 155, '2014')
plt.text(BJD_time_Kittipong_2020[3], 155, '2015')
plt.text(BJD_time_Kittipong_2020[5]-100, 155, '2016')
#plt.text(BJD_time_Kittipong_2020[7], 155, '2017')
plt.text(BJD_time_Kittipong_2020[10], 155, '2018')
#plt.text(BJD_time_Kittipong_2020[12], 155, '2019')
plt.text(BJD_time_Kittipong_2020[14], 155, '2020')

#plt.hlines(y= 0, xmin= x1, xmax= x2, colors='k', linestyles='dotted')
#plt.xlim(x1,x2)
plt.ylim(-150,150)
plt.xlabel('BJD - 2450000')
plt.ylabel('O-C (sec)')
plt.legend()
plt.grid(linestyle='dotted')
######plt.title('O-C diagram: DP Leo')
#output_filename = os.path.splitext(__file__)[0] + '.png'
#plt.savefig(output_filename, dpi=1000)
plt.savefig("Schwope_Beuermann_Kittipong_obs_rev1.png", dpi=1000)
plt.show()
