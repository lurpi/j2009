# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 09:52:34 2018

@author: Luca Urpi
"""


import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os as os

#%% define Mohr Coulomb parameters and joint paramenters

G=7.E6 # Shear modulus
K=1.e8 # Bulk modulus
csn=2.e3 #cohesion
frict_angle=40.0 #friction angle (degrees)
dilat_angle=0.0 #dilation angle (degrees)
tension=2.4e3 #tension limit

#G_joint=7.E6 # Shear modulus
#K_joint=1.e8 # Bulk modulus
csn_joint=1.e3 #cohesion
frict_angle_joint=30.0 #friction angle (degrees)
dilat_angle_joint=0.0 #dilation angle (degrees)
tension_joint=2e3 #tension limit

#%% find first occurence of rupture in output file
def find_ogs_failure_szz(path, name, chist):
    in_file=path+name
    p0 = pd.read_csv(in_file, skiprows=2, header=0, names=chist, delimiter=r"\s+")  
    time=p0["TIME"].values
    str_pls=p0["STRAIN_PLS"].values
    s_zz=p0["STRESS_ZZ"].values
    s_fail=[0]
    for i in range(len(time)):
        check = str_pls[i]
        if check==0:
            pass
        else:
            s_fail.append(s_zz[i])
    return s_fail[1]
#%% ----------------------------------------------

#%% find tested dip angles folder
def dip_angles(name):
    dip=[]
    explore=range(0,360)
    for i in explore:
        path=name+str(explore[i])
        path0=name+"0"+str(explore[i])
        if os.path.isdir(path) :
            dip.append(path[-2:])
        elif os.path.isdir(path0):
            dip.append(path0[-2:])
    return dip
#%% ----------------------------------------------


#%% initialize files and value readings
name_dir="dip_"
dips=dip_angles(name_dir)
weak_failure=[]
matrix_failure=[]
tot_failure=[]
ogs_failure=[]
#%% ----------------------------------------------
#%% iterate through solved directories 
for i_dip in dips:
    i=90-float(i_dip)
    # read value from OGS
    columns=("TIME","STRAIN_PLS","PRESSURE1","STRESS_XX","STRESS_YY","STRESS_ZZ","STRESS_XY","STRESS_XZ","STRESS_YZ","VELOCITY_X1","VELOCITY_Y1","VELOCITY_Z1","DISPLACEMENT_X1","DISPLACEMENT_Y1","DISPLACEMENT_Z1","p_(1st_Invariant)","q_(2nd_Invariant)","Effective_Strain")
    numf=find_ogs_failure_szz(name_dir+i_dip,"\\joint_tri_time_INJ_TOP.tec", columns)
    # calcualte analtical values
    #%% details of the procedure can be found in Jaeger (1960)
    k=1-np.tan(np.radians(frict_angle_joint))*np.tan(np.radians(i))
    if (k < 0):
        k=1e-16
    weakf=2*csn_joint/(k*np.sin(2*np.radians(i)))
    matrf=2*csn*np.sqrt((1+np.sin(np.radians(frict_angle)))/(1-np.sin(np.radians(frict_angle))))
    #%% find failure values and append them
    totf=min(weakf,matrf)
    weak_failure.append(weakf)
    matrix_failure.append(matrf)
    tot_failure.append(totf)
    ogs_failure.append(np.abs(numf))
    #
     

#%% plot
plt.figure(figsize=(12,5))
#plt.plot(dips, weak_failure)
plt.plot(dips, matrix_failure, 'b--', alpha=0.6, label="Matrix")
plt.plot(dips, tot_failure, 'r-', alpha=0.6, label="Analytical value")
plt.plot(dips, ogs_failure, 'mo', alpha=0.6, label="Numerical values")
#plt.title('Stress evolution in the fault, point %s' %r3, fontsize=24)
plt.xlabel('Dip angle', fontsize=20)
plt.ylabel('Compressive strength', fontsize=20)
plt.tick_params(axis='x', labelsize=20)
plt.tick_params(axis='y', labelsize=20)
#set_ylabel('treatment')
plt.legend()
plt.grid()
plt.show()

