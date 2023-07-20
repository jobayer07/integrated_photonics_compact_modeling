import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy import signal
import os
import glob
from scipy import interpolate
from scipy.optimize import curve_fit
import re

name='rfs_qc_cl_dc50p_se'
meas_filename='C:/Users\mh220218/Documents/My Assignments/Compact Model/Quantum/directional_coupler/Summary/dc_cl_si.csv'
plot_measurement=1
target_coupling_length=14.2  #in um
die_id=['2_6', '3_5', '5_7', '6_2', '6_6', '7_1', '8_4', '9_3'] # ['1_3', '2_2', '5_7'']#'alldie'

port1_fdtd_file='L13p9_port1.txt'
port2_fdtd_file='L13p9_port2.txt'

colors = plt.cm.hsv(np.linspace(0, 1, 10))

#For S-parameter
portnames = ['opt_1','opt_2', 'opt_4', 'opt_3']
portlocs  = ['LEFT', 'LEFT', 'RIGHT', 'RIGHT']
polarization='TE'
mode_id=1


#--------------------------------------------Functions---------------------------------------------
def loss_modeling(die_id, filename, coupling_length, port1_fdtd_file, port2_fdtd_file):
    df = pd.read_csv(filename, header=0)
    df=df[df['source'].str.contains(r'V:\\AIM_EPDA_Test\\EPDA Group\\Master_Data_Location\\ICT_Data\\Quantum\\JeepJK\\Directonal_Couplers_Si_CL/L1_JeepJK_2207AMPM001_000_')]
    
    df_L=df[df['L']==coupling_length]
    df_L=df_L[['x_y_die', 'lambda', 'loss1f', 'loss2f']]
    df_L=df_L.dropna()
    
    if (die_id == 'alldie'):
        die_id=df.groupby('x_y_die').apply(list).keys().values
        print('Dies:', die_id)
        
    for d in range(len(die_id)):
        try:
            df_d=df_L[df_L['x_y_die']==die_id[d]]
            df_port1=df_d[['lambda', 'loss1f']]
    
            print(die_id[d])
            
            df_port2=df_d[['lambda', 'loss2f']]
            
            plt.figure(1)
            if plot_measurement==1:
                plt.plot(df_port1[df_port1.columns[0]], df_port1[df_port1.columns[1]], '--', color=colors[d], label=die_id[d], linewidth=0.5)
                plt.plot(df_port2[df_port2.columns[0]], df_port2[df_port2.columns[1]], '--', color=colors[d], linewidth=0.5)
                #plt.plot(df_port2[df_port2.columns[0]], abs(df_port1[df_port1.columns[1]]-df_port2[df_port2.columns[1]]), '--k', linewidth=0.5)
                #ax0.set_ylim([11.7, 17.3])
        except:
            print(die_id[d] + ' is absent')
    
    port1_fdtd = pd.read_csv(port1_fdtd_file, sep=", ", engine='python')
    port2_fdtd = pd.read_csv(port2_fdtd_file, sep=", ", engine='python')
            
    lam=port1_fdtd['lambda(m)'].values
    port1=port1_fdtd['Y'].values
    port2=port2_fdtd['Y'].values

    # plt.plot(lam*1e9, port1, label='Port1_fdtd', linewidth=2)
    # plt.plot(lam*1e9, port2, label='Port2_fdtd', linewidth=2)
    # plt.legend()
 
    loss_model1=-10*np.log10(port1)+2
    loss_model2=-10*np.log10(port2)+2
    
    plt.plot(lam*1e9, loss_model1, label='Port1_fdtd', linewidth=2)
    plt.plot(lam*1e9, loss_model2, label='Port2_fdtd', linewidth=2)
    plt.xlabel('Wavelength (nm)', fontsize=14)
    plt.ylabel('Loss (dB/cm)', fontsize=14)
    plt.title('FDTD simulation @ L=' + str (str(port1_fdtd_file).split('_')[0].replace('p', '.').replace('L', '')) + ' um')
    plt.legend(frameon=False, handletextpad=0.1, ncol=5, columnspacing=0.7)
    return loss_model1, loss_model2, lam


def convert_to_s_parameter(loss_model1, loss_model2, lam):
    np.set_printoptions(suppress=True, formatter={'float_kind':'{:f}'.format})
    f=3e8/lam
    T4=10**(-loss_model1/20) #dB=20log(x), here 20 is for power
    T4_phase=T4-T4
    T3=10**(-loss_model2/20) #dB=20log(x), here 20 is for power
    T3_phase=T3-T3
    
    
    #Construct S parameters: S12=S21, S11=S22
    N=len(loss_model1)
    
    S41=np.zeros([N, 3], dtype=float)
    S41[:, 0]=f
    S41[:, 1]=T4
    S41[:, 2]=T4_phase
    
    S31=np.zeros([N, 3], dtype=float)
    S31[:, 0]=f
    S31[:, 1]=T3
    S31[:, 2]=T3_phase
    
    #Zeros
    S21=np.zeros([N, 3], dtype=float)
    S21[:, 0]=f
    S21[:, 1]=T4-T4
    S21[:, 2]=T4_phase-T4_phase
    
    S43=np.zeros([N, 3], dtype=float)
    S43[:, 0]=f
    S43[:, 1]=T4-T4
    S43[:, 2]=T4_phase-T4_phase
    
    #Symmetries
    S13=S31
    S23=S41
    S14=S41
    S24=S31
    S42=S31
    S32=S41

    #Making retun losses=0
    S11=np.zeros([N, 3], dtype=float)
    S11[:, 0]=f
    S11[:, 1]=T4-T4
    S11[:, 2]=T4_phase
    
    S22=S11
    S33=S11
    S44=S11
    
    #Making Cross coupling to zero
    S21=np.zeros([N, 3], dtype=float)
    S21[:, 0]=f
    S21[:, 1]=T4-T4
    S21[:, 2]=T4_phase
    
    S12=S21
    S34=S21
    S43=S21    
    
    #file_object  = open("filename", "mode")
    txt_file= open(str(name) +".txt","w+")
    
    #Locations
    txt_file.write('[\"' + portnames[0] + '\",\"' + portlocs[0] + '\"]\n')
    txt_file.write('[\"' + portnames[1] + '\",\"' + portlocs[1] + '\"]\n')
    txt_file.write('[\"' + portnames[2] + '\",\"' + portlocs[2] + '\"]\n')
    txt_file.write('[\"' + portnames[3] + '\",\"' + portlocs[3] + '\"]\n')
    
    #S11
    txt_file.write('(\"' + portnames[0] + '\",\"' +str(polarization)+ '\",' +str(mode_id) + ',\"' + portnames[0] +'\",' + str(mode_id) +',\"transmission\")\n')
    txt_file.write('('+ str(N) +',3)\n')
    txt_file.write(re.sub('( \[|\[|\])', '', str(S11)))
    
    #S21
    txt_file.write('\n(\"' + portnames[1] + '\",\"' +str(polarization)+ '\",' +str(mode_id) + ',\"' + portnames[0] +'\",' + str(mode_id) +',\"transmission\")\n')
    txt_file.write('('+ str(N) +',3)\n')
    txt_file.write(re.sub('( \[|\[|\])', '', str(S21)))
   
    #S41
    txt_file.write('\n(\"' + portnames[2] + '\",\"' +str(polarization)+ '\",' +str(mode_id) + ',\"' + portnames[0] +'\",' + str(mode_id) +',\"transmission\")\n')
    txt_file.write('('+ str(N) +',3)\n')
    txt_file.write(re.sub('( \[|\[|\])', '', str(S41))) 
   
    #S31
    txt_file.write('\n(\"' + portnames[3] + '\",\"' +str(polarization)+ '\",' +str(mode_id) + ',\"' + portnames[0] +'\",' + str(mode_id) +',\"transmission\")\n')
    txt_file.write('('+ str(N) +',3)\n')
    txt_file.write(re.sub('( \[|\[|\])', '', str(S31))) 
    
    #S33
    txt_file.write('\n(\"' + portnames[3] + '\",\"' +str(polarization)+ '\",' +str(mode_id) + ',\"' + portnames[3] +'\",' + str(mode_id) +',\"transmission\")\n')
    txt_file.write('('+ str(N) +',3)\n')
    txt_file.write(re.sub('( \[|\[|\])', '', str(S33))) 

    #S13
    txt_file.write('\n(\"' + portnames[0] + '\",\"' +str(polarization)+ '\",' +str(mode_id) + ',\"' + portnames[3] +'\",' + str(mode_id) +',\"transmission\")\n')
    txt_file.write('('+ str(N) +',3)\n')
    txt_file.write(re.sub('( \[|\[|\])', '', str(S13))) 

    #S23
    txt_file.write('\n(\"' + portnames[1] + '\",\"' +str(polarization)+ '\",' +str(mode_id) + ',\"' + portnames[3] +'\",' + str(mode_id) +',\"transmission\")\n')
    txt_file.write('('+ str(N) +',3)\n')
    txt_file.write(re.sub('( \[|\[|\])', '', str(S23))) 
    
    #S43
    txt_file.write('\n(\"' + portnames[2] + '\",\"' +str(polarization)+ '\",' +str(mode_id) + ',\"' + portnames[3] +'\",' + str(mode_id) +',\"transmission\")\n')
    txt_file.write('('+ str(N) +',3)\n')
    txt_file.write(re.sub('( \[|\[|\])', '', str(S43))) 
    
    #S44
    txt_file.write('\n(\"' + portnames[2] + '\",\"' +str(polarization)+ '\",' +str(mode_id) + ',\"' + portnames[2] +'\",' + str(mode_id) +',\"transmission\")\n')
    txt_file.write('('+ str(N) +',3)\n')
    txt_file.write(re.sub('( \[|\[|\])', '', str(S44))) 
    
    #S34
    txt_file.write('\n(\"' + portnames[3] + '\",\"' +str(polarization)+ '\",' +str(mode_id) + ',\"' + portnames[2] +'\",' + str(mode_id) +',\"transmission\")\n')
    txt_file.write('('+ str(N) +',3)\n')
    txt_file.write(re.sub('( \[|\[|\])', '', str(S34))) 
    
    #S14
    txt_file.write('\n(\"' + portnames[0] + '\",\"' +str(polarization)+ '\",' +str(mode_id) + ',\"' + portnames[2] +'\",' + str(mode_id) +',\"transmission\")\n')
    txt_file.write('('+ str(N) +',3)\n')
    txt_file.write(re.sub('( \[|\[|\])', '', str(S14))) 
    
    #S24
    txt_file.write('\n(\"' + portnames[1] + '\",\"' +str(polarization)+ '\",' +str(mode_id) + ',\"' + portnames[2] +'\",' + str(mode_id) +',\"transmission\")\n')
    txt_file.write('('+ str(N) +',3)\n')
    txt_file.write(re.sub('( \[|\[|\])', '', str(S24))) 
    
    #S22
    txt_file.write('\n(\"' + portnames[1] + '\",\"' +str(polarization)+ '\",' +str(mode_id) + ',\"' + portnames[1] +'\",' + str(mode_id) +',\"transmission\")\n')
    txt_file.write('('+ str(N) +',3)\n')
    txt_file.write(re.sub('( \[|\[|\])', '', str(S22)))
    
    #S42
    txt_file.write('\n(\"' + portnames[2] + '\",\"' +str(polarization)+ '\",' +str(mode_id) + ',\"' + portnames[1] +'\",' + str(mode_id) +',\"transmission\")\n')
    txt_file.write('('+ str(N) +',3)\n')
    txt_file.write(re.sub('( \[|\[|\])', '', str(S42)))
    
    #S32
    txt_file.write('\n(\"' + portnames[3] + '\",\"' +str(polarization)+ '\",' +str(mode_id) + ',\"' + portnames[1] +'\",' + str(mode_id) +',\"transmission\")\n')
    txt_file.write('('+ str(N) +',3)\n')
    txt_file.write(re.sub('( \[|\[|\])', '', str(S32)))
    
    #S12
    txt_file.write('\n(\"' + portnames[0] + '\",\"' +str(polarization)+ '\",' +str(mode_id) + ',\"' + portnames[1] +'\",' + str(mode_id) +',\"transmission\")\n')
    txt_file.write('('+ str(N) +',3)\n')
    txt_file.write(re.sub('( \[|\[|\])', '', str(S12)))
    
    txt_file.close()
    print('Successfully converted to S parameters')
    return


#---------------------------------------------------------Main program-----------------------------------------

port4_loss, port3_loss, lam=loss_modeling(die_id, meas_filename, target_coupling_length, port1_fdtd_file, port2_fdtd_file)
#convert_to_s_parameter(port4_loss, port3_loss, lam)



