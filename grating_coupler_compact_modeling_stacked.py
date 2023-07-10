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


name='rfs_qc_cl_gc_fn'
die_id='alldie' # ['1_3', '2_2', '5_7']#'alldie'
model_die='2_2'
target_wl=[1550] 
lam=np.arange(1485, 1600, 2)*1e-9
fig, ([ax0, ax1]) = plt.subplots(nrows=1, ncols=2, sharex=False, figsize=(20, 8))

#For S-parameter
portnames = ['opt_1_fbr','opt_2_fn']
portlocs  = ['LEFT', 'RIGHT']
polarization='TE'
mode_id=1

#----------------------------------------------------------------Functions--------------------------------------------
def wavelength_interpolation(x,y, xnew): 
    tck = interpolate.splrep(x, y, s=1)     #s=0 -> No smoothing/regression required
    ynew = interpolate.splev(xnew, tck, der=0)
    return ynew

def loss_modeling(filename, die_id):
    df = pd.read_csv(filename, header=0)
    
    if (die_id == 'alldie'):
        die_id=df.groupby('x_y_die').apply(list).keys()
        #print(die_id)
    
    my_dict=dict.fromkeys(set(target_wl), [])
        
    for d in range(len(die_id)):
        try:
            df_d=df[df['x_y_die']==die_id[d]]
            df_d=df_d[['lambda', 'loss1']]
            df_d=df_d.dropna()
            
            #print(die_id[d])
            
            ax0.plot(df[df.columns[0]], df[df.columns[1]], '--', label=die_id[d], linewidth=1)
            ax0.set_xlabel('Wavelength (nm)', fontsize=18)
            ax0.set_ylabel('Loss (dB/cm)', fontsize=18)
            ax0.set_ylim([11.7, 17.3])
            ax0.legend(frameon=False)
            
            for j in range (len(target_wl)):
                target_wl_loss=wavelength_interpolation(df_d[df.columns[0]].values,df_d[df.columns[1]].values, target_wl[j])
                my_dict[target_wl[j]].append(target_wl_loss)
            
        except:
            print(die_id[d] + ' is absent')
    
    df_1550 = pd.DataFrame(my_dict)
    ax1.boxplot(df_1550, 1)
    #ax1.set_title('Notched plot @ 1550 nm wavelength', fontsize=18)
    ax1.set_xticks([1], [target_wl[j]])
    ax1.tick_params(axis='both', which='major', labelsize=14)
    ax1.set_xlabel('Wavelength (nm)', fontsize=18)
    ax1.set_ylabel('Loss (dB/cm)', fontsize=18)

    model_df=df[df['x_y_die']==model_die]
    model_df=model_df[['lambda', 'loss1']]
    model_df=model_df.dropna()
    #print(model_df)
    model_coeff = np.polyfit(model_df[model_df.columns[0]].values*1e-9, model_df[model_df.columns[1]].values, 4)
    loss_model=np.polyval(model_coeff, lam)+0.15
    
    ax0.plot(lam*1e9, loss_model, '-', label='Model', linewidth=4, color='red')
    ax0.legend(frameon=False, fontsize=14)
    ax0.tick_params(axis='both', which='major', labelsize=14)
     
    return loss_model

def convert_to_s_parameter(loss_model):
    np.set_printoptions(suppress=True, formatter={'float_kind':'{:f}'.format})
    f=3e8/lam
    T=10**(-loss_model/20) #dB=20log(x), here 20 is for power
    T_phase=T-T
    
    #Construct S parameters: S12=S21, S11=S22
    N=len(loss_model)
    S21=np.zeros([N, 3], dtype=float)
    S21[:, 0]=f
    S21[:, 1]=T
    S21[:, 2]=T_phase
    
    S12=S21
    
    S11=np.zeros([N, 3], dtype=float)
    S11[:, 0]=f
    S11[:, 1]=T-T
    S11[:, 2]=T_phase
    
    S22=S11
    
    #file_object  = open("filename", "mode")
    txt_file= open(str(name) +".txt","w+")
    
    #S11
    txt_file.write('[\"' + portnames[0] + '\",\"' + portlocs[0] + '\"]\n')
    txt_file.write('[\"' + portnames[1] + '\",\"' + portlocs[1] + '\"]\n')
    txt_file.write('(\"' + portnames[0] + '\",\"' +str(polarization)+ '\",' +str(mode_id) + ',\"' + portnames[0] +'\",' + str(mode_id) +',\"transmission\")\n')
    txt_file.write('('+ str(N) +',3)\n')
    txt_file.write(re.sub('( \[|\[|\])', '', str(S11)))
    
    #S21
    txt_file.write('\n(\"' + portnames[1] + '\",\"' +str(polarization)+ '\",' +str(mode_id) + ',\"' + portnames[0] +'\",' + str(mode_id) +',\"transmission\")\n')
    txt_file.write('('+ str(N) +',3)\n')
    txt_file.write(re.sub('( \[|\[|\])', '', str(S21)))
    
    #S22
    txt_file.write('\n(\"' + portnames[1] + '\",\"' +str(polarization)+ '\",' +str(mode_id) + ',\"' + portnames[1] +'\",' + str(mode_id) +',\"transmission\")\n')
    txt_file.write('('+ str(N) +',3)\n')
    txt_file.write(re.sub('( \[|\[|\])', '', str(S22)))
    
    #S12
    txt_file.write('\n(\"' + portnames[0] + '\",\"' +str(polarization)+ '\",' +str(mode_id) + ',\"' + portnames[1] +'\",' + str(mode_id) +',\"transmission\")\n')
    txt_file.write('('+ str(N) +',3)\n')
    txt_file.write(re.sub('( \[|\[|\])', '', str(S12)))
    
    txt_file.close()
    print('Successfully converted to S parameters')
    return


#----------------------------------------------------------------Main-------------------------------------------------

filename='L8_JeepJK_2227AMPM002_000_7D5HE245SOG0_GC_PDA_pdgc.csv'
loss_model=loss_modeling(filename, die_id)
convert_to_s_parameter(loss_model)
