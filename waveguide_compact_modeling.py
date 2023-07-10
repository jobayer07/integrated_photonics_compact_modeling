import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy import signal
import os
import simfuncs as sf


wg_name='rfs_qc_o_wg_fn'
lam=np.arange(1240, 1380, 1)*1e-9
fig, ([ax0, ax1], [ax2, ax3]) = plt.subplots(nrows=2, ncols=2, sharex=False, figsize=(20, 12))

#----------------------------------------------------------------Functions--------------------------------------------
def loss_modeling(directory):
    dir_list = os.listdir(directory)
    Name= {'Name': dir_list}
    files = pd.DataFrame(data=Name)
    for i in range(len(files)):
        #print(filename_df.Name[i])
        df = pd.read_csv(str(directory)+'/'+files.Name[i], header=0)
        ax0.plot(df['Unnamed: 0'], df['slope'], '--', linewidth=1)
        ax0.set_xlabel('Wavelength (nm)')
        ax0.set_ylabel('Loss (dB/cm)')
        ax0.set_ylim([0, 0.6])
        
    l0=4.15
    l1=-3e+06
    loss_model=l0+l1*lam

    # model_coeff = np.polyfit(df['Unnamed: 0']*1e-9, df['slope'], 1)
    # print(model_coeff)
    # loss_model=np.polyval(model_coeff, lam)

    ax0.plot(lam*1e9, loss_model, '-', label='Model', linewidth=2, color='red')
    ax0.legend(frameon=False)
    ax0.tick_params(axis='both', which='major', labelsize=14)
    return loss_model

def neff_modeling(filename):
    neff_mode = pd.read_csv(filename)
    neff_mode=neff_mode.fillna('0')
    neff_mode[' Effective Index'] = neff_mode[' Effective Index'].replace({' Effective Index':'0'})
    neff_mode['Wavelength (microns)'] = neff_mode['Wavelength (microns)'].replace({'Effective Index':'0', 'Wavelength (microns)':'0'})
    neff_mode['Wavelength (microns)']=pd.to_numeric(neff_mode['Wavelength (microns)'], errors='coerce')
    neff_mode[' Effective Index']=pd.to_numeric(neff_mode[' Effective Index'], errors='coerce')
    neff_mode = neff_mode[(neff_mode[['Wavelength (microns)',' Effective Index']] != 0).all(axis=1)]
    
    n0=1.5756
    n1=-2.4e5
    n2=1.48e11
    neff_model=n0 + n1*(lam-lam0) + n2*(lam-lam0)**2
    ng_model = neff_model-lam*(np.gradient(neff_model,lam))
    D_model = (-2*np.pi/lam**2)*(np.gradient(ng_model,2*np.pi*f))
    
    ax1.plot(lam*1e9, neff_model, '-', label='Model', linewidth=2, color='red')
    ax1.plot(neff_mode['Wavelength (microns)']*1e3, neff_mode[' Effective Index'], '.', label='MODE simulation', linewidth=2, color='darkblue')
    ax1.set_xlabel('Wavelength (nm)')
    ax1.set_ylabel('Effective index')
    ax1.ticklabel_format(useOffset=False)
    ax1.legend(frameon=False)
    ax1.tick_params(axis='both', which='major', labelsize=14)
    
    ax2.plot(lam*1e9, ng_model, label='Model', linewidth=2, color='red')
    ax2.set_xlabel('Wavelength (nm)')
    ax2.set_ylabel('Group index')
    ax2.ticklabel_format(useOffset=False)
    ax2.legend(frameon=False)
    ax2.tick_params(axis='both', which='major', labelsize=14)
    
    ax3.plot(lam*1e9, D_model, label='Model', linewidth=2, color='red')
    ax3.set_xlabel('Wavelength (nm)')
    ax3.set_ylabel('Dispersion (s/m/m)')  #temporal spread (ps) per unit propagation distance (km), per unit pulse spectral width (nm)
    #ax3.set_ylim([-0.0014, -0.0011])
    ax3.set_xlim([1242, 1377])
    ax3.ticklabel_format(useOffset=False)
    ax3.legend(frameon=False)
    ax3.tick_params(axis='both', which='major', labelsize=14)
    
    plt.rc('axes', labelsize=20)
    plt.rc('legend', fontsize=18)
    print('..........Plotting DONE..........')
    return neff_model, ng_model, D_model
#----------------------------------------------------------------Main-----------------------------------------------------------

f=(3e8)/lam
lam0=1.31e-6
loss_directory='wg_loss_measured'
neff_filename='rfs_qc_o_wg_fn_MODE.txt'

loss_model=loss_modeling(loss_directory)
neff_model, ng_model, D_model=neff_modeling(neff_filename)



modeling=1
if (modeling):
    rrjson_src  = r'C:/Users/mh220218/Documents/My Assignments/Compact Model/Quantum/00model_cml/input_files/model_v0p1/source/' + str(wg_name)+ '/'+ 'input.json'
    rrjson_dst  =r'C:/Users/mh220218/Documents/My Assignments/Compact Model/Quantum/00model_cml/input_files/model_v0p1/source/'+ str(wg_name)+'/' + str(wg_name) + '.json'
    
    rrjsonparams = {'width_max': 4.9e-07}
    cmlfiledir    = r'C:/Users/mh220218/Documents/My Assignments/Compact Model/Quantum/00model_cml/input_files/model_v0p1'
    runcml  = False
    simintc = True
    
    dict_in    = sf.read_rr_json(rrjson_src)
    print("\033[H\033[J") #To clear screen
    print('THIS IS THE INPUT DICTIONARY')
    #print(rrdict)
    dict_out = dict_in
    dict_out['neff']['_data'][:]=neff_model
    dict_out['neff']['_size'][0]=len(lam)
    dict_out['ng']['_data'][:]=ng_model
    dict_out['ng']['_size'][0]=len(lam)
    dict_out['D']['_data'][:]=D_model
    dict_out['D']['_size'][0]=len(lam)
    dict_out['loss']['_data'][:]=loss_model*100   #dB/m
    dict_out['loss']['_size'][0]=len(lam)
    dict_out['wavelength_data']['_data'][:]=lam
    dict_out['wavelength_data']['_size'][0]=len(lam)
    sf.write_rr_json(rrjson_dst, dict_out)
    print('Modeling done!')
    
