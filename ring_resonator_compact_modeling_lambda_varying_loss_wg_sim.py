import time
start_time = time.time()
import pandas as pd
import sys, os
#sys.path.append(os.path.dirname(__file__)) #Current directory
#sys.path.append("/opt/lumerical/v222/api/python/") #Default windows lumapi path
import simfuncs as sf
#import lumapi
import numpy as np
import matplotlib.pyplot as plt
import subprocess


################## USER INPUT STARTS HERE #######################

wg_name='tlx_cl_te_mrr_airclad_fn_internal'

rrjson_src  = r'C:/Users/mh220218/Documents/My Assignments/Compact Model/Ring_resonator/NRL Data/rr_model_cml/input_files/rr_model_v0p1/source/' + str(wg_name)+ '/'+ 'input.json'
rrjson_dst  =r'C:/Users/mh220218/Documents/My Assignments/Compact Model/Ring_resonator/NRL Data/rr_model_cml/input_files/rr_model_v0p1/source/'+ str(wg_name)+'/' + str(wg_name) + '.json'

rrjsonparams = {'width_max': 4.9e-07}
cmlfiledir    = r'C:\Users\mh220218\Documents\My Assignments\Compact Model\Ring_resonator\NRL Data\rr_model_cml\input_files\rr_model_v0p1'
runcml  = False
simintc = True

test_df=pd.read_csv(r"C:\Users\mh220218\Documents\My Assignments\Compact Model\Ring_resonator\NRL Data\TLX_PDK1.0_MRR_C_L_band\tlx_cl_te_mrr_airclad_fn\tlx_cl_te_mrr_airclad_fn.txt", delimiter='\t')
test_analysis_df=pd.read_csv(r"C:\Users\mh220218\Documents\My Assignments\Compact Model\Ring_resonator\NRL Data\TLX_PDK1.0_MRR_C_L_band\tlx_cl_te_mrr_airclad_fn\tlx_cl_te_mrr_airclad_fn_Analysis.txt", delimiter='\t')
    
################## USER INPUT ENDS HERE #######################

# #### MODIFY JSON FILE HERE
dict_in    = sf.read_rr_json(rrjson_src)

print("\033[H\033[J") #To clear screen
print('THIS IS THE INPUT DICTIONARY')
#print(rrdict)

# METHOD 1
dict_out = dict_in
#rrdict_op['neff']['_data'] = 

# METHOD 2: NEEDS rrjsonparams
#rrdict_op = sf.modify_json(rrdict,rrjsonparams)
#print('THIS IS THE OUTPUT - MODIFIED DICTIONARY')
#print(rrdict_op)

#wg_neff_arr=GET FROM EQUATION
#dict_out['neff']['_data'][:]=wg_neff_arr

lam=dict_in['wavelength_data']['_data']
lam=np.array(lam)
f=(3e8)/lam
lam0=1.55e-6

# n0=1.728875
# n1=1.705
# n2=-4.5e5

n0=1.45033
n1=-1.65e5 #1.61
n2=2.8e11

neff_model=n0 + n1*(lam-lam0) + n2*(lam-lam0)**2
ng_model = neff_model-lam*(np.gradient(neff_model,lam))
D_model = (-2*np.pi/lam**2)*(np.gradient(ng_model,2*np.pi*f))

l0=5.6e2     #in dB/m
l1=4.5e9     
l2=-2e16     

loss_model=l0 + l1*(lam-lam0) + l2*(lam-lam0)**2

dict_out['neff']['_data'][:]=neff_model
dict_out['ng']['_data'][:]=ng_model
dict_out['D']['_data'][:]=D_model
dict_out['loss']['_data'][:]=loss_model   #dB/m

sf.write_rr_json(rrjson_dst, dict_out)

#-----------------------------------------------------------------------------------

runcml=1
if runcml:
    p = subprocess.run([r'C:\Program Files\Lumerical\v231\bin\cml-compiler.bat','all'],shell=True,cwd = cmlfiledir,check=True, stderr=subprocess.STDOUT)
    print('..........CML implementation DONE..........')

plot_wg=1
if plot_wg:
    fig, ([ax0, ax1], [ax2, ax3]) = plt.subplots(nrows=2, ncols=2, sharex=False, figsize=(20, 12))
    #plt.rcParams.update({'font.size': 10})
    
    model_json_file     = r'C:/Users/mh220218/Documents/My Assignments/Compact Model/Ring_resonator/NRL Data/rr_model_cml/input_files/rr_model_v0p1/source/'+ str(wg_name)+'/' + str(wg_name) + '.json'
    model_data = sf.read_wg_json(model_json_file)
                         
    # ax0.plot(test_df['Wavelength (nm)'], test_df['Normalized Transmission (dB)'], '-', label='Measured data', linewidth=2, color='darkblue')
    # ax0.set_xlabel('Wavelength (nm)')
    # ax0.set_ylabel('Transmission (dB)')
    # ax0.legend(frameon=False)    
    
    ax0.plot(test_analysis_df['Wavelength (nm)'], test_analysis_df['Loss (dB/cm)'], '.', label='Loss (measured)', linewidth=2, color='darkblue')
    ax0.plot(model_data['wavelength_data']*1e9, loss_model/100, label='Loss (model)', linewidth=2, color='red')
    ax0.set_xlabel('Wavelength (nm)')
    ax0.set_ylabel('Waveguide Loss (dB/cm)')  #temporal spread (ps) per unit propagation distance (km), per unit pulse spectral width (nm)
    ax0.ticklabel_format(useOffset=False)
    ax0.legend(frameon=False)
    ax0.tick_params(axis='both', which='major', labelsize=14)
    
    ax1.plot(model_data['wavelength_data']*1e9, model_data['neff'], '-', label='$n_{eff} (model)$', linewidth=2, color='red')
    ax1.set_xlabel('Wavelength (nm)')
    ax1.set_ylabel('Effective index')
    ax1.ticklabel_format(useOffset=False)
    ax1.legend(frameon=False)
    ax1.tick_params(axis='both', which='major', labelsize=14)
    
    ax2.plot(test_analysis_df['Wavelength (nm)'], test_analysis_df['ng'], '.', label='$n_g$ (measured)', linewidth=2, color='darkblue')
    ax2.plot(model_data['wavelength_data']*1e9, model_data['ng'], label='$n_g$ (model)', linewidth=2, color='red')
    ax2.set_xlabel('Wavelength (nm)')
    ax2.set_ylabel('Group index')
    ax2.ticklabel_format(useOffset=False)
    ax2.legend(frameon=False)
    ax2.tick_params(axis='both', which='major', labelsize=14)
    
    ax3.plot(model_data['wavelength_data']*1e9, model_data['D'], label='$D$ (model)', linewidth=2, color='red')
    ax3.set_xlabel('Wavelength (nm)')
    ax3.set_ylabel('Dispersion (s/m/m)')  #temporal spread (ps) per unit propagation distance (km), per unit pulse spectral width (nm)
    ax3.set_ylim([-0.0031, -0.0026])
    ax3.set_xlim([1510, 1632])
    ax3.ticklabel_format(useOffset=False)
    ax3.legend(frameon=False)
    ax3.tick_params(axis='both', which='major', labelsize=14)
    
    plt.rc('axes', labelsize=20)
    plt.rc('legend', fontsize=18)
    #plt.rc('axes', labelsize=12)
    print('..........Plotting DONE..........')