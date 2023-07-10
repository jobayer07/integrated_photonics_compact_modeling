import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy import interpolate

def pairwiseAvg(lst):
    size = len(lst)
    x=np.zeros(size-1)
    for i in range(len(lst)-1):
        x[i] = int((lst[i] + lst[i + 1])/2)
    return x

def find_index_of_nearest_value(array, value):
    idx = (np.abs(array - value)).argmin()
    return idx

def find_mid_index_er(d):
    nofitted_peak_position, nofitted_peak_value = find_peaks(d, height=10, distance=50)
    mid_index=pairwiseAvg(nofitted_peak_position)  
    
    pairwisePeak=pairwiseAvg(nofitted_peak_value.get('peak_heights'))
    mid_value=d[mid_index]
    
    er=pairwisePeak-mid_value                                   #extinction ratio
    return mid_index, er

def ng_calculation(d, mid_index, fsr, r, CL):
    fsr=fsr*1e-9
    lam=d[mid_index]*1e-9 -fsr/2 # By -fsr/2, we calculate at the resonance peaks
    L=2*np.pi*r + 2*CL
    ng=(lam**2)/(L*fsr)
    return ng

def wavelength_interpolation(x,y, xnew): 
    tck = interpolate.splrep(x, y, s=40)     #s=0 -> No smoothing/regression required
    ynew = interpolate.splev(xnew, tck, der=0)
    return ynew


def q_calculation(wl_all, loss, mid_index, peak_index):
    a=loss.values
    fwhm_line=loss[mid_index]+10*np.log10(2)           #Half power in dB =10log(2) = 3 dB. At mid position
    fwhm=np.zeros(len(fwhm_line))
    wl=np.zeros(len(fwhm_line))

    for i in range(len(fwhm_line)):
        y=a[peak_index[i]:int(mid_index[i])]
        x=wl_all[peak_index[i]:int(mid_index[i])]
        x=x.to_numpy()
        
        dwl_all=wl_all[1]-wl_all[0]
        xnew=np.arange(x[0], x[len(x)-1], dwl_all/100)
        ynew=wavelength_interpolation(x, y, xnew)
        
        start=wl_all[peak_index[i]]
        
        index=find_index_of_nearest_value(ynew, fwhm_line.values[i])       #Now found the index of the fwhm point
        end=xnew[index]
        fwhm[i]=2*(end-start)
    
        wl[i]=wl_all[peak_index[i]]

    q=wl/fwhm
    return q

def calculate_neff(r, CL, lam, m):
    L=2*np.pi*r+2*CL
    lam=lam.to_numpy(dtype=None, copy=False)
    lam_L=lam/L
    
    m_array=np.zeros(len(lam))
    for j in range(len(lam)):
        m_array[j]=m-j
    neff=m_array*lam_L
    
    return neff
#-----------------------------------------------   Main Code  -------------------------------------------------------------
ring_radius=3e-6
coupling_length=2e-6
df = pd.read_csv('3um_through_data_cleaning.csv', skiprows=0)

mid_pos, er=find_mid_index_er(df['loss'])

peak_position, peak_value = find_peaks(df['loss_fit'], height=10, distance=50)
delta_wl=(df['wl'][len(df['wl'])-1]-df['wl'][0])/(len(df['wl'])-1)#df['wl'][1]-df['wl'][0]                                    #Wavelength resolution

fsr=np.diff(peak_position)*delta_wl                                 #fsr calculation

ng=ng_calculation(df['wl'], mid_pos, fsr, ring_radius, coupling_length)

Q=q_calculation(df['wl'], df['loss_fit'], mid_pos, peak_position)

lam=df['wl'][peak_position]*1e-9
neff=calculate_neff(ring_radius, coupling_length, lam, 37)

#dneff=np.diff(neff)
#dlam=np.diff(lam)

#lam_used=lam[ : -1]
#neff_used=neff[ : -1]
#ng_calc=neff_used-lam_used*dneff/dlam
#print('neff_measured:', neff)
#print('ng_measured:', ng_calc)


#--------------------------------Compact modeling stage-------------------------------------------

def interconnect_data_cleaning(filename):
    df= pd.read_csv(str(filename), skiprows=0)
    df= df.rename({'wavelength(nm)': 'wl', ' Y': 'loss'}, axis='columns')
    df['loss']= -df['loss']
    
    df=df.sort_values(by=['wl'])
    df = df.reset_index(drop=True)
    return df

ic_file='3umRadius_2umCL.txt'

df_ic=interconnect_data_cleaning(ic_file)

mid_pos_ic, er_ic=find_mid_index_er(df_ic['loss'])

peak_position_ic, peak_value_ic = find_peaks(df_ic['loss'], height=5, distance=30)
delta_wl_ic=(df_ic['wl'][len(df_ic['wl'])-1]-df_ic['wl'][0])/(len(df_ic['wl'])-1)#df_ic['wl'][1]-df_ic['wl'][0]                                    #Wavelength resolution

fsr_ic=np.diff(peak_position_ic)*delta_wl_ic                                 #fsr calculation
ng_ic=ng_calculation(df_ic['wl'], mid_pos_ic, fsr_ic, ring_radius, coupling_length)
ng_ic=ng_ic.to_numpy(dtype=None, copy=False)
Q_ic=q_calculation(df_ic['wl'], df_ic['loss'], mid_pos_ic, peak_position_ic)

#print('ng_model:', ng_ic)

lam_ic=df_ic['wl'][peak_position_ic]*1e-9
neff_ic=calculate_neff(ring_radius, coupling_length, lam_ic, 38)



#-------------------------------------------------Saving-----------------------------------------------
param_sum = {'wl': df['wl'][mid_pos], 'fsr': fsr, 'ng': ng, 'Q': Q}
param_sum = pd.DataFrame(data=param_sum)
param_sum.to_csv('3um_through_parameter.csv', index=False)

param_sum_model = {'wl': df_ic['wl'][mid_pos], 'fsr': fsr, 'ng': ng, 'Q': Q}
param_sum_model = pd.DataFrame(data=param_sum_model)
param_sum_model.to_csv('3um_through_model_parameter.csv', index=False)

#------------------------------------------------Plotting----------------------------------------------

fig, ([ax0, ax1, ax2], [ax3, ax4, ax5]) = plt.subplots(nrows=2, ncols=3, sharex=False, figsize=(18, 12))

ax0.plot(df['wl'], df['loss'], '--', label='Baseline corrected (measured data)', linewidth=2, color='darkblue')
ax0.set_xlabel('Wavelength (nm)')
ax0.set_ylabel('Loss (dB)')
#ax0.plot(df['wl'][peak_position], peak_value['peak_heights'], 'x', linewidth=2, color='darkblue')

ax1.plot(df['wl'][mid_pos], fsr, 'o--', label='FSR (measured data)', linewidth=2, color='darkblue')
ax1.set_xlabel('Wavelength (nm)')
ax1.set_ylabel('fsr (nm)')
#ax1.set_ylim([10, 18])
ax1.legend(frameon=False)

ax2.plot(df['wl'][mid_pos], ng, 'o--', label='$n_g$ (measured data)', linewidth=2, color='darkblue')
ax2.set_xlabel('Wavelength (nm)')
ax2.set_ylabel('Refractive index')
#ax2.set_ylim([1, 6])
ax2.legend(frameon=False)

ax3.plot(lam*1e9, neff, 'o-.', label='$n_{eff}$ (measured data)', linewidth=2, color='darkblue')
ax3.set_xlabel('Wavelength (nm)')
ax3.set_ylabel('Refractive index')
#ax3.set_ylim([1, 6])
ax4.legend(frameon=False)

ax4.plot(df['wl'][mid_pos], er, 'o--', label='Extinction ratio (measured data)', linewidth=2, color='darkblue')
ax4.set_xlabel('Wavelength (nm)')
ax4.set_ylabel('extinction ratio (dB)')
#ax4.set_ylim([20, 50])
ax4.legend(frameon=False)

ax5.plot(df['wl'][mid_pos], Q, 'o--', label='Q factor (measured data)', linewidth=2, color='darkblue')
ax5.set_xlabel('Wavelength (nm)')
#ax5.set_ylim([250, 500])
ax5.set_ylabel('Q factor')


#-----------------------------------------Plotting Compact Model----------------------------------------------

#fig, ([ax0, ax1, ax2], [ax3, ax4, ax5]) = plt.subplots(nrows=2, ncols=3, sharex=False, figsize=(18, 12))

ax0.plot(df_ic['wl'], df_ic['loss'], label='Model', linewidth=2, color='red')
ax0.set_xlabel('Wavelength (nm)')
ax0.set_ylabel('Loss (dB)')
#ax0.plot(df_ic['wl'][peak_position_ic], peak_value_ic['peak_heights'], 'x', linewidth=2, color='red')
ax0.legend(frameon=False)

ax1.plot(df_ic['wl'][mid_pos_ic], fsr_ic, 'o-', label='FSR (model)', linewidth=2, color='red')
ax1.set_xlabel('Wavelength (nm)')
ax1.set_ylabel('fsr (nm)')
#ax1.set_ylim([10, 18])
ax1.legend(frameon=False)

ax2.plot(df_ic['wl'][mid_pos_ic], ng_ic, 'o-', label='$n_g$ (model)', linewidth=2, color='red')
ax2.set_xlabel('Wavelength (nm)')
ax2.set_ylabel('Group index')
ax2.set_ylim([3, 5])
ax2.legend(frameon=False)

ax3.plot(lam_ic*1e9, neff_ic, 'o-.', label='$n_{eff}$ (model)', linewidth=2, color='red')
ax3.set_xlabel('Wavelength (nm)')
ax3.set_ylabel('Effective index')
#ax3.set_ylim([0, 5])
ax3.legend(frameon=False)

ax4.plot(df_ic['wl'][mid_pos_ic], er_ic, 'o-', label='Extinction ratio (model)', linewidth=2, color='red')
ax4.set_xlabel('Wavelength (nm)')
ax4.set_ylabel('extinction ratio (dB)')
ax4.set_ylim([20, 40])
ax4.legend(frameon=False)

ax5.plot(df_ic['wl'][mid_pos_ic], Q_ic, 'o-', label='Q factor (model)', linewidth=2, color='red')
ax5.set_xlabel('Wavelength (nm)')
#ax5.set_ylim([250, 500])
ax5.set_ylabel('Q factor')
ax5.legend(frameon=False)

fig.savefig('3um_through_parameter_summary.png', bbox_inches='tight', dpi=200)



