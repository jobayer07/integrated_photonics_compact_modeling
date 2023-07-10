import numpy as np
import matplotlib.pyplot as plt
import imp
import pandas as pd
from scipy.signal import find_peaks


def pairwiseAvg(lst):
    size = len(lst)
    x=np.zeros(size-1)
    for i in range(len(lst)-1):
        x[i] = int((lst[i] + lst[i + 1])/2)
    return x

def find_index_of_nearest_value(array, value):
    idx = (np.abs(array - value)).argmin()
    return idx

def baseline_correction(df, deg):
    loss=-df['Normalized Transmission (dB)']
    poly = np.polyfit(df['Wavelength (nm)'], loss, deg=deg)
    loss_trend = np.polyval(poly, df['Wavelength (nm)'])        #finding the trend using 3rd order polynomial fitting

    trend_adjustment=min(loss_trend)-min(loss)
    loss_trend_shifted =  loss_trend-trend_adjustment     #adjusting the losser trend
    loss_base_corr=loss-loss_trend_shifted
    df['Normalized Transmission (dB)']=-loss_base_corr
    return df

def find_mid_index_er(d):
    peak_position, peak_value = find_peaks(d, height=2, distance=50)
    mid_index=pairwiseAvg(peak_position)  
    
    pairwisePeak=pairwiseAvg(peak_value.get('peak_heights'))
    mid_value=d[mid_index]
    
    er=pairwisePeak-mid_value                                   #extinction ratio
    return mid_index, er

def q_calculation(wl_all, loss, mid_index, peak_index):
    a=loss.values
    wl_all=model_df['wl']
    power_linear=10**(-a/10)
    fwhm_value_linear=(power_linear[mid_index.astype(int)]+power_linear[peak_index.astype(int)])/2
    fwhm_line=-10*np.log10(fwhm_value_linear)
    
    fwhm=np.zeros(len(fwhm_line))
    wl=np.zeros(len(fwhm_line))
    
    for i in range(len(fwhm_line)):
        y=a[peak_index[i]:int(mid_index[i])]
        x=wl_all[peak_index[i]:int(mid_index[i])]
        x=x.to_numpy()
        
        start=wl_all[peak_index[i]]
        index=find_index_of_nearest_value(y, fwhm_line[i])       #Now found the index of the fwhm point
        end=x[index]
        fwhm[i]=2*abs(end-start)
    
        wl[i]=wl_all[peak_index[i]]
    
    q=wl/fwhm
    return q

#----------------------------User Input--------------------
wg_name='tlx_cl_te_mrr_airclad_fn_internal'
N=100001
ring_radius=(225)*1e-6
ring_cir=2*np.pi*ring_radius

measured_df=pd.read_csv(r"C:\Users\mh220218\Documents\My Assignments\Compact Model\Ring_resonator\NRL Data\TLX_PDK1.0_MRR_C_L_band\tlx_cl_te_mrr_airclad_fn\tlx_cl_te_mrr_airclad_fn.txt", delimiter='\t')
measured_Analysis_df=pd.read_csv(r"C:\Users\mh220218\Documents\My Assignments\Compact Model\Ring_resonator\NRL Data\TLX_PDK1.0_MRR_C_L_band\tlx_cl_te_mrr_airclad_fn\tlx_cl_te_mrr_airclad_fn_Analysis.txt", delimiter='\t')

#-----------------------


directional_coupler=1
if directional_coupler:
    #fdtd_coupling = pd.read_csv(r"TLX_PDK1.0_MRR_C_L_band\tlx_cl_te_mrr_airclad_fn\tlx_fn_ring_arirclad_crossPort_T.txt", sep=", ")
    #lam=fdtd_coupling ['lambda(m)']
    lam=measured_Analysis_df['Wavelength (nm)']*1e-9
    
    lam0=1550e-9
    m0=0.6
    m1=1.25e7
    m2=0

    x=(lam-lam0)*8.71e7
    coupling_power_model= 0.77/(1+0.382*np.exp(-x)) +0.08 #m0 + m1*(lam-lam0) + m2*(lam)**2   0.77/(1+0.36*np.exp(-x)) +0.08

    
    plt.figure(1)    
    #plt.plot(fdtd_coupling ['lambda(m)']*1e9, fdtd_coupling ['Y'], '.b', markerfacecolor = 'k', label='Power coupling FDTD', linewidth=2)
    plt.plot(measured_Analysis_df['Wavelength (nm)'], measured_Analysis_df['Ring Power Coupling (fraction)'], 'ob', markerfacecolor = 'k', label='Power coupling (measurement)', linewidth=1)
    plt.plot(lam*1e9, coupling_power_model, '-r', linewidth=2, label='Power coupling Model')
    plt.xlabel('Wavelength ($\mu$m)')
    plt.ylabel('Power coupling (fraction)')
    plt.legend(frameon=False)

    model_dc_power_coupling_coefficient={'f':3e8/lam.values, 'pow': coupling_power_model.values}
    model_dc_power_coupling_coefficient = pd.DataFrame(model_dc_power_coupling_coefficient)
    model_dc_power_coupling_coefficient.to_csv('model_dc_power_coupling_coefficient.txt', header=None, index=None, sep=' ') 
    
#---------------------------------------------------------Interconnect------------------------------------------
ring_simulation=1
if ring_simulation:
    lumapi=imp.load_source("lumapi", "C:\\Program Files\\Lumerical\\v231\\api\\python\\lumapi.py")
    inc = lumapi.INTERCONNECT(filename="rr_schematic_single_bus.icp") #hide=True
    
    #inc.switchtolayout()
    inc.switchtodesign()
    
    #Updating the directional couplers
    inc.select('DC1')
    inc.set('measurement filename 1','model_dc_power_coupling_coefficient.txt')
    
    #Install design kit for waveguide
    inc.installdesignkit (r"C:\Users\mh220218\Documents\My Assignments\Compact Model\Ring_resonator\NRL Data\rr_model_cml\input_files\rr_model_v0p1\artifacts\interconnect\rr.cml", "C:/Users/mh220218/Documents/My Assignments/Compact Model/Ring_resonator/NRL Data/rr_model_cml/installation_not_needed", True);
    
    #inc.setnamed('::Root Element','start frequency', 1591)
    #inc.addelement("Optical Oscilloscope")
    
    
    inc.select("::Root Element::ONA_1")
    inc.set("number of points", N)
    
    #Adding waveguides
    inc.addelement(str(wg_name))
    inc.set("name","WG_SiN_1")
    inc.set("wg_length", ring_cir)

    #Connection
    inc.connect("WG_SiN_1","opt_1","DC1","port 2")
    inc.connect("WG_SiN_1","opt_2","DC1","port 4")

    inc.run(3)                           #1:single processor mode. 2: single processor mode, Pop-up dialogs no longer take focus. 3: parallel mode as defined in the resource manager
    ring_model = inc.getresult("ONA_1", "input 1/mode 1/gain")
    
    fsr_model= inc.getresult("ONA_1", "input 1/mode 1/peak/free spectral range")
    q_interconnect= inc.getresult("ONA_1", "input 1/mode 1/peak/quality factor")
    
    # ------After simulation is done--------
    # inc.switchtodesign()
    # inc.select("WG_SiN_1")
    # inc.delete()
    # inc.select("WG_SiN_2")
    # inc.delete()
    # inc.uninstalldesignkit("rr")
    # inc.save()
    # inc.close()
    
    '''
    # set('coupling coefficients 1', ([1,2], [3,4]));
    #set('measurement filename 1','test.txt');
    '''

#------------------------------------Calculation----------------
measured_df=baseline_correction(measured_df, deg=3)
mid_index_measure, er_measure=find_mid_index_er(-measured_df['Normalized Transmission (dB)'])
peak_index_measure, peak_value_measure = find_peaks(-measured_df['Normalized Transmission (dB)'], height=1, distance=5)


model={'wl':ring_model['wavelength'].flatten(), 'loss': -ring_model['TE gain (dB)'].flatten()}
model_df=pd.DataFrame.from_dict(model)
model_df.to_csv('model_transmission_df.csv')

peak_index_model, peak_value_model = find_peaks(model_df['loss'], height=1, distance=1)
mid_index_model, er_model=find_mid_index_er(model_df['loss'])
q_model=q_calculation(model_df['wl'], model_df['loss'], mid_index_model, peak_index_model[:-1])


#----------------------------------Plotting--------------------------------
    
plot=1
if plot:
    fig, ([ax0, ax1], [ax2, ax3]) = plt.subplots(nrows=2, ncols=2, sharex=False, figsize=(20, 12))
    
    ax0.plot(ring_model['wavelength']*1e9, ring_model['TE gain (dB)'], linewidth=2, color='red', label='Model')
    ax0.plot(measured_df['Wavelength (nm)'], measured_df['Normalized Transmission (dB)'], '-', label='Measured data', linewidth=2, color='darkblue')
    ax0.set_xlabel('Wavelength (nm)')
    ax0.set_ylabel('Through port power (dB)')
    #ax0.set_xlim([1500, 1590])
    #ax0.set_xlim([1540, 1560])
    #ax0.set_ylim([-5, 0])
    ax0.legend(frameon=False)
    ax0.tick_params(axis='both', which='major', labelsize=14)
    #ax0.set_yticks([0, -1, -3, -5, -10, -15, -20], [0, -1, -3, -5, -10, -15, -20])
    
    ax1.plot(measured_Analysis_df['Wavelength (nm)'], measured_Analysis_df['FSR (nm)'], '.', linewidth=2, color='darkblue', label='FSR (measured)')
    ax1.plot(fsr_model['wavelength']*1e9, fsr_model['TE free spectral range (m)']*1e9, '-', linewidth=2, color='red', label='FSR (Model)')
    ax1.set_xlabel('Wavelength (nm)')
    ax1.set_ylabel('FSR (nm)')
    ax1.legend(frameon=False)
    ax1.tick_params(axis='both', which='major', labelsize=14)

    ax2.plot(measured_df['Wavelength (nm)'][mid_index_measure], er_measure, '.', linewidth=2, color='darkblue', label='Extinction ratio (measured)')
    ax2.plot(model_df['wl'][mid_index_model]*1e9, er_model, '-', linewidth=2, color='red', label='Extinction ratio (model)')
    ax2.set_xlabel('Wavelength (nm)')
    ax2.set_ylabel('Extinction ratio')
    ax2.legend(frameon=False)
    ax2.set_ylim([0, 25])
    ax2.tick_params(axis='both', which='major', labelsize=14)
    
    ax3.plot(measured_Analysis_df['Wavelength (nm)'], measured_Analysis_df['Q-factor'], '.', linewidth=2, color='darkblue', label='Q (measured)')
    #ax3.plot(wl*1e9, q, color='red')
    ax3.plot(model_df['wl'][mid_index_model]*1e9, q_model, '-', linewidth=2, color='red', label='Q (Model)')
    #ax3.plot(q_interconnect['wavelength']*1e9, q_interconnect['TE quality factor'], '-', linewidth=2, color='red', label='Q (Interconnect)')
    ax3.set_xlabel('Wavelength (nm)')
    ax3.set_ylabel('Q Factor')
    ax3.set_ylim([0, 5e4])
    ax3.legend(frameon=False)
    ax3.tick_params(axis='both', which='major', labelsize=14)
    
    plt.rc('axes', labelsize=20)
    plt.rc('legend', fontsize=18)
    
    plt.savefig('tlx_cl_te_mrr_airclad_fn2.tiff', bbox_inches='tight')

# plt.figure(3)
# plt.plot(measured_df['Wavelength (nm)'][peak_index_measure], model_df['wl'][peak_index_model])        
