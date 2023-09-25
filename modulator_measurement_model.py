import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy import signal
import os
import imp, sys

#---------------------------------------------------------------Main Code----------------------------------------------------------------
def raw_data_clean(raw_data_file):
    df = pd.read_csv('Single BUS Carrier Injection 15umR/'+str(raw_data_file), skiprows=5311, header=None)
    df2=df.iloc[:, [0,5]]                                       #slice only the wavelength and TE_through port
    df3=df2.dropna()                                            #drop all rows having nan values
    df3 = df3.rename({0: 'wl', 5: 'loss'}, axis='columns')      #changing the column name
    df3['wl']=df3['wl']*1e9
    return df3

def grating_correction(df3):
    poly = np.polyfit(df3['wl'], df3['loss'], deg=3.8)
    df3['loss_trend'] = np.polyval(poly, df3['wl'])             #finding the trend using 3rd order polynomial fitting

    trend_adjustment=min(df3['loss_trend'])-min(df3['loss'])
    df3['loss_trend'] =  df3['loss_trend']-trend_adjustment     #adjusting the losser trend

    df3['loss_base_corr']=df3['loss']-df3['loss_trend']
    
    x = df3['wl'].values
    y = df3['loss_base_corr'].values
    return x, y

def base_trend_correction(y):
    b, a = signal.butter(3, 0.01, btype='lowpass', analog=False)
    base_trend = signal.filtfilt(b, a, y)

    condi=np.zeros(len(base_trend ))
    for i in range (len(base_trend )):
        if y[i]<base_trend [i]:
            condi[i]=base_trend [i]
        else:
            condi[i]=y[i]
    return condi

def high_pass_filtering(condi):
    b,a=signal.butter(3, 0.01, btype='highpass', analog=False)
    filtered = signal.filtfilt(b,a, condi)
    return filtered

def peak_masked_cleaning(y_high_pass_filtered, data_points_on_peak):
    peak_position, peak_value = find_peaks(y_high_pass_filtered, height=1, distance=90)
    mask=np.zeros(len(y_high_pass_filtered))
    for j in range (data_points_on_peak+1):
        mask[peak_position-j]=1
        mask[peak_position+j]=1
        
    clean_signal=y_high_pass_filtered*mask

    masked_index=np.argwhere(mask==0)
    masked_values=y_high_pass_filtered[masked_index]
    med=np.ma.median(masked_values)

    clean_signal[masked_index]=med
    clean_signal[clean_signal<med]=med
    return clean_signal

def bottom_zeroing(y):
    y=y-min(y)
    return y

def output_file_creation(x, y, x2, y2, x3, y3, x4, y4):
    clean_data = {'wl': x, 'loss': y, 'wl2': x2, 'loss2': y2, 'wl3': x3, 'loss3': y3, 'wl4': x4, 'loss4': y4}
    clean_data = pd.DataFrame(data=clean_data)
    clean_data.to_csv('clean_data_multiple.csv', index=False)
    return

def data_cleanup_one_modulator(raw_data_file):
    clean_raw=raw_data_clean(raw_data_file)
    x,y=grating_correction(clean_raw)
    y_base_corrected=base_trend_correction(y)
    #y_high_pass_filtered=high_pass_filtering(y_base_corrected)
    #y_masked_clean=peak_masked_cleaning(y_high_pass_filtered, 9)    #number of datapoints in the peak mask=7*2+1=15
    #y_bottom_zeroed=bottom_zeroing(y_masked_clean)
    return x, y_base_corrected

def filename_processing(directory):
    dir_list = os.listdir(directory)
    Name= {'Name': dir_list}
    files = pd.DataFrame(data=Name)

    files[['device', 'lot_id', 'lot_id2', 'wafer_id', 'test_algorithm', 'die', 'die2', 'voltage', 'current']] = files.Name.str.split("_", expand = True)
    files['voltage'] = files['voltage'].str.replace(r'p', '.')
    files['current'] = files['current'].str.replace(r'p', '.')
    files['current'] = files['current'].str.replace(r'.csv', '')
    #files['current'] = files['current'].str.replace(r'-0.000', '0.000')
    
    files['voltage']=pd.to_numeric(files['voltage'])
    files['current']=pd.to_numeric(files['current'])
    files['current']=abs(files['current'])

    files["lot_id"] = files["lot_id"].str.cat(files["lot_id2"], sep = "_")
    files["die"] = files["die"].str.cat(files["die2"], sep = "_")

    files=files.drop('die2', axis=1)
    files=files.drop('lot_id2', axis=1)
    return files

def data_cleanup_all(filename_df):
    x={}
    y={}
    for i in range(len(filename_df)):
        x['wl'+str(i)], y['loss'+str(i)]=data_cleanup_one_modulator(filename_df['Name'][i])
        
    x=pd.DataFrame(data=x)
    y=pd.DataFrame(data=y)
    output=pd.concat([x, y], axis=1).reindex(x.index)
    return output

#----------------------------------------------------------------Measurement-----------------------------------------------------------
directory=str('Single BUS Carrier Injection 15umR')

filename_df=filename_processing(directory)
clean_data=data_cleanup_all(filename_df)
clean_data.to_csv('clean_data_multiple.csv', index=False)

fig, ([ax0, ax1]) = plt.subplots(nrows=1, ncols=2, sharex=False, figsize=(12, 6))
ax0.set_xlabel('Wavelength (nm)')
ax0.set_ylabel('Loss (dB)')
ax0.set_xlim([1540, 1560])

ax1.plot(filename_df['voltage'], filename_df['current'], '-o', linewidth=2, color='darkblue')
ax1.set_xlabel('Voltage (V)')
ax1.set_ylabel('Current (mA)')
ax1.legend(frameon=False)

voltages=filename_df['voltage'].values
#len(voltages)
for i in range(8, len(voltages)):
    x, y=data_cleanup_one_modulator(filename_df['Name'][i])
    ax0.plot(x, y-1, '--', color='darkblue', linewidth=0.5, label= str(voltages[i])+' V') #color='darkblue',
    ax0.legend(frameon=False)    

#-------------------------------------------------------------------Model----------------------------------------------

sys.path.append("C:\\Program Files\\Lumerical\\v231\\api\\python\\") 
import importlib.util
spec_win = importlib.util.spec_from_file_location('lumapi', 'C:\\Program Files\\Lumerical\\v231\\api\\python\\lumapi.py')
lumapi = importlib.util.module_from_spec(spec_win) #windows
spec_win.loader.exec_module(lumapi)

#------------------------------------------------------------------------------------------------

r=15e-6

loss=220        #dB/m
coupling=[0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.07, 0.15, 0.17]
neff=[2.3, 2.3, 2.3, 2.3, 2.3, 2.3, 2.3, 2.3, 2.4206, 2.4244, 2.42765, 2.4282]
ng=3.86
D=0#-0.0192
dD=0

inc = lumapi.INTERCONNECT(filename="modulator_sim.icp")

inc.select("::Root Element::RING_2")
inc.set("length", 2*np.pi*r)
inc.set("loss", 120)
inc.set("group index", ng)
inc.set("dispersion", D)
inc.set("dispersion slope", dD)



for i in range(8, len(voltages)):
    inc.switchtodesign()
    
    inc.select("::Root Element::DC_1")
    inc.set("amplitude", voltages[i])
    
    inc.select("::Root Element::RING_2")
    inc.set("effective index", neff[i])
    inc.set("coupling coefficient 1", coupling[i])
    inc.set("coupling coefficient 2", coupling[i])
    
    inc.run(3)  
    
    gain_model = inc.getresult("ONA_1", "input 1/mode 1/gain")
    
    wl=gain_model['wavelength']
    gain=gain_model['TE gain (dB)']
    
    fsr_model= inc.getresult("ONA_1", "input 1/mode 1/peak/free spectral range")
    q_model= inc.getresult("ONA_1", "input 1/mode 1/peak/quality factor")
    
    wl_peak=fsr_model['wavelength']
    fsr=fsr_model['TE free spectral range (m)']
    q=q_model['TE quality factor']
    
    ax0.plot(wl[:, 0]*1e9, -gain, linewidth=1, label= str(voltages[i])+' V (model)')
    ax0.legend(frameon=False) 
    
inc.close()



















