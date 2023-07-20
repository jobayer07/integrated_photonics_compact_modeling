import pandas as pd
import numpy as np

portnames = ['opt_1_fbr','opt_2_fn']
portlocs  = ['LEFT', 'RIGHT']
polarization='TE'
mode_id=1

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

convert_to_s_parameter(loss_model)
