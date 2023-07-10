import pandas as pd
import numpy as np

portnames = ['opt_1_fbr','opt_2_fn']
portlocs  = ['LEFT', 'RIGHT']

N=76 #wavelength points
polarization='TE'
mode_id=1

#file_object  = open("filename", "mode")
f= open("first_s.txt","w+")

#S11
f.write('[\"' + portnames[0] + '\",\"' + portlocs[0] + '\"]\n')
f.write('[\"' + portnames[1] + '\",\"' + portlocs[1] + '\"]\n')
f.write('(\"' + portnames[0] + '\",\"' +str(polarization)+ '\",' +str(mode_id) + ',\"' + portnames[0] +'\",' + str(mode_id) +',\"transmission\")\n')
f.write('('+ str(N) +',3)')
f.close()