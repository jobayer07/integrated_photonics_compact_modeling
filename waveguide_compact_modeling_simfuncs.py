import pandas as pd
import itertools
import sys, os
#sys.path.append("/opt/lumerical/v222/api/python") #Default windows lumapi path
#sys.path.append(os.path.dirname(__file__)) #Current directory
#import lumapi
import json
import numpy as np
from shutil import copyfile
import matplotlib.pyplot as plt

def flat(data):
    return list(itertools.chain(*data))

def savedf(dataframe=pd.DataFrame(),datadir="./",datafile="foo.csv"):
    if not os.path.isdir(datadir):
        print(datadir," DOES NOT EXIST...CREATING")
        os.mkdir(datadir)
    fullpath = os.path.join(datadir, datafile)
    print('Saving: ', fullpath)
    dataframe.to_csv(fullpath,index=False)

def onarrsimdata(lumhandle,ona="ONA_1",onainput="input 1"):
    # data              = s.getresult("ONA_1",'input 1/mode 1/gain');
    onadatastr = onainput + '/mode 1/gain'
    data = lumhandle.getresult(ona,onadatastr);
    df = pd.DataFrame()
    # df['wavelength']  = sf.flat(data['wavelength'])
    df['gain'] = data['TE gain (dB)']
    # Wavelength is not a flat list. Either add column after gain or use sf.flat to flatten first
    df['wavelength'] = data['wavelength']

    # Available parameters
    # input 1/mode 1/peak/transmission     # input 1/mode 1/peak/frequency     # input 1/mode 1/peak/angle
    # input 1/mode 1/peak/group delay    # input 1/mode 1/peak/group velocity    # input 1/mode 1/peak/dispersion
    # input 1/mode 1/peak/loss    # input 1/mode 1/peak/gain    # input 1/mode 1/peak/bandwidth
    # input 1/mode 1/peak/free spectral range    # input 1/mode 1/peak/quality factor    # input 1/mode 1/peak/finesse
    # input 1/mode 1/peak/extinction ratio
    pkonaparam     = ["frequency","quality factor","extinction ratio","gain","free spectral range"]
    pkonaparam_var = ["wavelength", "Q", "ER", "gain","fsr"]
    dfpk = pd.DataFrame()

    for param,var in zip(pkonaparam,pkonaparam_var):
        pkdata    = lumhandle.getresultdata(ona, onainput + '/mode 1/peak/' + param)
        if isinstance(pkdata,float) or isinstance(pkdata,int):
            print('NOTE: NO PEAK FOUND OR JUST ONE PEAK FOUND FOR ',var, " IN ",onainput)
            dfpk[var] = pkdata
        else:
            dfpk[var] = flat(pkdata)


    #dfpk['wavelength'] = sf.flat(s.getresultdata("ONA_1", 'input 1/mode 1/peak/frequency'));
    #dfpk['Q'] = sf.flat(s.getresultdata("ONA_1", 'input 1/mode 1/peak/quality factor'));
    #dfpk['ER'] = sf.flat(s.getresultdata("ONA_1", 'input 1/mode 1/peak/extinction ratio'));
    #dfpk['gain'] = sf.flat(s.getresultdata("ONA_1", 'input 1/mode 1/peak/gain'));

    return df, dfpk

def read_wg_json(ipjson):
    # Import the input json file
    f = open(ipjson)
    data = json.load(f);

    # Modify parameters
    # 1) neff real and imaginary 2) width 3) temperature
    # 4) D 5) loss vs. wavelength 6) mode_data 7) ng

    # First change parmeters that have a single value
    data_op = data.copy();

    wavelength_data = data['wavelength_data']['_data']
    neff            = data['neff']['_data']
    D               = data['D']['_data']
    ng              = data['ng']['_data']


    f.close()

    #rtndict = {'wavelength_data':wavelength_data,'neff':neff,'ng':ng,'D':D}

    df = pd.DataFrame()
    df['wavelength_data'] = data['wavelength_data']['_data']
    df['neff']            = data['neff']['_data']
    df['ng']              = data['ng']['_data']
    df['D']               = data['D']['_data']

    return df

def list_raised_to(list,pow):
    return [i ** pow for i in list]

def diff2(x,y):
  n = len(x)-1

  # Taking linear extrapolation of last 2 points
  dydx = np.diff(y)/np.diff(x)

  nd = len(dydx)-1

  m_slope = (dydx[nd] - dydx[nd-1])/(x[nd] - x[nd-1])
  y_int   = dydx[nd] - m_slope*x[nd]

  print(x[n])
  print(dydx.shape)
  print(x.shape)
  print(m_slope)
  dydx = np.append(dydx,m_slope*x[n] + y_int)

  return dydx

def neff_func(w,w0,n1,n2,n3,opt=0):
    neff_calc = n1 + n2*(w - w0) + n3*(w - w0)**2
    if opt:
        return neff_calc
    else:
        return neff_calc

def medianp(xarray):
    return(100*(xarray - np.median(xarray))/np.median(xarray))


def installdk(lumhandle,dkpath='foo.cml',dkinstpath='.'):
    lumhandle.switchtolayout()
    lumhandle.installdesignkit(dkpath, '.', True)

def copy_file_temp(srcfile='foo.txt',tmpext='_work'):
    # foo.txt is now saved as foo_work.txt
    dstfile = os.path.join(os.path.dirname(srcfile), os.path.basename(srcfile).replace('.', '_work.'))
    print("Copying: ",srcfile, '  to  ', dstfile)
    copyfile(src=srcfile, dst=dstfile)
    return dstfile


def read_rr_json(ipjson):
    # Import the input json file
    f = open(ipjson)
    print("READING JSON FILE: ", ipjson)
    data = json.load(f);
    print(data)
    print("....DONE...")

    # wavelength_data = data['wavelength_data']['_data']
    # neff            = data['neff']['_data']
    # D               = data['D']['_data']
    # ng              = data['ng']['_data']


    f.close()

    #rtndict = {'wavelength_data':wavelength_data,'neff':neff,'ng':ng,'D':D}

    # df = pd.DataFrame()
    # df['wavelength_data'] = data['wavelength_data']['_data']
    # df['neff']            = data['neff']['_data']
    # df['ng']              = data['ng']['_data']
    # df['D']               = data['D']['_data']

    return data

def write_rr_json(opjson,data_op):
    f_op = open(opjson, 'w')
    print('WRITING json FILE: ',opjson)
    json.dump(data_op, f_op, indent=4, ensure_ascii=False)
    f_op.close()

# def modify_json_v0p0(jsondict,jsonparams):
#     cnt = 0
#     for k1, v1 in jsonparams.items():
#         # Only proceed if the key exists in jsondict
#         if not is_key_present(k1,jsondict):
#             continue
#         print("KEY: ",k1)
#         if isinstance(v1, dict):
#             print("FOUND dictionary: ",v1)
#             for k2, v2 in v1.items():
#                 if not is_key_present(k2, jsondict[k1]):
#                     continue
#                 cnt = cnt + 1
#                 print("#",cnt," Setting: ", k1, ' ', k2, ' to ', v1[k2])
#                 jsondict[k1][k2] = v1[k2]
#         else:
#             cnt = cnt + 1
#             print("#",cnt," Setting: ", k1, ' to ', v1)
#             jsondict[k1] = v1
#
#     return jsondict


def modify_json(jsondict,jsonparams):
    print('MODIFYING json DICTIONARY')
    cnt = 0
    if len(jsonparams) == 0:
        print("EMPTY PARAMETER DICTIONARY FOUND .. NOT CHANGING ANYTHING")
        return jsondict
    for k1, v1 in jsonparams.items():
        # Only proceed if the key exists in jsondict
        if not is_key_present(k1,jsondict):
            print("KEY: ", k1," NOT IN json file ... skipping")
            continue
        else:
            print("CHANGING VALUE of KEY: ",k1,' to ',v1)
            jsondict[k1] = v1
    return jsondict

def is_key_present(input_key,thisdict):
    if input_key in thisdict.keys():
        return True
    else:
        print("Key: ", input_key, " is not present in ", thisdict)
        return False


def setuprrsim(intcfile):
    # Runs only for a single mode. Set input to one port and get outputs from other ports
    model_name = row_df['model_name'];
    # Connects inputs from right ot left of list connected to input 1, input 2, input 3 etc
    portlist = row_df['ports'].split('|')
    #portnums = range((len(portlist),0,-1)
    wavelengths = [float(x) for x in row_df['lambda_start_stop'].split("|")]
    inc = lumapi.INTERCONNECT(hide=True)
    print("Setting up ONA...")
    # Add network analyzer connect inputs and outputs
    inc.addelement('Optical Network Analyzer')
    inc.setnamed("ONA_1", "number of input ports", len(portlist));
    inc.setnamed("ONA_1", "input parameter", "start and stop");
    inc.setnamed("ONA_1", "start frequency", c * 1e9 / wavelengths[0]);
    inc.setnamed("ONA_1", "stop frequency", c * 1e9 / wavelengths[1]);
    inc.setnamed("ONA_1", "number of points", 2000);
    inc.setnamed("ONA_1", "plot kind", "wavelength");
    # Set the mode for ONA based
    if (row_df['mode'] == 'TM'):
        print("Setting ONA mode to: ",row_df['mode'])
        inc.setnamed("ONA_1",'orthogonal identifier',2)
    inc.setposition("ONA_1", 200, 100);
    print("ONA setup complete")
    # Now instantiate the device
    print("Instantiate " + model_name)
    inc.addelement(row_df['model_name'])

    instname = row_df['inst_prefix'] + "_1"
    inc.setposition(instname, 250, 300);
    if (row_df['template'].startswith('wg_te')):
        wg_length = 0.01 # 1 cm
        print("Waveguide element detected: "+model_name)
        print("Setting waveguide length = "+ str(wg_length) +' cm')
        inc.setnamed(instname, "wg_length", wg_length);
    print("Connecting ports...","Input: ",inport)
    inc.connect("ONA_1", "output", instname, inport)
    # ONA Input array from string represenation of '1' to 'numports'
    # output is ['input 1','input 2'...] which the ONA needs
    onaip        = ['input ' + str(i + 1) for i in range(len(portlist))]
    #ona_port_str = ['input ' + str(i + 1) + '_' + portlist[i] for i in range(len(portlist))]
    for cnt,op in zip(onaip,portlist):
        print("Connecting ports: ",op)
        #inc.connect("ONA_1", "input " + cnt, instname, op)
        inc.connect("ONA_1",cnt, instname, op)

    qaschematic =os.path.join( row_df['kit_sim_dir'],model_name+'_'+row_df['mode']+'_qa.icp')
    print("Running simulation for: "+model_name+"........")
    inc.run()

    # Create a string called "input 1/mode 1/transmission"
    # And check which ports have data available
    onadatastr = [x + '/mode 1/transmission' for x in onaip]

    # Check which ports have results in them. Lumerical ONA does not have any data for Sxx=0
    #onadatachk = [inc.haveresult("ONA_1",x) for x in onadatastr]

    # Now get data for the ports that have data
    # Prepare a dictionary of ONA inputs and model ports
    #portdict = dict(zip(portlist,onadatastr))

    dfsim = pd.DataFrame()
    for thisport,thisinput in zip(portlist,onadatastr):
        if inc.haveresult("ONA_1", thisinput):
            # Can use getresultdata to a specific data instead of dictionary
            # However that results in several calls for wavelength, transmission, angle etc
            #Tdata  = inc.getresultdata("ONA_1", "input 1/mode 1/transmission");

            print("Getting Data for: " + thisinput)
            data      = inc.getresult("ONA_1", thisinput);

            #dflocalsim  = pd.DataFrame(
            #    {
            #        'wavelength':data['wavelength'],
            #        'T'         :abs(data['TE transmission']),
            #        'ph'        :inc.angle(data['TE transmission'])
            #    }
            #
            #)
            dflocalsim = pd.DataFrame()
            # Transmission should be added to the data frame first.
            # Else use np.ravel to flatten the array
            # The array is returned as [[xxx][xxxx].. which has to be flattened
            # T is returned as flat but not sure. And after T other variables are coerced.
            mode_tr_str = row_df['mode'] + ' transmission'
            dflocalsim['wavelength']   = np.ravel(data['wavelength'])
            #dflocalsim['T']            = abs(data['TE transmission'])
            #dflocalsim['ph']           = inc.angle(data['TE transmission'])
            dflocalsim['T']            = abs(data[mode_tr_str])
            dflocalsim['T_r']          = np.real(data[mode_tr_str])
            dflocalsim['T_i']          = np.imag(data[mode_tr_str])
            dflocalsim['ph']           = inc.angle(data[mode_tr_str])
            dflocalsim['gain']         =  10 * np.log10(abs(dflocalsim['T']) ** 2)
            dflocalsim['input']        = inport
            dflocalsim['output']       = thisport
            dflocalsim['ona_port_str'] = thisport + '_' + thisinput
            dflocalsim['model_name']   = model_name
            dflocalsim['mode']         = row_df['mode']
            dfsim = dfsim.append(dflocalsim)

    inc.switchtolayout()
    print("Saving schematic: " + qaschematic)
    inc.save(qaschematic);
    inc.close()
    return(dfsim)

def pltsubplots(data_frame=pd.DataFrame(),pltmatrix=[2,2],xparam='wavelength',yparams=['gain','Q'],pltstyle='--o'):
    p2 = plt
    p2.figure()
    cnt = 0;
    for param in yparams:
        cnt = cnt + 1
        p2.subplot(pltmatrix[0], pltmatrix[1], cnt)
        if isinstance(data_frame,pd.DataFrame):
            p2.plot(data_frame[xparam],data_frame[param],pltstyle)
        else:
            # Plot all the data-frames in the list
            for thisdf in data_frame:
                p2.plot(thisdf[xparam],thisdf[param],pltstyle)

        p2.grid()
        p2.title(param)
    return p2

def changedevparam(lumhandle,devname,param,val):
    #lumhandle.setnamed("ONA_1", "number of input ports", len(portlist));
    lumhandle.setnamed(devname,param,val);
    return lumhandle

def switch2layoutandrun(lumhandle):
    lumhandle.switchtolayout();lumhandle.run()
    return lumhandle
# class pltdfsubplot:
#     def __init__(self,data_frame,pltmatrix,xparam,yparams):
#         self.data_frame = data_frame
