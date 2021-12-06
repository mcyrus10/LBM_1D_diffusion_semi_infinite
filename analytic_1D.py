#!/home/cyrus/miniconda3/bin/python
"""
This is the basic 1D diffusion solver from Alexiades' course. This script plots
the analytic solution with the LBM overlay
"""

import numpy as np
from scipy.special import erfc
import matplotlib.pyplot as plt

def read_xml(   f_name : str, fields : list) -> dict:
    """
    parser for the xml parameters file
    Parameters
    ----------
    f_name : string
        file name of parameters file (supply .xml extension)

    fields : list
        list of tuples with the first element as a string identifying the item
        and the second as a datatype (e.g. float, int, etc.)

    Returns:
    --------
    ret_dict : dictionary
        dictionary with keys as the entries from the parsed .xml file

    """
    retDict = {}
    with open(f_name,'r') as f:
        text = f.read().split("\n")
    for line in text:
        for field in fields:
            if "<{}>".format(field[0]) in line:
                temp = list(filter(None,line.split(" ")))
                retDict[field[0]] = field[1](temp[1])
    return retDict

def read_lbm_file(f_name : str) -> np.array:
    """
    This function reads the *.dat file that Palabos writes on the

    Parameters
    ----------
    f_name : string
        file name of c++ output

    Returns
    -------
    lbm_out : np.array
        array of c++ output

    """
    with open(f_name,'r') as f:
        text = f.read().split(" ")
    lbm_out = np.array([float(j) for j in text[0:-1]]).reshape(1,len(text)-1)
    return lbm_out
    
def u(x : np.array,t : float,D : float) -> np.array:
    """
    Analytic solution to 1D diffusion with zero concentration initially then
    constant boundary of 1 at time > 0

    Parameters
    ----------
    x : array-like or float
        x-coordinate
    
    t : float
        time

    D : float
        Diffusivity

    Returns
    -------
    solution to the PDE (array-like or float)
   
    """
    return erfc(x/(2*(D*t)**(1/2)))

if __name__ == "__main__":
    fields = [
                ('lx',int),
                ('ly',int),
                ('resolution',int),
                ('tau_ad',float),
                ('c_init',float),
            ]
    params = read_xml("params.xml",fields)

    resolution = params['resolution']
    nx = int(params['lx']*resolution)
    ny = int(params['ly']*resolution)
    tau = params['tau_ad']

    b = 1.0
    mm = resolution
    D = (1 / 3) * (tau - 0.5)
    t_end = 1
    x = np.linspace(0,b,mm)

    lbm_data=read_lbm_file("concentration_final.dat").flatten()
    x2 = np.linspace(0,b,len(lbm_data))

    plt.figure(0)
    ax = plt.gca()
    ax.plot(x,u(x,t_end,D),linewidth = .5, label='analytic')
    ax.plot(x2,lbm_data,'k.',linewidth = .5, label ='lbm')
    print("value at max x value = ",u(b,t_end,D))
    ax.legend()
    ax.set_xlabel('x')
    ax.set_ylabel('u(x,1)')
    plt.title("D = {:.2f}".format(D))
    plt.savefig("1d_diffusion.png",dpi = 300)
    plt.show()
