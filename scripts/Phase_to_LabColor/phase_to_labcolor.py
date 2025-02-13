import numpy as np
from colormath.color_objects import LabColor, sRGBColor
from colormath.color_conversions import convert_color

def phase_to_labcolor( φ:float, amp=[], th_amp=1, a_f=90, b_f=80 ,a_0=10 , b_0=20, l_p=45, l_m=25, obs='2', ill='d65'):

    # Function to convert phase to Lab color
    # The function converts a phase signal to a Lab color signal
    # 
    # Input:
    # φ: phase in radians
    # amp: amplitude of the signal
    # th_amp: threshold amplitude
    # a_f: factor for a
    # b_f: factor for b
    # a_0: offset for a
    # b_0: offset for b
    # l_p: luminance for positive values
    # l_m: luminance for negative values
    # obs: observer (2 or 10)
    # ill: illuminant (d65 or d50)
    #
    # Returns: Lab color
    
    # get dimension
    n = φ.shape[0]

    # convert phase to Lab color
    a = (np.cos(φ))*a_f + a_0
    b = (np.sin(φ))*b_f + b_0
    L = np.ones(n)
    idx = np.sin(φ)>0
    L[idx] = 53+np.sin(φ[idx])*l_p
    L[~idx] = 53-np.sin(3*φ[~idx])*l_m

    # convert Lab color to sRGB
    lab = np.full([n,3],np.nan)
    for i in range(n):
        col = convert_color(LabColor(lab_l=L[i],lab_a=a[i],lab_b=b[i],observer=obs, illuminant=ill),sRGBColor)
        lab[i,:] = np.array(col.get_value_tuple())

    # For low amplitude values, reduce the saturation of the color by scaling the color vector by the amplitude value below the threshold divided by the threshold amplitude
    if not len(amp) == 0:
        idx_lt1 = amp<th_amp
        factor = amp[idx_lt1]/th_amp
        lab[idx_lt1,:] = lab[idx_lt1,:]*np.repeat(factor[:,None],3,axis=1)

    # Normalize the color values to be between 0 and 1
    if np.any(lab>1):
        for j in range(3):
            if np.nanmax(lab[:,j])>1:
                lab[:,j] /= np.nanmax(lab[:,j])
    
    # Return the Lab color
    return lab



