import numpy as np
from colormath.color_objects import LabColor, sRGBColor
from colormath.color_conversions import convert_color

def phase_to_labcolor( φ:float, a_f=90, b_f=80 ,a_0=10 , b_0=20, l_p=45, l_m=25, obs='2', ill='d65'):

    n = φ.shape[0]

    a = (np.cos(φ))*a_f + a_0
    b = (np.sin(φ))*b_f + b_0
    L = np.ones(n)
    idx = np.sin(φ)>0
    L[idx] = 53+np.sin(φ[idx])*l_p
    L[~idx] = 53-np.sin(3*φ[~idx])*l_m

    lab = np.full([n,3],np.nan)
    for i in range(n):
        col = convert_color(LabColor(lab_l=L[i],lab_a=a[i],lab_b=b[i],observer=obs, illuminant=ill),sRGBColor)
        lab[i,:] = np.array(col.get_value_tuple())
    lab = lab/lab.max(axis=0)
    
    return lab