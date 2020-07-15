# -*- coding: utf-8 -*-
#===============================================================================
import pandas as pd
from main_file import *
import numpy as np

def create_data_frame(data):
    data_frame = pd.DataFrame({'angle': data[:,0]*2/np.pi,
                           'eps': data[:,1],
                           't_1': data[:,2],
                           't_2': data[:,3],
                           't_3': data[:,4],
                           't_4': data[:,5]})
    data_frame.to_csv("diffusion_time.csv")

create_data_frame(data_time)  
