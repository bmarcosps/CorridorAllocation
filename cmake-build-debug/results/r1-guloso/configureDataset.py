import sys
import numpy as np 
import pandas as pd

import os
import glob

df_mean = pd.DataFrame(columns=['mean_score','mean_time'])
df_best = pd.DataFrame(columns=['best_score','best_time'])

for filename in glob.glob('*.txt'):
    print(filename)
    df = pd.read_csv(filename, header=None)
    df1 = pd.DataFrame({"mean_score":[df[0].mean(0)], 
                        "mean_time":[df[1].mean(0)]}) 
    df2 = pd.DataFrame({"best_score":[df[0].min(0)], 
                        "best_time":[df[1].min(0)]}) 

    df_mean= df_mean.append(df1, ignore_index=True)
    df_best= df_best.append(df2, ignore_index=True)
    
    
   #with open(os.path.join(os.cwd(), filename), 'r') as f:

print(df_mean)
print(df_best)


df_mean.to_csv('.\\mean.txt', ',', index=False)
df_best.to_csv('.\\best.txt', ',', index=False)