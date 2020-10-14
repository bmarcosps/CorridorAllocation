import sys
import numpy as np 
import pandas as pd

import os
import glob

df_all = pd.DataFrame()

for filename in glob.glob('*.txt'):
    df = pd.read_csv(filename)
    df_all.append(pd)


print(df_all)
   #with open(os.path.join(os.cwd(), filename), 'r') as f:
        
    # open in readonly mode
    # do your stuff






#print(df)

# drop column by index
#df = df.drop(df.columns[0], axis=1)


# drop column by name
#df = df.drop('Subject index', 1)


#print(newDf)

#newDf.to_csv('.\\sensors2.txt', ' ', index=False)