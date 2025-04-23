# -*- coding: utf-8 -*-
import os
import glob
import pandas as pd

folders = glob.glob("/data2/lishuhan/GSEdir/wget*")

dfs = []
for folder in folders:
    files = glob.glob(os.path.join(folder, "*.txt"))
    
    for file in files:
        with open(file, "r") as f:
            lines = f.readlines()
        
        data = {}
        for line in lines:
            if "=" in line:
                key, value = line.split("=", 1)
                if "Sample_characteristics_ch1" in key:
                    key, value = key + value.split(":")[0] + ":", value.split(":")[1] if ":" in value else value
                data[key.strip()] = value.strip()

        df = pd.DataFrame(data, index=[0])
        df['SAMPLE'] = os.path.basename(file).replace(".txt", "")
        dfs.append(df)

df_final = pd.concat(dfs, ignore_index=True)
df_final.to_csv("/data2/lishuhan/GSEdir/filter_data.csv", index=False)


