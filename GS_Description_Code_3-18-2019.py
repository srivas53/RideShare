# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 18:22:50 2019

@author: Varun
"""

import pandas as pd
import math
from time import time
import os


#read the table
greenscore_table = pd.read_csv("C:\\Users\\Varun Srivastava\\Desktop\\GC_Backend\\GC_Backend\\GreenScore Tables and Consolidation Code\\GS_GreenScores_all_systems_Nominated_brands_05-07.csv")

#greenscore_table =greenscore_table.sample(frac=0.12, random_state=123)
description_table = pd.read_csv("C:\\Users\\Varun Srivastava\\Desktop\\GC_Backend\\GC_Backend\\GS_Description_Model\\Supporting Data\\GS_Descriptions_Template_03_14.csv",encoding ='latin1')

systems=list(greenscore_table)
systems.remove('product_code_gc')
#systems.remove('Product_Code')
df = {'product_code_gc':[ ],'gs_system':[ ],'description_type':[ ],'description':[ ]}
Description = pd.DataFrame(df)


start1=time()
def Descriptions(x):
    for i in range(0,len(systems)):
        if math.isnan(x[systems[i]]):
            continue
        description_table2=description_table[description_table.gs_sub_system==systems[i]]
        sdlkmsdlkdsm
        kljlmm
        jksdn
        jfjjkdsssddssd
        for row2 in description_table2.itertuples(): 
            global Description
            if row2.score_lower<=x[systems[i]]<=row2.score_upper:
                df2= {'product_code_gc':x['product_code_gc'],'gs_system':row2.gs_system,'description_type':row2.description_type,'description':row2.description}
                df2 = pd.Series(df2)
                Description = Description.append(df2,ignore_index=True)        
    return (x)



greenscore_table=greenscore_table.apply(lambda x: Descriptions(x),axis=1)      
end1=time()
total=end1-start1       
    
Description=Description.drop_duplicates(['product_code_gc','description'],keep='first')
Description.to_csv('gs_description_table_all_except_nutrition_nominated_brands_5_08.csv',index = False)

path1=r'C:\Users\Varun Srivastava\Desktop\GC_Backend\GC_Backend\GS_Description_Model\GS_description_output_tables'
nut=pd.read_csv(os.path.join(path1,'nutrition_manual_data_04_29_descriptions.csv'))
#nut=nut.drop_duplicates(['product_code_gc','description'],keep='first')

nut = nut.reindex(sorted(nut.columns), axis=1)
Description = Description.reindex(sorted(Description.columns), axis=1)
final=pd.concat([Description,nut], axis=0)

final = final[pd.notnull(final['product_code_gc'])]
final=final.drop_duplicates(['product_code_gc','description'],keep='first')

final.to_csv(os.path.join(path1,'gs_description_output_nominated_brands_05_08.csv'))





















