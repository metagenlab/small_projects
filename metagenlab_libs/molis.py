#!/usr/bin/env python

import cx_Oracle
import os
import datetime 
from datetime import date


# Collection date       DW_MOLIS.T_MI_06 DATE_PREL
# Location
# Gender                DW_MOLIS.T_MI_06 SEXE
# Patient age           DW_MOLIS.T_MI_06 DATE_NAISS
# Patient status
# Specimen source       DW_MOLIS.RES


class Molis():
    
    def __init__(self, 
                 ):
        import sys
        # os.environ['MOLISPSW']
        pwd = 'larmes'
        # connection to the database
        self.connection = cx_Oracle.connect("etude_usr", pwd, "dw_molis_db")
        self.cursor = self.connection.cursor()

    def molis2patient_data(self,
                           molis_id,
                           num_bon=False):
        '''
        Either molis id or num bon
        '''
        LID = molis_id[0:2]
        LPERIOD = molis_id[2:6]
        ORDNB = molis_id[6:10]
        colnames = ["RES","NOM_PRENOM","DATE_NAISS","NO_HISTORIQUE", "DATE_PREL", "SEXE", "LID", "LPERIOD", "ORDNB", "NO_DEMANDEUR", "NOM_DEMANDEUR"] 
        updated_colnames = ["SMP_NAME","NOM_PRENOM","DATE_NAISS","PAT_ID", "DATE_PREL", "SEXE", "LID", "LPERIOD", "ORDNB", "NO_DEMANDEUR", "NOM_DEMANDEUR"]
        col_filter = ','.join(colnames)
        if num_bon:
            sql = f"""SELECT {col_filter} FROM DW_MOLIS.T_MI_06 
                      WHERE NO_BON='{molis_id}' and MC='PREL'"""

        else:
            sql = f"""SELECT {col_filter} FROM DW_MOLIS.T_MI_06 
                    WHERE LID='{LID}' and LPERIOD='{LPERIOD}' and ORDNB='{ORDNB}' and MC='PREL'"""
        #print(sql)
        self.cursor.execute(sql)
        data = {}
        for values in self.cursor:
            lst = [str(i) for i in values]
            for x,y in zip(updated_colnames, lst):
                data[x] = y
        
        data["DATE_NAISS"] = datetime.datetime.strptime(data["DATE_NAISS"], '%Y-%m-%d %H:%M:%S').strftime('%Y-%m-%d')
        
        # bon number should have 4 digits
        data["MOLIS_ID"] = f'{data["LID"]}{data["LPERIOD"]}{data["ORDNB"].zfill(4)}'

        year = '20' + data["LID"]
        month = data["LPERIOD"][0:2]
        day = data["LPERIOD"][2:5]

        data["DATE_REGISTERED"] = datetime.datetime.strptime(f"{year}-{month}-{day}", '%Y-%m-%d').strftime('%Y-%m-%d')

        # DATE_PREL might be missing
        if data["DATE_PREL"] not in [None, "None"] :
            print('PREL is NOT in [None "None"]', data["DATE_PREL"], type(data["DATE_PREL"]))
            data["DATE_PREL"] = datetime.datetime.strptime(data["DATE_PREL"], '%Y-%m-%d %H:%M:%S').strftime('%Y-%m-%d')
        else:
            print("PREL is NONE--------")

        return data



m = Molis()
print(m.molis2patient_data("8601654775", num_bon=True))
