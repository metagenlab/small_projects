

import pandas
from metagenlab_libs import gendb_utils
import datetime 
from django.db import IntegrityError

def calculate_age(birth_date, prel_date):
    from dateutil.relativedelta import relativedelta
    rdelta = relativedelta(prel_date, birth_date)
    age_years = rdelta.years 
    if age_years == 0:
        age_months = rdelta.months
        if age_months == 0:
            age_weeks = rdelta.weeks
            return f'{age_weeks} weeks'
        else:
            return f'{age_months} months'
    else:
        return age_years

def parse_molis_xml(XML_TABLE, SQL_DB_PATH):
    '''
<Cell ss:StyleID="th1"><Data ss:Type="Spath_listtring">N° de demande</Data></Cell> 
 <Cell ss:StyleID="th4"><Data ss:Type="String">Période</Data></Cell> 
 <Cell ss:StyleID="th3"><Data ss:Type="String">Numéro de demande</Data></Cell> 
 <Cell ss:StyleID="th4"><Data ss:Type="String">Patient hospitalisé</Data></Cell> 
 <Cell ss:StyleID="th3"><Data ss:Type="String">Demandeur</Data></Cell> 
 <Cell ss:StyleID="th4"><Data ss:Type="String">Sexe</Data></Cell> 
 <Cell ss:StyleID="th3"><Data ss:Type="String">Date de naissance</Data></Cell> 
 <Cell ss:StyleID="th4"><Data ss:Type="String">Numéro de patient</Data></Cell> 
 <Cell ss:StyleID="th3"><Data ss:Type="String">Unité de soins</Data></Cell> 
 <Cell ss:StyleID="th4"><Data ss:Type="String">Numéro de projet</Data></Cell> 
 <Cell ss:StyleID="th3"><Data ss:Type="String">Patient</Data></Cell> 
 <Cell ss:StyleID="th4"><Data ss:Type="String">Numéro alias</Data></Cell> 
 <Cell ss:StyleID="th3"><Data ss:Type="String">Référence externe</Data></Cell> 
 <Cell ss:StyleID="th4"><Data ss:Type="String">Numéro patient externe</Data></Cell> 
 <Cell ss:StyleID="th3"><Data ss:Type="String">Date saisie dem.</Data></Cell> 
 <Cell ss:StyleID="th4"><Data ss:Type="String">Heure de saisie de la demande</Data></Cell> 
 <Cell ss:StyleID="th3"><Data ss:Type="String">Date de réception</Data></Cell> 
 <Cell ss:StyleID="th4"><Data ss:Type="String">Heure de réception</Data></Cell> 
 <Cell ss:StyleID="th3"><Data ss:Type="String">Date prélèvement</Data></Cell> 
 <Cell ss:StyleID="th4"><Data ss:Type="String">Heure de prélèvement</Data></Cell> 
 <Cell ss:StyleID="th3"><Data ss:Type="String">Numéro de séjour</Data></Cell> 
 <Cell ss:StyleID="th4"><Data ss:Type="String">Date dern. édition</Data></Cell> 
 <Cell ss:StyleID="th3"><Data ss:Type="String">L&apos;heure du compte-rendu</Data></Cell> 
 <Cell ss:StyleID="th4"><Data ss:Type="String">Remarque interne</Data></Cell> 
 <Cell ss:StyleID="th3"><Data ss:Type="String">Remarque sur compte rendu</Data></Cell> 
 <Cell ss:StyleID="th4"><Data ss:Type="String">Renseignements cliniques</Data></Cell> 
 <Cell ss:StyleID="th3"><Data ss:Type="String">Statut validation niv. 2 (Méd.)</Data></Cell> 
 <Cell ss:StyleID="th4"><Data ss:Type="String">Matériel</Data></Cell> 
 <Cell ss:StyleID="th3"><Data ss:Type="String">Adresse</Data></Cell> 
 <Cell ss:StyleID="th4"><Data ss:Type="String">CodePostal</Data></Cell> 
 <Cell ss:StyleID="th3"><Data ss:Type="String">Ville</Data></Cell> 
 <Cell ss:StyleID="th4"><Data ss:Type="String">Canton-Pays</Data></Cell> 
 <Cell ss:StyleID="th3"><Data ss:Type="String">COVTYP</Data></Cell>
 <Cell ss:StyleID="th3"><Data ss:Type="String">PREL</Data></Cell>
    '''
    GEN_DB = gendb_utils.DB(SQL_DB_PATH)
    from xml.dom import minidom
    xmldoc = minidom.parse(XML_TABLE)

    itemlist=xmldoc.getElementsByTagName('Row')

    row_list = []
    for n,rows in enumerate(itemlist):
        item=rows.getElementsByTagName('Cell')
        if n == 0:
            columns = [cells.childNodes[0].childNodes[0].nodeValue if len(cells.childNodes[0].childNodes) > 0 else '' for cells in item]
        else:
            row_list.append([cells.childNodes[0].childNodes[0].nodeValue if len(cells.childNodes[0].childNodes) > 0 else '' for cells in item])

    df = pandas.DataFrame(row_list) 
    df.columns = columns
    
    return df 
    
    
