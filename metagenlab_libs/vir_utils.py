#!/usr/bin/env python
from VFapp import settings
import os
import sys
import pandas

from BioSQL import BioSeqDatabase
from Bio import SeqIO
from Bio.SeqUtils import GC

import sqlite3

CLUSTERING_TABLE = settings.CLUSTERING_TABLE 
CLUSTERING2SPECIES_TABLE = settings.CLUSTERING2SPECIES_TABLE 

def quote(v):
    return f"\"{v}\""

class DB:
    def __init__(self, conn, cursor):
        self.server = cursor
        self.conn = conn
        self.db_path = settings.DB_PATH

    def load_db():

        db = settings.DB_PATH
        print("connecting to ", db)
        conn = sqlite3.connect(db)
        cursor = conn.cursor()
        
        return DB(conn, cursor)
    
    
    def get_VF_conservation_matrix(taxon_id, identity_cutoff=90, coverage_cutoff=80):
        
        sql = '''
        
        
        
        '''
    
        # convert long table to matrix with pandas
        
    def create_term_table(self,):
        sql = 'create table if not exists terms (term_id INTEGER PRIMARY KEY, term_name varchar(200))'
        self.server.execute(sql,)


    def add_term(self,term_name):
        sql = f'select term_id from terms where term_name="{term_name}"'
        try:
            term_id = self.server.execute(sql,).fetchall()[0][0]
        except IndexError:
            sql = f'insert into terms (term_name) values ("{term_name}")'
            self.server.execute(sql,)
            self.conn.commit()
            term_id = self.server.lastrowid
        return term_id

    def get_uniparc_entry_annotation(self, uniparc_id_list, add_source_db=True):
        
        uniparc_id_filter = ','.join([str(i) for i in uniparc_id_list])
        
        if add_source_db:
            # will duplicate row if the same uniparc entry is present in multiple dbs
            sql = f'''
            select distinct t3.db_name,t1.uniparc_accession,t4.cluster_name,t5.description from uniparc_entry t1
            inner join VF_table t2 on t1.uniparc_id=t2.uniparc_id
            inner join VF_databases t3 on t2.db_id=t3.db_id
            inner join {CLUSTERING_TABLE} t4 on t1.uniparc_id=t4.uniparc_id
            inner join uniparc_consensus_annotation t5 on t1.uniparc_id=t5.uniparc_id 
            where t1.uniparc_id in ({uniparc_id_filter})
            '''
        else:
             sql = f'''
            select distinct t1.uniparc_accession,t2.cluster_name,t3.description from uniparc_entry t1
            inner join {CLUSTERING_TABLE} t2 on t1.uniparc_id=t2.uniparc_id
            inner join uniparc_consensus_annotation t3 on t1.uniparc_id=t3.uniparc_id
            where t1.uniparc_id in ({uniparc_id_filter})
            '''           

        return self.server.execute(sql,).fetchall()
    
    def get_pmid_data(self,
                      pmid_list):
        
        
        pmid_filter = ','.join([str(i) for i in pmid_list])
        
        sql = '''
                select pmid,title,authors pmid2data where pmid in (?);
        '''
    
        return self.server.execute(sql,[pmid_filter]).fetchall()
    
    def get_VF_list_from_pmid(self, 
                              pmid):
        
        sql = f'''
                select distinct t3.db_name,t4.uniparc_accession,t5.cluster_name,t6.description from VF_id2pmid t1 
                inner join VF_table t2 on t1.VF_id=t2.VF_id 
                inner join VF_databases t3 on t2.db_id=t3.db_id
                inner join uniparc_entry t4 on t2.uniparc_id=t4.uniparc_id
                inner join {CLUSTERING_TABLE} t5 on t2.uniparc_id=t5.uniparc_id
                inner join uniparc_consensus_annotation t6 on t2.uniparc_id=t6.uniparc_id
                where t1.pmid=?;
        '''
        
        return list(self.server.execute(sql,[pmid]).fetchall())
    
    
    def get_pmid_info(self, 
                      pmid_list):
        
        pmid_filter = ','.join(['?']*len(pmid_list))
        
        sql = f'''
                select title,abstract,authors,source from pmid2data
                where pmid in ({pmid_filter});
        '''
        
        return list(self.server.execute(sql,pmid_list).fetchall())
    
    
    def get_species_taxids_with_comparative_data(self,):
        
        sql = 'select taxon_id from species_phylogeny'
        
        return [i[0] for i in self.server.execute(sql,).fetchall()]
    
    
    def get_db_VF_stats(self):
               
        sql = 'select db_name, count(*) from ' \
            ' (select distinct t2.db_name, t1.db_id, t2.db_name,t1.VF_id from VF_table t1 ' \
            ' inner join VF_databases t2 on t1.db_id=t2.db_id ' \
            ' inner join uniparc2species t3 on t1.uniparc_id=t3.uniparc_id where superkingdom_taxon_id=2) A' \
            ' group by A.db_id;'
        
        return {i[0]:i[1] for i in self.server.execute(sql,).fetchall()}

    def get_total_nr_VFs(self,):
        
        sql = 'select count(*) from (select distinct t1.uniparc_id from VF_table t1 ' \
        ' inner join uniparc2species t2 on t1.uniparc_id=t2.uniparc_id) A;'
        
        return self.server.execute(sql,).fetchall()[0][0]
        
    def get_total_nr_VFs_bacteria(self,):
        
        sql = '''select count(*) from (select distinct uniparc_id from (select * from VF_table t1 
               inner join uniparc2species t2 on t1.uniparc_id=t2.uniparc_id 
               where superkingdom_taxon_id=2) A  ) AA;'''
        
        return self.server.execute(sql,).fetchall()[0][0]


    def get_total_clusters(self,):
        
        sql = f'''select count(*) from (select distinct cluster_id from (select * from VF_table t1 
              inner join {CLUSTERING_TABLE} t2 on t1.uniparc_id=t2.uniparc_id) A ) AA; '''
        
        return self.server.execute(sql,).fetchall()[0][0]

    def get_total_clusters_bacteria(self,):
        
        sql = f'''select count(*) from (select distinct cluster_id from (select * from VF_table t1
             inner join {CLUSTERING_TABLE} t2 on t1.uniparc_id=t2.uniparc_id
             inner join uniparc2species t3 on t1.uniparc_id=t3.uniparc_id where superkingdom_taxon_id=2) A ) AA;
           ''' 
        
        return self.server.execute(sql,).fetchall()[0][0]

    def get_db_UP_stats(self):
               
        sql = 'select db_name, count(*) from (select distinct t2.db_name,t1.uniparc_id from VF_table t1 ' \
              ' inner join VF_databases t2 on t1.db_id=t2.db_id ' \
              ' inner join uniparc2species t3 on t1.uniparc_id=t3.uniparc_id where superkingdom_taxon_id=2' \
              ' group by t2.db_name,t1.uniparc_id)A group by db_name;'
        
        return {i[0]:i[1] for i in self.server.execute(sql,).fetchall()}
 
    def get_db_cluster_stats(self):
        
        sql = f'''select AA.db_name, count(*) from (select C.db_name,cluster_id from (
                   select distinct t1.uniparc_id,cluster_id from {CLUSTERING_TABLE} t1
                   inner join uniparc2species t2 on t1.uniparc_id=t2.uniparc_id where superkingdom_taxon_id=2) A
                   inner join VF_table B on A.uniparc_id=B.uniparc_id
                   inner join VF_databases C on B.db_id=C.db_id group by db_name,cluster_id) AA group by AA.db_name;
                   ''' 
        
        return {i[0]:int(i[1]) for i in self.server.execute(sql,).fetchall()}
 
        
    def get_taxon_id2n_VF_clusters(self, rank='species'):
        
        sql = f'''select {rank}_taxon_id, count(*) as n from (
          select distinct cluster_id, {rank}_taxon_id from VF_table t1
          inner join uniparc2species t2 on t1.uniparc_id=t2.uniparc_id
          inner join {CLUSTERING_TABLE} t3 on t1.uniparc_id=t3.uniparc_id 
          where species_name not like "%% sp.%%" 
          and species_name not like "%%uncultured%%"
          ) A
          group by A.{rank}_taxon_id order by n DESC;
          '''
        
        return {str(i[0]):i[1] for i in self.server.execute(sql,).fetchall()}


    def get_taxon_id2n_VF_clusters_separate_dbs(self, rank='species'):
        
        sql = f'''select db_name,{rank}_taxon_id, count(*) as n from (
          select distinct db_name,t3.cluster_id, t1.{rank}_taxon_id from uniparc2species t1 
          inner join {CLUSTERING_TABLE} t2 on t1.uniparc_id=t2.uniparc_id 
          inner join {CLUSTERING_TABLE} t3 on t2.cluster_id=t3.cluster_id
          inner join VF_table t4 on t3.uniparc_id=t4.uniparc_id
          inner join VF_databases t5 on t4.db_id=t5.db_id
          where t1.species_name not like "%% sp.%%" 
          and t1.species_name not like "%%uncultured%%"
          ) A
          group by A.db_name,A.{rank}_taxon_id order by n DESC;
          '''
        
        db2taxid2cout = {}
        for row in  self.server.execute(sql,).fetchall():
            db_name, taxon_id, count = row
            if db_name not in db2taxid2cout:
                db2taxid2cout[db_name] = {}
            db2taxid2cout[db_name][str(taxon_id)] = count
        
        return db2taxid2cout

    def get_db2taxon_id2VF_freq(self, rank='species', min_identity=90, min_coverage=60):
        
        sql = f'''
                    select AA.db_name,AA.species_taxon_id,AA.cluster_id,ifnull(BB.n, 0) from (
                select distinct db_name,cluster_id,species_taxon_id from uniparc2species t1
                inner join VF_table t2 on t1.uniparc_id=t2.uniparc_id
                inner join VF_databases t3 on t2.db_id=t3.db_id
                inner join {CLUSTERING_TABLE} t4 on t1.uniparc_id=t4.uniparc_id
            ) AA
            left join 
                (
                select B.taxon_id,B.cluster_id, count(*) as n from 
                    (select distinct t7.taxon_id,t5.assembly_id,cluster_id from homology_search t1
                    inner join homology_search_db_xrefs t2 on t1.search_tool_id=t2.id
                    inner join homology_search_db_xrefs t3 on t1.substitution_matrix_id=t3.id
                    inner join homology_search_db_xrefs t4 on t1.filtering_id=t4.id
                    inner join genome_assembly_table t5 on t1.assembly_id=t5.assembly_id
                    inner join uniparc2mmseqs_90_80 t6 on t1.uniparc_id=t6.uniparc_id
                    inner join genome_assembly_table2taxon_id t7 on t5.assembly_id=t7.assembly_id
                    where t1.percent_identity>={min_identity}
                    and uniparc_seq_coverage>={min_coverage}
                    and t2.name in ("ssearch", "tfastx") 
                    and t3.name="BL50" 
                    and t4.name="no_filtering_lowcomplexity"
                    ) B group by B.taxon_id, B.cluster_id
                ) BB
            on AA.cluster_id=BB.cluster_id
            and AA.species_taxon_id=BB.taxon_id;
          '''
        print(sql)
        db2taxid2cout = {}
        for row in  self.server.execute(sql,).fetchall():
            db_name, taxon_id, cluster_id, count = row
            if db_name not in db2taxid2cout:
                db2taxid2cout[db_name] = {}
            if str(taxon_id) not in  db2taxid2cout[db_name]:
                db2taxid2cout[db_name][str(taxon_id)] = {}
            db2taxid2cout[db_name][str(taxon_id)][cluster_id] = count
        
        return db2taxid2cout


    def get_VF_cluster2frequency_within_species2(self,
                                                 taxon_id_list, 
                                                 min_identity=90, 
                                                 min_coverage=60,
                                                 percentages=False):
        '''
        Work with species and genus taxids
        '''
        taxon_id_list = [str(i) for i in taxon_id_list]
        
        if percentages:
            taxid2n_genomes = {}
            for taxon_id in taxon_id_list:
                taxid2n_genomes[taxon_id] = self.get_n_genomes(taxon_id)
        
        print("taxid2n_genomes", taxid2n_genomes)
        
        taxid_filter = ','.join(taxon_id_list)

        # get cluster frequency
        # start from spacies table to count 
        # VF clusters with 0 hits on all genomes
        sql = f'''
            select AA.species_taxon_id,AA.cluster_name,ifnull(BB.n, 0) as n from (
                select distinct species_taxon_id,cluster_name from uniparc2species t1 
                inner join uniparc2mmseqs_90_80 t2 on t1.uniparc_id=t2.uniparc_id
                where (t1.species_taxon_id in ({taxid_filter}) or t1.genus_taxon_id in ({taxid_filter}))
            ) AA
            left join 
                (
                select taxon_id,cluster_name, count(*) as n from 
                    (select distinct t7.taxon_id,t5.accession,cluster_name from homology_search t1
                    inner join homology_search_db_xrefs t2 on t1.search_tool_id=t2.id
                    inner join homology_search_db_xrefs t3 on t1.substitution_matrix_id=t3.id
                    inner join homology_search_db_xrefs t4 on t1.filtering_id=t4.id
                    inner join genome_assembly_table t5 on t1.assembly_id=t5.assembly_id
                    inner join uniparc2mmseqs_90_80 t6 on t1.uniparc_id=t6.uniparc_id
                    inner join genome_assembly_table2taxon_id t7 on t5.assembly_id=t7.assembly_id
                    where t1.percent_identity>={min_identity}
                    and uniparc_seq_coverage>={min_coverage}
                    and t2.name in ("ssearch", "tfastx") 
                    and t3.name="BL50" 
                    and t4.name="no_filtering_lowcomplexity"
                    and t7.taxon_id in ({taxid_filter})) B group by taxon_id,cluster_name
                ) BB
            on AA.cluster_name=BB.cluster_name and AA.species_taxon_id=BB.taxon_id;
        '''

        df = pandas.read_sql(sql, self.conn)
        print("DF shape", df.shape)
        # calculate frequency
        def calculate(s):
            taxid = str(s["species_taxon_id"])
            n_genomes = taxid2n_genomes[taxid]
            s["freq"] = round((float(s["n"])/n_genomes) * 100, 2)
            return s

        df_with_freq = df.apply(calculate, axis=1)
        print("DF shape freq", df_with_freq.shape)
        return df_with_freq


  
    def get_VF_cluster2frequency_within_species(self,
                                             taxon_id, 
                                             min_identity=90, 
                                             min_coverage=60,
                                             percentages=False):
        '''
        Work with species and genus taxids
        '''
        
        
        if percentages:
            n_genomes = self.get_n_genomes(taxon_id)
        # get cluster list
        sql = f'''select distinct cluster_name from {CLUSTERING_TABLE} t1 
               inner join uniparc2species t2 on t1.uniparc_id=t2.uniparc_id where species_taxon_id={taxon_id}
               '''
        cluster_list = [i[0] for i in self.server.execute(sql,).fetchall()]

        # get cluster frequency
        # start from spacies table to count 
        # VF clusters with 0 hits on all genomes
        sql = f'''
            select AA.cluster_name,ifnull(BB.n, 0) from (
                select distinct cluster_name from uniparc2species t1 
                inner join uniparc2mmseqs_90_80 t2 on t1.uniparc_id=t2.uniparc_id
                where (t1.species_taxon_id={taxon_id} or t1.genus_taxon_id={taxon_id})
            ) AA
            left join 
                (
                select cluster_name, count(*) as n from 
                    (select distinct t5.accession,cluster_name from homology_search t1
                    inner join homology_search_db_xrefs t2 on t1.search_tool_id=t2.id
                    inner join homology_search_db_xrefs t3 on t1.substitution_matrix_id=t3.id
                    inner join homology_search_db_xrefs t4 on t1.filtering_id=t4.id
                    inner join genome_assembly_table t5 on t1.assembly_id=t5.assembly_id
                    inner join uniparc2mmseqs_90_80 t6 on t1.uniparc_id=t6.uniparc_id
                    inner join genome_assembly_table2taxon_id t7 on t5.assembly_id=t7.assembly_id
                    where t1.percent_identity>={min_identity}
                    and uniparc_seq_coverage>={min_coverage}
                    and t2.name in ("ssearch", "tfastx") 
                    and t3.name="BL50" 
                    and t4.name="no_filtering_lowcomplexity"
                    and t7.taxon_id={taxon_id}) B group by cluster_name
                ) BB
            on AA.cluster_name=BB.cluster_name;
        '''
        if not percentages:
            hsh = {i[0]:i[1] for i in self.server.execute(sql,).fetchall()}
            for i in cluster_list:
                if i not in hsh:
                    hsh[i] = 0
        else:
            hsh = {i[0]:(float(i[1])/n_genomes)*100 for i in self.server.execute(sql,).fetchall()}
            for i in cluster_list:
                if i not in hsh:
                    hsh[i] = 0
        return hsh
    
    def get_n_genomes(self, taxon_id=False):
        if not taxon_id:
            sql = f'''select taxon_id,count(*) as n from genome_assembly_table t1
            inner join genome_assembly_table2taxon_id t2 on t1.assembly_id=t2.assembly_id
            group by taxon_id
            '''
            return {str(i[0]):i[1] for i in self.server.execute(sql,).fetchall()} 
        else:
            sql = f'''select count(*) from genome_assembly_table t1
            inner join genome_assembly_table2taxon_id t2 on t1.assembly_id=t2.assembly_id
            where taxon_id={taxon_id}'''
            n_genomes = int(self.server.execute(sql,).fetchall()[0][0])
            return n_genomes
    
    def get_assembly_id(self, assembly_accession):
        sql = f'''select assembly_id from genome_assembly_table
        where accession="{assembly_accession}"'''
        print(sql)
        return self.server.execute(sql,).fetchall()[0][0]
 
    def get_assembly_accession2assembly_id(self, taxon_id):
        sql = f'''select accession,t1.assembly_id from genome_assembly_table t1
        inner join genome_assembly_table2taxon_id t2 on t1.assembly_id=t2.assembly_id
        where taxon_id={taxon_id}'''
        return {i[0]:i[1] for i in self.server.execute(sql,).fetchall()}
    
    def get_VF_cluster2frequency_within_species_separated_dbs(self,
                                                            taxon_id, 
                                                            min_identity=90, 
                                                            min_coverage=80,
                                                            percentages=False):
        
        if percentages:
            n_genomes = self.get_n_genomes(taxon_id)
               
        # get cluster frequency
        # start from spcies table to count 
        # VF clusters with 0 hits on all genomes
        sql = f'''
                select distinct EE.db_name,AA.cluster_name,ifnull(BB.n, 0) from (
                select distinct cluster_name,cluster_id from uniparc2species t1 
                inner join uniparc2mmseqs_90_80 t2 on t1.uniparc_id=t2.uniparc_id
                where t1.species_taxon_id={taxon_id}
                ) AA
                left join 
                (
                select cluster_name,cluster_id, count(*) as n from 
                (select distinct t5.accession,cluster_name,cluster_id from homology_search t1
                inner join homology_search_db_xrefs t2 on t1.search_tool_id=t2.id
                inner join homology_search_db_xrefs t3 on t1.substitution_matrix_id=t3.id
                inner join homology_search_db_xrefs t4 on t1.filtering_id=t4.id
                inner join genome_assembly_table t5 on t1.assembly_id=t5.assembly_id
                inner join uniparc2mmseqs_90_80 t6 on t1.uniparc_id=t6.uniparc_id
                inner join genome_assembly_table2taxon_id t8 on t5.assembly_id=t8.assembly_id
                where t1.percent_identity>={min_identity}
                and uniparc_seq_coverage>={min_coverage}
                and t2.name in ("ssearch", "tfastx") 
                and t3.name="BL50" 
                and t4.name="no_filtering_lowcomplexity"
                and t8.taxon_id={taxon_id}) B group by cluster_name,cluster_id
                ) BB
                on AA.cluster_id=BB.cluster_id
                inner join uniparc2mmseqs_90_80 CC on AA.cluster_id=CC.cluster_id 
                inner join VF_table DD on CC.uniparc_id = DD.uniparc_id 
                inner join VF_databases EE on DD.db_id =EE.db_id
                ;
        '''
        
        db2counts = {}
        for row in self.server.execute(sql,).fetchall():
            db = row[0]
            if percentages:
                val = (int(row[2])/n_genomes)*100
            else:
                val = int(row[2])
            if db not in db2counts:
                db2counts[db] = [val]
            else:
                db2counts[db].append(val)
        
        return db2counts

    def get_genome_assembly_id2_n_proteins(self, taxon_id):
        
        sql = f'''
        select t1.assembly_id,count(*) from genome_assembly_table t1 
        inner join genome_assembly_table2taxon_id t2 on t1.assembly_id=t2.assembly_id
        inner join circos_genome2proteins t3 on t1.assembly_id=t3.assembly_id 
        and t2.taxon_id={taxon_id} group by t1.assembly_id;
        '''
        
        return {i[0]:i[1] for i in self.server.execute(sql,).fetchall()}

    def get_genome_accession2description(self,
                                         taxon_id):
        
        sql = f'''
        select accession,description from genome_assembly_table t1
        inner join genome_assembly_table2taxon_id t2 on t1.assembly_id=t2.assembly_id
        where taxon_id={taxon_id}
        '''
        
        return {i[0]:i[1] for i in self.server.execute(sql,).fetchall()}
    
    def get_genome_accession2n_frameshifts(self,
                                           taxon_id,
                                           min_identity):
        
        sql = f'''select accession,count(*) as n from 
        (select t2.accession from homology_search t1 
        inner join genome_assembly_table t2 on t1.assembly_id=t2.assembly_id
        inner join genome_assembly_table2taxon_id t3 on t2.assembly_id=t3.assembly_id
        where (alignment_id is not NULL 
        and percent_identity>={min_identity} 
        and t3.taxon_id={taxon_id})) A 
        group by accession order by n DESC;
        '''
        
        return {i[0]:i[1] for i in self.server.execute(sql,).fetchall()}
           
        
    def get_genome_accession2statistics(self,
                                         taxon_id):
        
        sql = f'''select accession, GC, cumulated_size, n_pseudo from genome_assembly_table t1
        inner join genome_assembly_table2taxon_id t2 on t1.assembly_id=t2.assembly_id
        where taxon_id={taxon_id};
        '''
        genome_data = self.server.execute(sql,).fetchall()
        accession2gc = {i[0]:i[1] for i in genome_data}
        accession2genome_size = {i[0]:i[2] for i in genome_data}
        accession2pseudo = {i[0]:i[3] for i in genome_data}
        
        return accession2gc, accession2genome_size, accession2pseudo
        
        
    def get_genome_accession2mlst(self,
                                  taxon_id):
        
        sql = f'''select accession,st from mlst t1 
             inner join genome_assembly_table t2 on t1.assembly_id=t2.assembly_id
             inner join genome_assembly_table2taxon_id t3 on t2.assembly_id=t3.assembly_id
             where taxon_id={taxon_id};
          '''
        
        return {i[0]:i[1] for i in self.server.execute(sql,).fetchall()}
    
    def get_genome_accession2ani(self,
                                 taxon_id,
                                 reference_assembly_accession):
        
        sql = f'''select t3.accession,ani from ani t1 
        inner join genome_assembly_table t2 on t1.ref_assembly_id=t2.assembly_id 
        inner join genome_assembly_table t3 on t1.query_assembly_id=t3.assembly_id
        inner join genome_assembly_table2taxon_id t4 on t1.ref_assembly_id=t4.assembly_id
        inner join genome_assembly_table2taxon_id t5 on t1.query_assembly_id=t5.assembly_id
        where t1.taxon_id={taxon_id} 
        and t2.accession="{reference_assembly_accession}" 
        and t4.taxon_id={taxon_id}
        and t5.taxon_id={taxon_id};
          '''
        #print(sql)
        return {i[0]:i[1] for i in self.server.execute(sql,).fetchall()}


    def get_genome_n_VF_clusters(self,
                                 taxon_id, 
                                 min_identity=90, 
                                 min_coverage=80):
        '''
        Get count of number of VF clusters with significant hit(s) in each genome
        No distinct is made if a single VF cluster has one or multiple hits
        '''
        
        # get cluster frequency
        sql = f'''select A.accession,ifnull(n_VFs, 0) from (
            select t1.accession,t1.assembly_id from genome_assembly_table t1 
            inner join genome_assembly_table2taxon_id t2 on t1.assembly_id=t2.assembly_id
            where t2.taxon_id={taxon_id}
            ) A
            left join 
            (
            select B.assembly_id,B.accession,count(*) as n_VFs from (            
            select distinct t5.assembly_id,t5.accession,t5.description,hit_accession from homology_search t1
            inner join homology_search_db_xrefs t2 on t1.search_tool_id=t2.id
            inner join homology_search_db_xrefs t3 on t1.substitution_matrix_id=t3.id
            inner join homology_search_db_xrefs t4 on t1.filtering_id=t4.id
            inner join genome_assembly_table t5 on t1.assembly_id=t5.assembly_id
            inner join uniparc_entry t6 on t1.uniparc_id=t6.uniparc_id
            inner join genome_assembly_table2taxon_id t7 on t5.assembly_id=t7.assembly_id
            where t7.taxon_id={taxon_id}
            and t1.percent_identity>={min_identity}
            and uniparc_seq_coverage>={min_coverage}
            and t2.name in ("ssearch", "tfastx")
            and t3.name="BL50" 
            and t4.name="no_filtering_lowcomplexity"
            group by t5.accession,t1.hit_accession
            ) B group by B.assembly_id,B.accession) C
            on A.assembly_id=C.assembly_id
            ;'''
        #print(sql)
        return {i[0]:i[1] for i in self.server.execute(sql,).fetchall()}
    

    def get_db_cluster_accession_list(self,
                                      db_name):
        
        sql = f'''select distinct cluster_name, cluster_id from {CLUSTERING_TABLE} t1
                  inner join VF_table t2 on t1.uniparc_id=t2.uniparc_id
                  inner join VF_databases as t3 on t2.db_id=t3.db_id
                  inner join uniparc2species t4 on t1.uniparc_id=t4.uniparc_id
                  where t3.db_name="{db_name}" and t4.superkingdom_taxon_id=2 
              '''
        
        return [i[0] for i in self.server.execute(sql,).fetchall()]



    def get_species_cluster_accession_list(self,
                                           taxon_id, 
                                           db_name,
                                           rank='species'):
        
        sql = f'''
        select distinct cluster_name from (
        select distinct cluster_id from uniparc2species t1
        inner join uniparc2mmseqs_90_80 t2 on t1.uniparc_id=t2.uniparc_id 
        where t1.{rank}_taxon_id = {taxon_id}) A 
        inner join uniparc2mmseqs_90_80 B on A.cluster_id=B.cluster_id 
        inner join VF_table C on B.uniparc_id=C.uniparc_id 
        inner join VF_databases D on C.db_id= D.db_id 
        where D.db_name = "{db_name}"
        '''
        
        return [i[0] for i in self.server.execute(sql,).fetchall()]


    def get_db_cluster_id_list(self,
                               taxon_id, 
                               db_name,
                               rank='species'):
        
        sql = f'''
        select distinct A.cluster_id from (
        select distinct cluster_id from uniparc2species t1
        inner join uniparc2mmseqs_90_80 t2 on t1.uniparc_id=t2.uniparc_id 
        where t1.{rank}_taxon_id = {taxon_id}) A 
        inner join uniparc2mmseqs_90_80 B on A.cluster_id=B.cluster_id 
        inner join VF_table C on B.uniparc_id=C.uniparc_id 
        inner join VF_databases D on C.db_id= D.db_id 
        where D.db_name = "{db_name}"
        '''
        
        return [i[0] for i in self.server.execute(sql,).fetchall()]
   
   
    def get_VF_mmseqs_clusters(self, taxon_id):
        
        '''
        join on protein_accession?
                
        '''
        
        sql = ''
        
        return pandas.read_sql(sql, self.conn)
        
    def pan_genome_clusters(self,taxon_id):
        

        sql = f'''select distinct t3.accession,t2.cluster_id,count(*) as n from circos_genome2proteins t1 
        inner join circos_proteins2mmseqs t2 on t1.protein_id=t2.protein_id 
        inner join genome_assembly_table t3 on t1.assembly_id=t3.assembly_id
        inner join genome_assembly_table2taxon_id t4 on t3.assembly_id=t4.assembly_id
        where t4.taxon_id={taxon_id}
        and t2.id_cutoff=80
        and t2.cov_cutoff=80
        group by t1.assembly_id,t2.cluster_id;
        '''
    
        return pandas.read_sql(sql, self.conn)
   
   
    def get_genome_accession2ani(self,
                                 taxon_id):
        
        sql = f'''select t2.accession,t3.accession,ani from ani t1 
        inner join genome_assembly_table t2 on t1.ref_assembly_id=t2.assembly_id 
        inner join genome_assembly_table t3 on t1.query_assembly_id=t3.assembly_id
        inner join genome_assembly_table2taxon_id t4 on t1.ref_assembly_id=t4.assembly_id
        inner join genome_assembly_table2taxon_id t5 on t1.query_assembly_id=t5.assembly_id
        where t1.taxon_id={taxon_id} 
        and t4.taxon_id={taxon_id}
        and t5.taxon_id={taxon_id};
          '''
        #print(sql)
        return pandas.read_sql(sql, self.conn)
        
        
    def VF_cluster_genome_counts(self,
                                 taxon_id, 
                                 min_identity=90, 
                                 min_coverage=60,
                                 search_tools=["ssearch","tfastx"],
                                 matrices=["BL50"]):
        '''
        Get number of hits for each cluster in each genome assembly
        Only possible for ssearch and not tfastx (would need to check position)
        '''
        
        tool_filter = '","'.join(search_tools)
        matrices_filter = '","'.join(matrices)
        # get cluster frequency
        sql_VF_freq = f'''select t6.cluster_name, t5.accession, count(*) as n from homology_search t1 
                    inner join homology_search_db_xrefs t2 on t1.search_tool_id=t2.id 
                    inner join homology_search_db_xrefs t3 on t1.substitution_matrix_id=t3.id 
                    inner join homology_search_db_xrefs t4 on t1.filtering_id=t4.id 
                    inner join genome_assembly_table t5 on t1.assembly_id=t5.assembly_id 
                    inner join {CLUSTERING_TABLE} t6 on t1.uniparc_id=t6.uniparc_id 
                    inner join {CLUSTERING2SPECIES_TABLE} t7 on t6.cluster_id=t7.cluster_id
                    inner join genome_assembly_table2taxon_id t8 on t5.assembly_id=t8.assembly_id
                    where t8.taxon_id={taxon_id}
                    and t7.species_taxon_id={taxon_id}
                    and percent_identity>={min_identity}
                    and uniparc_seq_coverage>={min_coverage}
                    and t2.name in ("{tool_filter}") 
                    and t3.name in ("{matrices_filter}")
                    and t4.name="no_filtering_lowcomplexity"
                    group by t1.assembly_id,t6.cluster_id;'''
        
        return self.server.execute(sql_VF_freq,).fetchall()
        

    def VF_cluster_genome_2mmseqs(self,
                                  taxon_id, 
                                  min_identity=90, 
                                  min_coverage=60,
                                  search_tools=["ssearch"],
                                  matrices=["BL50"]):
        '''
        Get number of hits for each cluster in each genome assembly
        Only possible for ssearch and not tfastx (would need to check position)
        '''
        
        tool_filter = '","'.join(search_tools)
        matrices_filter = '","'.join(matrices)
        # get cluster frequency
        sql_VF_freq = f'''select distinct t5.accession,t10.cluster_id  as n from homology_search t1 
                    inner join homology_search_db_xrefs t2 on t1.search_tool_id=t2.id 
                    inner join homology_search_db_xrefs t3 on t1.substitution_matrix_id=t3.id 
                    inner join homology_search_db_xrefs t4 on t1.filtering_id=t4.id 
                    inner join genome_assembly_table t5 on t1.assembly_id=t5.assembly_id 
                    inner join {CLUSTERING_TABLE} t6 on t1.uniparc_id=t6.uniparc_id 
                    inner join {CLUSTERING2SPECIES_TABLE} t7 on t6.cluster_id=t7.cluster_id
                    inner join genome_assembly_table2taxon_id t8 on t5.assembly_id=t8.assembly_id
                    inner join circos_genome2proteins t9 on t1.hit_accession=t9.accession
                    inner join circos_proteins2mmseqs t10 on t9.protein_id=t10.protein_id 
                    where t8.taxon_id={taxon_id}
                    and t10.taxon_id={taxon_id}
                    and t7.species_taxon_id={taxon_id}
                    and percent_identity>={min_identity}
                    and uniparc_seq_coverage>={min_coverage}
                    and t2.name in ("{tool_filter}") 
                    and t3.name in ("{matrices_filter}")
                    and t4.name="no_filtering_lowcomplexity"'''
        
        print(sql)
        
        return self.server.execute(sql_VF_freq,).fetchall()
        


    def get_species_name_from_taxid(self,taxon_id):
        
        sql = f'select species from ncbi_taxonomy where tax_id={taxon_id};'
        
        return self.server.execute(sql,).fetchall()[0][0]

    def get_species_VF_UP_count(self,
                                taxon_id,
                                rank="species"):
        
        sql = f'''select count(*) from (select distinct t1.uniparc_id from {CLUSTERING_TABLE} t1 
        inner join uniparc2species t2 on t1.uniparc_id=t2.uniparc_id where {rank}_taxon_id={taxon_id}) A;'''
        
        return self.server.execute(sql,).fetchall()[0][0]
       
        
    def get_combined_dbs_cluster_count(self, rank):


        sql2 = f'''
            select A.{rank}_name, count(*) as n_cluster from
            (select t1.{rank}_name,{rank}_taxon_id,t2.cluster_id from uniparc2species t1
            inner join {CLUSTERING_TABLE} t2 on t1.uniparc_id =t2.uniparc_id 
            where t1.{rank}_name not like "%%%% sp. %%%%" 
            and t1.{rank}_name not like "%%%%uncultured%%%%" 
            and t1.{rank}_name not like "%%%%genomosp%%%%" 
            and t1.superkingdom_taxon_id=2
            group by t1.{rank}_name, t2.cluster_id) A group by A.{rank}_name,A.{rank}_taxon_id;
            '''

        rank2count_all = {i[0]:i[1] for i in self.server.execute(sql2,).fetchall()}
        
        return rank2count_all

    def get_cluster_id(self, cluster_name_list):
        
        cluster_filter = '","'.join(cluster_name_list)
        
        sql = f'select cluster_name,cluster_id from {CLUSTERING_TABLE} where cluster_name in ("{cluster_filter}")'
        
        return {i[0]:i[1] for i in self.server.execute(sql,).fetchall()}

    def get_uniparc_id(self, uniparc_accession_list):
        
        up_filter = '","'.join(uniparc_accession_list)
        
        sql = f'select uniparc_accession,uniparc_id from uniparc_entry where uniparc_accession in ("{up_filter}")'
        
        return {i[0]:i[1] for i in self.server.execute(sql,).fetchall()}

    
    def get_clusters_from_gene_names(self,
                                     rank_name, 
                                     taxon_id,
                                     gene_names_list):
        
        gene_list = [i[0] for i in gene_names_list if i[0] is not None]
        gene_filter = '","'.join(gene_list)
        sql_gene2cluster = f'select distinct cluster_name,gene from uniparc_consensus_annotation t1 ' \
                         f' inner join uniparc2species t2 on t1.uniparc_id=t2.uniparc_id ' \
                         f' inner join {CLUSTERING_TABLE} t3 on t1.uniparc_id=t3.uniparc_id ' \
                         f' inner join VF_table t4 on t1.uniparc_id=t4.uniparc_id ' \
                         f'where {rank_name}_taxon_id={taxon_id} and gene in ("{gene_filter}");'
        
        return self.server.execute(sql_gene2cluster,).fetchall()


    def get_db_2pmid_list(self,):
        
        sql = 'select distinct db_name,pmid from VF_id2pmid t1 inner join VF_table t2 on t1.VF_id=t2.VF_id inner join VF_databases t3 on t2.db_id=t3.db_id;'
        data = self.server.execute(sql,).fetchall()
        db2pmid_list = {}
        for row in data:
            if row[0] not in db2pmid_list:
                db2pmid_list[row[0]] = []
            db2pmid_list[row[0]].append(str(row[1]))
        
        return db2pmid_list
 
    def get_pmid2description(self,):
        
        sql = 'select t1.pmid,title from VF_id2pmid t1 inner join pmid2data t2 on t1.pmid=t2.pmid group by t1.pmid,title;'
        
        pmi2descr = {}
        for row in self.server.execute(sql,).fetchall():
            title_format = row[1].replace("'","").replace('"',"")
            pmi2descr[str(row[0])] = f"<td>{title_format}</td>"

        # 
        return pmi2descr
 
    def get_mmseq_cluster_frequency(self, 
                                    taxon_id,
                                    percentage=False):
        
        sql = f'''
        select cluster_id,count(*) as n from (select distinct t1.assembly_id,cluster_id from circos_genome2proteins t1 
        inner join circos_proteins2mmseqs t2 on t1.protein_id=t2.protein_id 
        inner join genome_assembly_table t3 on t1.assembly_id=t3.assembly_id
        inner join genome_assembly_table2taxon_id t4 on t3.assembly_id=t4.assembly_id
        where t4.taxon_id={taxon_id}
        and t2.id_cutoff=80
        and t2.cov_cutoff=80
        ) A group by A.cluster_id;
        '''
        
        cluster2count = {i[0]:i[1] for i in self.server.execute(sql,).fetchall()}
 
        if percentage:
            n_genomes = self.get_n_genomes(taxon_id)
            cluster2count_percent = {}
            for cluster, count in cluster2count.items():
                cluster2count_percent[cluster] = (count/float(n_genomes))*100
            return cluster2count_percent
        else:
            return cluster2count

    def get_VF_cluster_hits(self, 
                            cluster_id,
                            taxon_id,
                            search_tool='ssearch',
                            matrix='BL50',
                            filtering='no_filtering_lowcomplexity',
                            identity_cutoff=90,
                            coverage_cutoff=60):
        
        # TODO: add filter on VF taxnonomy
        sql = f'''select t2.accession as assembly_accession,t1.hit_accession, t5.uniparc_accession,
                    t1.percent_identity,uniparc_seq_coverage,assembly_seq_coverage, 
                    t1.evalue, t1.bitscore,t6.description
                    from homology_search t1
                    inner join genome_assembly_table t2 on t1.assembly_id=t2.assembly_id
                    inner join {CLUSTERING_TABLE} t4 on t1.uniparc_id=t4.uniparc_id
                    inner join uniparc_entry t5 on t1.uniparc_id=t5.uniparc_id
                    inner join uniparc_consensus_annotation t6 on t1.uniparc_id=t6.uniparc_id
                    inner join homology_search_db_xrefs t7 on t1.search_tool_id=t7.id
                    inner join homology_search_db_xrefs t8 on t1.substitution_matrix_id=t8.id
                    inner join homology_search_db_xrefs t9 on t1.filtering_id=t9.id
                    inner join genome_assembly_table2taxon_id t10 on t2.assembly_id=t10.assembly_id
                    inner join mmseqs_90_80_2species t11 on t4.cluster_id=t11.cluster_id 
                    where t4.cluster_id="{cluster_id}"
                    and percent_identity>={identity_cutoff}
                    and uniparc_seq_coverage>={coverage_cutoff} 
                    and t7.name="{search_tool}" 
                    and t8.name="{matrix}" 
                    and t9.name="{filtering}"
                    and t10.taxon_id={taxon_id}
                    and (t11.species_taxon_id={taxon_id} or t11.genus_taxon_id={taxon_id})
                    '''

        return pandas.read_sql(sql, self.conn)

    def get_VF_cluster_best_hits(self, 
                                 cluster_id,
                                 taxon_id,
                                 search_tool='ssearch',
                                 matrix='BL50',
                                 filtering='no_filtering_lowcomplexity',
                                 identity_cutoff=90,
                                 coverage_cutoff=60):

        all_hits_df = self.get_VF_cluster_hits(cluster_id=cluster_id,
                                            taxon_id=taxon_id,
                                            search_tool=search_tool,
                                            matrix=matrix,
                                            filtering=filtering,
                                            identity_cutoff=identity_cutoff,
                                            coverage_cutoff=identity_cutoff)
        
        print("all_hits_df", all_hits_df.head())
        
        assembly_accession2best_hit = {}
        for n, row in all_hits_df.iterrows():
            if row["assembly_accession"] not in assembly_accession2best_hit:
                assembly_accession2best_hit[row["assembly_accession"]] = {}
                assembly_accession2best_hit[row["assembly_accession"]]["hit_accession"] = row["hit_accession"]
                assembly_accession2best_hit[row["assembly_accession"]]["uniparc_accession"] = row["uniparc_accession"]
                assembly_accession2best_hit[row["assembly_accession"]]["percent_identity"] = row["percent_identity"]
                assembly_accession2best_hit[row["assembly_accession"]]["uniparc_seq_coverage"] = row["uniparc_seq_coverage"]
                assembly_accession2best_hit[row["assembly_accession"]]["assembly_seq_coverage"] = row["assembly_seq_coverage"]
                assembly_accession2best_hit[row["assembly_accession"]]["evalue"] = row["evalue"]
                assembly_accession2best_hit[row["assembly_accession"]]["bitscore"] = row["bitscore"]
                assembly_accession2best_hit[row["assembly_accession"]]["description"] = row["description"]

        return assembly_accession2best_hit
    
    def get_assembly_VF_list(self, assembly_accession, identity_cutoff, taxon_id):
        
        # TODO: add filter on VF taxnonomy
        sql = f'''select t1.hit_accession, t4.cluster_name, t5.uniparc_accession,
                    t1.percent_identity,uniparc_seq_coverage,assembly_seq_coverage, 
                    t1.evalue, t1.bitscore,t6.description
                    from homology_search t1
                    inner join genome_assembly_table t2 on t1.assembly_id=t2.assembly_id
                    inner join {CLUSTERING_TABLE} t4 on t1.uniparc_id=t4.uniparc_id
                    inner join uniparc_entry t5 on t1.uniparc_id=t5.uniparc_id
                    inner join uniparc_consensus_annotation t6 on t1.uniparc_id=t6.uniparc_id
                    inner join homology_search_db_xrefs t7 on t1.search_tool_id=t7.id
                    inner join homology_search_db_xrefs t8 on t1.substitution_matrix_id=t8.id
                    inner join homology_search_db_xrefs t9 on t1.filtering_id=t9.id
                    inner join genome_assembly_table2taxon_id t10 on t2.assembly_id=t10.assembly_id
                    inner join mmseqs_90_80_2species t11 on t4.cluster_id=t11.cluster_id 
                    where accession="{assembly_accession}" 
                    and percent_identity>={identity_cutoff}
                    and uniparc_seq_coverage>=60 
                    and t7.name="ssearch" 
                    and t8.name="BL50" 
                    and t9.name="no_filtering_lowcomplexity"
                    and t10.taxon_id={taxon_id}
                    and (t11.species_taxon_id={taxon_id} or t11.genus_taxon_id={taxon_id})
                    '''
        #print(sql)           
        return pandas.read_sql(sql, self.conn)
    
    
    def get_proportion_core_genes(self,):
        
        sql = '''
        select A.taxon_id,(CAST(A.n_core as FLOAT)/B.proteins)*100 as percentage from 
        (select taxon_id,assembly_id,value as n_core from assembly_statistics t1 
        inner join terms t2 on t1.term_id=t2.term_id where term_name="n_core_90") A 
        inner join (select assembly_id,value as proteins from assembly_statistics t1 
        inner join terms t2 on t1.term_id=t2.term_id where term_name="n_proteins") B 
        on A.assembly_id=B.assembly_id;
        '''
                    
        return pandas.read_sql(sql, self.conn)    

    def get_VFs_assembly_stats(self,percentages=False):
        
        sql = '''
        select taxon_id,value as n_core_90 from assembly_statistics t1 
        inner join terms t2 on t1.term_id=t2.term_id where term_name="n_VF_core_90";
        '''
        species_median_N_core = pandas.read_sql(sql, self.conn)
        #print(species_median_N_core.head())
        species_median_N_core["n_core_90"] = pandas.to_numeric(species_median_N_core["n_core_90"]) 
        species_median_N_core["taxon_id"] = species_median_N_core["taxon_id"].astype(str) 
        taxid2core_VF_median = species_median_N_core.groupby(["taxon_id"]).mean().to_dict()["n_core_90"]

        sql = '''
        select taxon_id,value as n_pan_90 from assembly_statistics t1 
        inner join terms t2 on t1.term_id=t2.term_id where term_name="n_VF_pan_90";
        '''
        species_median_N_pan = pandas.read_sql(sql, self.conn)
        species_median_N_pan["n_pan_90"] = pandas.to_numeric(species_median_N_pan["n_pan_90"])
        species_median_N_pan["taxon_id"] = species_median_N_pan["taxon_id"].astype(str) 
        taxid2pan_VF_median = species_median_N_pan.groupby(["taxon_id"]).mean().to_dict()["n_pan_90"]
        
        taxid_list = set(list(taxid2pan_VF_median.keys()) + list(taxid2core_VF_median.keys()))
        for taxid in taxid_list:
            if taxid not in taxid2core_VF_median:
                taxid2core_VF_median[taxid] = 0
            if taxid not in taxid2pan_VF_median:
                taxid2pan_VF_median[taxid] = 0 
                
        if not percentages:
            return taxid2core_VF_median, taxid2pan_VF_median
        else:
            print("taxon-id", taxid)
            taxid2proportion_core = {}
            taxid2proportion_pan = {} 
            for taxid in taxid_list:
                combined = (taxid2core_VF_median[taxid] + taxid2pan_VF_median[taxid])
                print("taxid",taxid,"n combined", combined)
                print("core", taxid2core_VF_median[taxid])
                print("pan", taxid2pan_VF_median[taxid])
                taxid2proportion_core[taxid] = round((taxid2core_VF_median[taxid] / combined) * 100, 2)
                taxid2proportion_pan[taxid] = round((taxid2pan_VF_median[taxid] / combined) * 100, 2)
            return taxid2proportion_core, taxid2proportion_pan

    def get_assemblies_size(self,):
        
        sql = '''
        select taxon_id,CAST(cumulated_size as FLOAT)/1000000 as genome_size from genome_assembly_table t1 
        inner join genome_assembly_table2taxon_id t2 on t1.assembly_id=t2.assembly_id;
        '''
                    
        return pandas.read_sql(sql, self.conn)    
    
    
    
    def write_circos_files(self, assembly_accession, taxon_id, output_path):
        
        import metagenlab_libs.gbk2circos as gbk2circos
        
        assembly_acc_no_version = assembly_accession.split(".")[0]
        
        genes_minus = open(os.path.join(output_path, f"circos_{assembly_acc_no_version}_minus.txt"), "w")
        genes_plus = open(os.path.join(output_path, f"circos_{assembly_acc_no_version}_plus.txt"), "w")
        conservation = open(os.path.join(output_path, f"circos_{assembly_acc_no_version}_conservation.txt"), "w")
        caryotype = open(os.path.join(output_path, f"circos_{assembly_acc_no_version}_caryotype.txt"), "w")
        config = open(os.path.join(output_path, f"circos_{assembly_acc_no_version}.config"), "w")
        rare_VFs = open(os.path.join(output_path, f"circos_{assembly_acc_no_version}_rare_VFs.txt"), "w")
        
        genome_data = self.get_genome_data(assembly_accession, taxon_id)
        cluster_id2frequency = self.get_mmseq_cluster_frequency(taxon_id)
        VF_cluster2freq = self.get_VF_cluster2frequency_within_species(taxon_id, 
                                                                       min_identity=80, 
                                                                       min_coverage=60,
                                                                       percentages=True)
        
        n_genomes = self.get_n_genomes(taxon_id)
        
        record_list = genome_data["record_id"].unique()
        spacing_pairs = []
        for n, record in enumerate(record_list):
            if n > 0:
                spacing_pairs.append([record_list[n], record_list[n-1]])
            if n % 2 == 0:
                col = 'spectral-5-div-3'
            else:
                col = 'spectral-5-div-4'
            start = 0
            stop = max(genome_data["end"][genome_data["record_id"] == record])
            caryotype.write(f"chr - {record} {record} 0 {stop} {col}\n")
        if len(record_list) > 2:
            spacing_pairs.append([record_list[0], record_list[-1]])

        print("chr_spacing_list", spacing_pairs)
        circos_conf = gbk2circos.Circos_config(caryotype_file=caryotype.name.split("/")[-1],
                                                radius=0.7,
                                                show_tick_labels="yes",
                                                chr_spacing_list=spacing_pairs)
        
        circos_conf.add_highlight(genes_plus.name.split("/")[-1], fill_color="grey_a1", r1="0.88r", r0="0.84r")
        circos_conf.add_highlight(genes_minus.name.split("/")[-1], fill_color="grey_a1", r1="0.84r", r0="0.80r")
        supp = f'''

<rules>
<rule>
condition          = var(value) < 90
fill_color         = lred
</rule>

<rule>
condition          = var(value) > 90
fill_color         = lgreen
</rule>

</rules>
'''

        circos_conf.add_plot(conservation.name.split("/")[-1], type="line", r0="0.89r", r1="0.99r",color="black", rules=supp,
                            thickness="0.5p")


        supp = '''
                label_snuggle             = yes
                max_snuggle_distance            = 20r
                show_links     = yes
                link_dims      = 10p,168p,30p,4p,4p
                link_thickness = 2p
                link_color     = blue
                label_size   = 24p
                label_font   = condensed
                padding  = 0p
                rpadding = 0p
                '''

        circos_conf.add_plot(rare_VFs.name.split("/")[-1] ,type="text", r0="1r", r1="1.3r",color="black", rules=supp)



        config.write(circos_conf.get_file())
        config.close()

        VF_detail = self.get_assembly_VF_list(assembly_accession, 90, taxon_id)
        cluster_filter = '","'.join(VF_detail["cluster_name"].unique())
        sql = f'''select distinct cluster_name,gene from uniparc_consensus_annotation t1 
        inner join {CLUSTERING_TABLE} t2 on t1.uniparc_id=t2.uniparc_id where t2.cluster_name in ("{cluster_filter}")'''
        #print(sql)
        cluster2gene = {i[0]:i[1] for i in self.server.execute(sql,).fetchall()}
        
        VF_nr_list = VF_detail["hit_accession"].unique()
        
        for n, row in genome_data.iterrows():
            if abs(row["end"] - row["start"]) > 10000:
                continue
            if row["accession"] in VF_nr_list:
                # get first cluster name
                vf_cluster = VF_detail[VF_detail["hit_accession"]==row["accession"]]["cluster_name"].unique()[0]
                try:
                    vf_cluster_freq = VF_cluster2freq[vf_cluster]
                except:
                    vf_cluster_freq = 0
                if vf_cluster_freq < 80:
                    rare_VFs.write(f'{row["record_id"]}\t{row["start"]}\t{row["end"]}\t{cluster2gene[vf_cluster]}\n')
                col = 'fill_color=piyg-5-div-1,'
            else:
                col =''
            if row["strand"] == 1:
                genes_plus.write(f'{row["record_id"]}\t{row["start"]}\t{row["end"]}\t{col}id={row["accession"]}\n')
            if row["strand"] == -1:
                genes_minus.write(f'{row["record_id"]}\t{row["start"]}\t{row["end"]}\t{col}id={row["accession"]}\n')
            try:
                freq = (float(cluster_id2frequency[row["cluster_id"]])/n_genomes)*100
            except:
                freq = 0
            
            conservation.write(f'{row["record_id"]}\t{row["start"]}\t{row["end"]}\t{freq}\n')
            
        return config.name

    def get_genome_data(self, assembly_accession, taxon_id):
        
        sql = f'''select record_id,t1.accession,start,end,strand,plasmid,t3.cluster_id from circos_genome2proteins t1 
        inner join genome_assembly_table t2 on t1.assembly_id=t2.assembly_id
        inner join circos_proteins2mmseqs t3 on t1.protein_id=t3.protein_id
        inner join genome_assembly_table2taxon_id t4 on t2.assembly_id=t4.assembly_id
        where t2.accession="{assembly_accession}" 
        and t4.taxon_id={taxon_id}
        and t3.taxon_id={taxon_id}
        and t3.id_cutoff=80
        and t3.cov_cutoff=80
        order by record_id,start;
        '''
        print(sql)
        table = pandas.read_sql(sql, self.conn)
        #print(table.head())
 
        return table

    def get_rank(self, taxon_id):
        from ete3 import NCBITaxa
        ncbi = NCBITaxa()
        print(ncbi.get_rank([taxon_id]))
        rank = ncbi.get_rank([taxon_id])[int(taxon_id)]
        return rank

    def get_rank_name(self, rank, taxon_id):
    
        sql = f'select distinct {rank}_name from uniparc2species where {rank}_taxon_id={taxon_id}'
    
        return self.server.execute(sql,).fetchall()[0][0]

    def get_species_phylogeny(self, taxon_id):
    
        sql = f'select phylogeny from species_phylogeny where taxon_id={taxon_id}'
    
        return self.server.execute(sql,).fetchall()[0][0]

    def get_genome_distance(self, taxon_id, distance='AAI_median'):
        
        sql = f'''
            select t3.accession,t4.accession, distance from genome_assembly_distance t1 
            inner join terms t2 on t1.term_id=t2.term_id 
            inner join genome_assembly_table t3 on t1.assembly_b=t3.assembly_id 
            inner join genome_assembly_table t4 on t1.assembly_a=t4.assembly_id 
            inner join genome_assembly_table2taxon_id t5 on t3.assembly_id=t5.assembly_id
            inner join genome_assembly_table2taxon_id t6 on t4.assembly_id=t6.assembly_id
            where t5.taxon_id={taxon_id}
            and t6.taxon_id={taxon_id} 
            and term_name="{distance}"
        '''
        return pandas.read_sql(sql, self.conn)
        
    def get_accession2genome_distance(self, assembly_id, distance="AAI"):
        
        sql = f'''select accession,distance from genome_assembly_distance t1 
                 inner join terms t2 on t1.term_id=t2.term_id 
                 inner join genome_assembly_table t3 on t1.assembly_b=t3.assembly_id 
                 where assembly_a={assembly_id} and term_name="{distance}"
                 union 
                 select accession,distance from genome_assembly_distance t1 
                 inner join terms t2 on t1.term_id=t2.term_id 
                 inner join genome_assembly_table t3 on t1.assembly_a=t3.assembly_id 
                 where assembly_b={assembly_id} and term_name="{distance}"
                 ;
             '''
        
        return {i[0]:round(i[1],1) for i in self.server.execute(sql,).fetchall()}
        

    def get_ducplicate_gene_names(self,
                                  rank,
                                  taxon_id):
        
        
        sql = f'''select * from (select gene,count(*) as n from
                 (select distinct cluster_name,gene from uniparc_consensus_annotation t1
                 inner join uniparc2species t2 on t1.uniparc_id=t2.uniparc_id
                 inner join {CLUSTERING_TABLE} t3 on t1.uniparc_id=t3.uniparc_id
                 inner join VF_table t4 on t1.uniparc_id=t4.uniparc_id
                 where {rank}_taxon_id={taxon_id}) A group by gene) B where n >1;
               ''' # and db_id!=2

        return self.server.execute(sql,).fetchall()
        

    def cluster2annotation(self, 
                           cluster_id_list, 
                           accessions=False):
        
        if accessions:
            column = 'cluster_name'
            cluster_fam_filter = '"%s"' % '","'.join([str(i) for i in cluster_id_list])
        else:
            column = 'cluster_id'
            cluster_fam_filter = ','.join([str(i) for i in cluster_id_list])

            
        print("NUMBER OF CLUSTERS:", len (cluster_id_list))

        sql = f'select cluster_name,uniparc_accession,description from uniparc_entry t1' \
            f' inner join {CLUSTERING_TABLE} t2 on t1.uniparc_id=t2.uniparc_id ' \
            f' inner join uniparc_consensus_annotation t3 on t1.uniparc_id=t3.uniparc_id where t2.{column} in ({cluster_fam_filter});'

        sql2 = f'select distinct cluster_name,phylum_name,species_name from {CLUSTERING_TABLE} t1 ' \
            f' inner join uniparc2species t2 on t1.uniparc_id=t2.uniparc_id where t1.{column} in ({cluster_fam_filter}) and species_name not like "%% sp.%%";'

        sql3 = f'select cluster_name,count(*) as n_db from ' \
            f' (select distinct cluster_name,db_id from uniparc_entry t1 ' \
            f' inner join {CLUSTERING_TABLE} t2 on t1.uniparc_id=t2.uniparc_id ' \
            f' inner join VF_table t3 on t1.uniparc_id=t3.uniparc_id where t2.{column} in ({cluster_fam_filter})) A group by cluster_name'

        self.server.execute(sql,)
        annotation = self.server.fetchall()
        self.server.execute(sql2,)
        species = self.server.fetchall()
        self.server.execute(sql3,)
        cluster_name2n_db = {i[0]:i[1] for i in self.server.fetchall()}

        cluster2annotations = {}
        cluster2species = {}
        for row in species:
            if row[0] not in cluster2species:
                cluster2species[row[0]] = ["%s (%s)" % (row[2], row[1])]
            else:
                cluster2species[row[0]].append("%s (%s)" % (row[2], row[1]))

        for row in annotation:
            if row[0] not in cluster2annotations:
                # species (phylum)
                cluster2annotations[row[0]] = [row[1:]]
            else:
                cluster2annotations[row[0]].append(row[1:])
        return cluster2annotations, cluster2species, cluster_name2n_db


      
def get_cluster_profile_matrix(taxon_id, min_identity=95):

    # get assembly_accession 2 cluster
    # reformat to build a matrix of presence/absence
    # genome            VF_1    VF_2    vf_3
    # GCF_000009665     1       1       0
    # GCF_000011525     0       1       0
    #

    import os
    import sqlite3
    db = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + "/virulence.db"
    conn = sqlite3.connect(db)

    cursor = conn.cursor()

    cursor = conn.cursor()
    sql = '''select distinct accession,cluster_name from
           (select distinct t2.accession,t1.VF_id from genome_blast_results t1 
           inner join genome_assembly_table t2 on t1.assembly_id=t2.assembly_id
           inner join genome_assembly_table2taxon_id t3 on t2.assembly_id=t3.assembly_id
           where percent_identity>=%s and t3.taxon_id=%s) A 
           inner join cluster_90 B on A.VF_id=B.VF_id;' % (min_identity, taxon_id)
        '''
    cursor.execute(sql,)
    assembly_clusters = cursor.fetchall()
    cluster_list = list(set([i[1] for i in assembly_clusters]))
    assembly_list = list(set([i[0] for i in assembly_clusters]))

    assembly2cluster_count = {}
    cluster2count = {}
    for row in assembly_clusters:
        if row[0] not in assembly2cluster_count:
            assembly2cluster_count[row[0]] = {}
            assembly2cluster_count[row[0]][row[1]] = 1
        else:
            assembly2cluster_count[row[0]][row[1]] = 1
        # count clusters
        if row[1] not in cluster2count:
            cluster2count[row[1]] =1
        else:
            cluster2count[row[1]] +=1

    cluster_non_universal = []
    for cluster in cluster_list:
        #print cluster2count[cluster], len(assembly2cluster_count), (len(assembly2cluster_count)*0.95)
        if cluster2count[cluster] < (len(assembly2cluster_count)*1):
            cluster_non_universal.append(cluster)

    #print assembly2cluster_count

    for assembly in assembly_list:
        for cluster in cluster_non_universal:
            if cluster not in assembly2cluster_count[assembly]:
                assembly2cluster_count[assembly][cluster] = 0
    #print 'N clusters: %s' % len(cluster_list)
    #print 'N non universal: %s' % len(cluster_non_universal)



    #print '\t'.join(cluster_list)
    for assembly in assembly_list:
        data = [str(assembly2cluster_count[assembly][cluster]) for cluster in cluster_non_universal]

#get_cluster_profile_matrix(1280, min_identity=90)




def make_div(figure_or_data, include_plotlyjs=False, show_link=False, div_id=None):
    from plotly import offline
    div = offline.plot(
        figure_or_data,
        include_plotlyjs=include_plotlyjs,
        show_link=show_link,
        output_type="div",
    )
    if ".then(function ()" in div:
        div = """{div.partition(".then(function ()")[0]}</script>"""
    if div_id:
        import re

        try:
            existing_id = re.findall(r'id="(.*?)"|$', div)[0]
            #print(existing_id, div_id)
            div = div.replace(existing_id, div_id)
        except IndexError:
            pass
    return div




#cluster_fam_freq(95)
