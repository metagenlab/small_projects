#!/usr/bin/env python
from django.conf import settings
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
        self.db_name = settings.DB_NAME

    def load_db():

        db = settings.DB_NAME
        conn = sqlite3.connect(db)
        cursor = conn.cursor()
        
        return DB(conn, cursor)
    
    
    def get_VF_conservation_matrix(taxon_id, identity_cutoff=90, coverage_cutoff=80):
        
        sql = '''
        
        
        
        '''
    
        # convert long table to matrix with pandas
        
        
    
    
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
        
        return {i[0]:i[1] for i in self.server.execute(sql,).fetchall()}
 
        
    def get_cluster2frequency_within_species(self,
                                             taxon_id, 
                                             min_identity=90, 
                                             min_coverage=60,
                                             percentages=False):
        
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
                where t1.species_taxon_id={taxon_id}
            ) AA
            left join 
                (
                select cluster_name, count(*) as n from 
                    (select distinct t5.accession,cluster_name from homology_search t1
                    inner join homology_search_db_xrefs t2 on t1.search_tool_id=t2.id
                    inner join homology_search_db_xrefs t3 on t1.substitution_matrix_id=t3.id
                    inner join homology_search_db_xrefs t4 on t1.filtering_id=t4.id
                    inner join genome_assembly_table t5 on t1.assembly_id=t5.assembly_id
                    inner join uniparc2mmseqs_90_80 t3 on t1.uniparc_id=t3.uniparc_id
                    where t1.percent_identity>={min_identity}
                    and uniparc_seq_coverage>={min_coverage}
                    and t2.name in ("ssearch", "tfastx") 
                    and t3.name="BL50" 
                    and t4.name="no_filtering_lowcomplexity"
                    and t5.taxon_id={taxon_id}) B group by cluster_name
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
    
    def get_n_genomes(self, taxon_id):
        sql = f'select count(*) from genome_assembly_table where taxon_id={taxon_id}'
        n_genomes = int(self.server.execute(sql,).fetchall()[0][0])
        return n_genomes
    
    
    def get_cluster2frequency_within_species_separated_dbs(self,
                                                           taxon_id, 
                                                           min_identity=90, 
                                                           min_coverage=80,
                                                           percentages=False):
        
        if percentages:
            n_genomes = self.get_n_genomes(taxon_id)
               
        # get cluster frequency
        # start from spacies table to count 
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
                inner join uniparc2mmseqs_90_80 t3 on t1.uniparc_id=t3.uniparc_id
                where t1.percent_identity>={min_identity}
                and uniparc_seq_coverage>={min_coverage}
                and t2.name in ("ssearch", "tfastx") 
                and t3.name="BL50" 
                and t4.name="no_filtering_lowcomplexity"
                and t5.taxon_id={taxon_id}) B group by cluster_name,cluster_id
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

    def get_genome_accession2description(self,
                                         taxon_id):
        
        sql = f'''
        select accession,description from genome_assembly_table
        where taxon_id={taxon_id}
        '''
        
        return {i[0]:i[1] for i in self.server.execute(sql,).fetchall()}
    
    def get_genome_accession2n_frameshifts(self,
                                           taxon_id,
                                           min_identity):
        
        sql = f'''select accession,count(*) as n from 
        (select t2.accession from homology_search t1 
        inner join genome_assembly_table t2 on t1.assembly_id=t2.assembly_id 
        where (alignment_id is not NULL 
        and percent_identity>={min_identity} 
        and taxon_id={taxon_id})) A 
        group by accession order by n DESC;
        '''
        
        return {i[0]:i[1] for i in self.server.execute(sql,).fetchall()}
           
        
    def get_genome_accession2statistics(self,
                                         taxon_id):
        
        sql = f'select accession, GC, cumulated_size, n_pseudo  from genome_assembly_table where taxon_id={taxon_id};'
        genome_data = self.server.execute(sql,).fetchall()
        accession2gc = {i[0]:i[1] for i in genome_data}
        accession2genome_size = {i[0]:i[2] for i in genome_data}
        accession2pseudo = {i[0]:i[3] for i in genome_data}
        
        return accession2gc, accession2genome_size, accession2pseudo
        
        
    def get_genome_accession2mlst(self,
                                  taxon_id):
        
        sql = f'''select accession,st from mlst t1 
             inner join genome_assembly_table t2 on t1.assembly_id=t2.assembly_id 
             where taxon_id={taxon_id};
          '''
        
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
        sql = f'''select accession,count(*) from (
            select distinct t5.accession,t5.description,hit_accession from homology_search t1
            inner join homology_search_db_xrefs t2 on t1.search_tool_id=t2.id
            inner join homology_search_db_xrefs t3 on t1.substitution_matrix_id=t3.id
            inner join homology_search_db_xrefs t4 on t1.filtering_id=t4.id
            inner join genome_assembly_table t5 on t1.assembly_id=t5.assembly_id
            inner join uniparc_entry t6 on t1.uniparc_id=t6.uniparc_id
            where t5.taxon_id={taxon_id}
            and t1.percent_identity>={min_identity}
            and uniparc_seq_coverage>={min_coverage}
            and t2.name in ("ssearch", "tfastx")
            and t3.name="BL50" 
            and t4.name="no_filtering_lowcomplexity"
            group by t5.accession,t1.hit_accession) A 
            group by accession;'''
        
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
                                           db_name):
        
        sql = f'''
        select distinct cluster_name from (
        select distinct cluster_id from uniparc2species t1
        inner join uniparc2mmseqs_90_80 t2 on t1.uniparc_id=t2.uniparc_id 
        where t1.species_taxon_id = {taxon_id}) A 
        inner join uniparc2mmseqs_90_80 B on A.cluster_id=B.cluster_id 
        inner join VF_table C on B.uniparc_id=C.uniparc_id 
        inner join VF_databases D on C.db_id= D.db_id 
        where D.db_name = "{db_name}"
        '''
        
        return [i[0] for i in self.server.execute(sql,).fetchall()]


    def get_db_cluster_id_list(self,
                               taxon_id, 
                               db_name):
        
        sql = f'''
        select distinct A.cluster_id from (
        select distinct cluster_id from uniparc2species t1
        inner join uniparc2mmseqs_90_80 t2 on t1.uniparc_id=t2.uniparc_id 
        where t1.species_taxon_id = {taxon_id}) A 
        inner join uniparc2mmseqs_90_80 B on A.cluster_id=B.cluster_id 
        inner join VF_table C on B.uniparc_id=C.uniparc_id 
        inner join VF_databases D on C.db_id= D.db_id 
        where D.db_name = "{db_name}"
        '''
        
        return [i[0] for i in self.server.execute(sql,).fetchall()]
   
        
        
    def VF_cluster_genome_counts(self,
                                 taxon_id, 
                                 min_identity=90, 
                                 min_coverage=80):
        '''
        Get number of hits for each cluster in each genome assembly
        Only possible for ssearch and not tfastx (would need to check position)
        '''
        
        # get cluster frequency
        sql_VF_freq = f'''select t6.cluster_name, t5.accession, count(*) as n from homology_search t1 
                    inner join homology_search_db_xrefs t2 on t1.search_tool_id=t2.id 
                    inner join homology_search_db_xrefs t3 on t1.substitution_matrix_id=t3.id 
                    inner join homology_search_db_xrefs t4 on t1.filtering_id=t4.id 
                    inner join genome_assembly_table t5 on t1.assembly_id=t5.assembly_id 
                    inner join {CLUSTERING_TABLE} t6 on t1.uniparc_id=t6.uniparc_id 
                    inner join {CLUSTERING2SPECIES_TABLE} t7 on t6.cluster_id=t7.cluster_id
                    where t5.taxon_id={taxon_id}
                    and t7.species_taxon_id={taxon_id}
                    and percent_identity>={min_identity}
                    and uniparc_seq_coverage>={min_coverage}
                    and t2.name in("ssearch") 
                    and t3.name="BL50" 
                    and t4.name="no_filtering_lowcomplexity"
                    group by t1.assembly_id,t6.cluster_id;'''
        
        return self.server.execute(sql_VF_freq,).fetchall()
        

    def get_species_name_from_taxid(self,taxon_id):
        
        sql = f'select species from ncbi_taxonomy where tax_id={taxon_id};'
        
        return self.server.execute(sql,).fetchall()[0][0]

    def get_species_VF_UP_count(self,
                                taxon_id):
        
        sql = f'''select count(*) from (select distinct t1.uniparc_id from {CLUSTERING_TABLE} t1 
        inner join uniparc2species t2 on t1.uniparc_id=t2.uniparc_id where species_taxon_id={taxon_id}) A;'''
        
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
 
    def get_mmseq_cluster_frequency(self, taxon_id):
        
        sql = f'''
        select cluster_id,count(*) as n from (select distinct t1.assembly_id,cluster_id from circos_genome2proteins t1 
        inner join circos_proteins2mmseqs_80_80 t2 on t1.protein_id=t2.protein_id 
        inner join genome_assembly_table t3 on t1.assembly_id=t3.assembly_id 
        where t3.taxon_id={taxon_id}) A group by A.cluster_id;
        '''
 
        return {i[0]:i[1] for i in self.server.execute(sql,).fetchall()}
    
    
    def get_assembly_VF_list(self, assembly_accession, identity_cutoff, taxon_id):
        
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
                    where accession="{assembly_accession}" 
                    and percent_identity>={identity_cutoff}
                    and uniparc_seq_coverage>=60 
                    and t7.name="ssearch" 
                    and t8.name="BL50" 
                    and t9.name="no_filtering_lowcomplexity"
                    and t2.taxon_id={taxon_id}
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
        VF_cluster2freq = self.get_cluster2frequency_within_species(taxon_id, 
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
        sql = f'select distinct cluster_name,gene from uniparc_consensus_annotation t1 inner join {CLUSTERING_TABLE} t2 on t1.uniparc_id=t2.uniparc_id where t2.cluster_name in ("{cluster_filter}")'
        print(sql)
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
        inner join circos_proteins2mmseqs_80_80 t3 on t1.protein_id=t3.protein_id
        where t2.accession="{assembly_accession}" and t2.taxon_id={taxon_id} order by record_id,start;
        '''
        print(sql)

        return pandas.read_sql(sql, self.conn)

    def get_rank_name(self, rank, taxon_id):
    
        sql = f'select distinct {rank}_name from uniparc2species where {rank}_taxon_id={taxon_id}'
    
        return self.server.execute(sql,).fetchall()[0][0]

    def get_species_phylogeny(self, taxon_id):
    
        sql = f'select phylogeny from species_phylogeny where taxon_id={taxon_id}'
    
        return self.server.execute(sql,).fetchall()[0][0]

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
    sql = 'select distinct accession,cluster_name from ' \
          ' (select distinct t2.accession,t1.VF_id ' \
          ' from genome_blast_results t1 inner join genome_assembly_table t2 ' \
          ' on t1.assembly_id=t2.assembly_id ' \
          ' where percent_identity>=%s and taxon_id=%s) A ' \
          ' inner join cluster_90 B on A.VF_id=B.VF_id;' % (min_identity, taxon_id)
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


def cluster2annotation(cluster_id_list, accessions=False):
    import os
    import sqlite3

    db = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + "/virulence.db"
    conn = sqlite3.connect(db)

    cursor = conn.cursor()
    
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

    cursor.execute(sql,)
    annotation = cursor.fetchall()
    cursor.execute(sql2,)
    species = cursor.fetchall()
    cursor.execute(sql3,)
    cluster_name2n_db = {i[0]:i[1] for i in cursor.fetchall()}

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
            print(existing_id, div_id)
            div = div.replace(existing_id, div_id)
        except IndexError:
            pass
    return div




#cluster_fam_freq(95)
