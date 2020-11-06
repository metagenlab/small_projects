import os
import sys

from BioSQL import BioSeqDatabase
from Bio import SeqIO
from Bio.SeqUtils import GC

import sqlite3
import pandas as pd


# This file defines a class DB, that encapsulates all the SQL requests
# necessary to create the chlamdb database.
# In the future, the goal is to import all database queries needed by the
# chlamdb website as methods of this class.
#
# This improves code readability by removing SQL queries from the main python
# code and more importantly, it would allow to change of database without having
# to modify chlamdb's code.


# to litteral
# encases the string into quotes
def quote(v):
    return f"\"{v}\""

class DB:

    def __init__(self, server, db_name):
        self.server = server
        self.db_name = db_name
        self.conn_ncbi_taxonomy = None
        self.conn_refseq = None

    def create_indices_on_cds(self):
        sql_index1 = 'create index ftgcga on feature_tables_genomes_cds(genome_accession)'    
        sql_index2 = 'create index ftgctx on feature_tables_genomes_cds(taxon_id)'
        sql_index3 = 'create index ftgrga on feature_tables_genomes_rrna(taxon_id)'
        sql_index4 = 'create index ftgrtx on feature_tables_genomes_rrna(genome_accession)'    
        sql_index5 = 'create index ftcain on feature_tables_cds_accessions(id_name)'
        sql_index6 = 'create index ftcait on feature_tables_cds_accessions(id_type)'
        
        self.server.adaptor.execute(sql_index1,)
        self.server.adaptor.execute(sql_index2,)
        self.server.adaptor.execute(sql_index3,)
        self.server.adaptor.execute(sql_index4,)
        self.server.adaptor.execute(sql_index5,)
        self.server.adaptor.execute(sql_index6,)

    def get_taxid_from_accession(self, accession):
        sql = (
            f"SELECT taxon_id "
            f"FROM bioentry t1 INNER JOIN biodatabase t2 "
            f"ON t1.biodatabase_id = t2.biodatabase_id "
            f"WHERE t2.name={quote(self.db_name)} AND t1.accession={quote(accession)}"
        )
        return self.server.adaptor.execute_and_fetchall(sql,)[0][0]

    def create_seq_hash_to_seqid(self, to_load):
        sql = (
            f"CREATE TABLE sequence_hash_dictionnary (hsh INTEGER, seqid INTEGER," 
            " PRIMARY KEY(hsh, seqid), "
            " FOREIGN KEY(seqid) REFERENCES seqfeature(seqfeature_id));"
        )
        self.server.adaptor.execute(sql,)
        self.load_data_into_table("sequence_hash_dictionnary", to_load)

        sql = "CREATE INDEX shd_hsh on sequence_hash_dictionnary (hsh)"
        self.server.adaptor.execute(sql)
        sql = "CREATE INDEX shd_seqid on sequence_hash_dictionnary (seqid)"
        self.server.adaptor.execute(sql)

    def get_seq_hash_to_seqid(self, to_load):
        query = (
            "SELECT * FROM sequence_hash_dictionnary;"
        )
        results = self.server.adaptor.execute_and_fetchall(query)
        hsh_results = {}
        for line in results:
            hsh_results[line[0]] = line[1]
        return hsh_results

    # Returns all refseq hits associated with a given orthogroup
    def get_diamond_match_for_og(self, og):
        query = (
            "SELECT match_id.accession, match_id.taxid, hit.hit_count "
            "FROM og_hits " 
            "INNER JOIN sequence_hash_dictionnary AS hsh ON hsh.seqid = og_hits.seqid "
            "INNER JOIN diamond_refseq AS hit ON hit.seq_hash = hsh.hsh "
            "INNER JOIN diamond_refseq_match_id AS match_id ON hit.sseqid = match_id.match_id "
            f"WHERE og_hits.orthogroup={og} "
            "GROUP BY match_id.accession, match_id.taxid ORDER BY hit.hit_count ASC;" 
        )
        results = self.server.adaptor.execute_and_fetchall(query)
        gen = ((line[0], line[1]) for line in results)
        return gen

    def get_all_orthogroups(self, min_size=None):
        query = (
            "SELECT orthogroup, COUNT(*) "
            "FROM og_hits "
            "GROUP BY orthogroup;" 
        )
        results = self.server.adaptor.execute_and_fetchall(query)
        acc = (line[0] for line in results if min_size==None or line[1]>=min_size)
        return acc

    def get_all_sequences_for_orthogroup(self, orthogroup):
        from Bio import SeqRecord
        from Bio.Seq import Seq

        query = (
            "SELECT locus.value, seq.value "
            "FROM og_hits AS ortho "
            "INNER JOIN seqfeature_qualifier_value AS seq ON seq.seqfeature_id=ortho.seqid "
            " INNER JOIN term AS seq_term "
            "   ON seq_term.term_id=seq.term_id AND seq_term.name=\"translation\""
            "INNER JOIN seqfeature_qualifier_value locus ON locus.seqfeature_id=ortho.seqid "
            " INNER JOIN term as locus_term "
            "   ON locus_term.term_id=locus.term_id AND locus_term.name=\"locus_tag\" "
            f"WHERE ortho.orthogroup = {orthogroup};"
        )
        results = self.server.adaptor.execute_and_fetchall(query,)
        tab = []
        for result in results:
            record = SeqRecord.SeqRecord(Seq(result[1]), id=result[0])
            tab.append(record)
        return tab

    def create_diamond_refseq_match_id(self):
        query = (
            "CREATE TABLE diamond_refseq_match_id ( "
            "match_id INT, accession VARCHAR(200), taxid INT, description TEXT, length INT, "
            "PRIMARY KEY(match_id));"
        )
        self.server.adaptor.execute(query,)

    def load_diamond_refseq_match_id(self, data):
        self.load_data_into_table("diamond_refseq_match_id", data)

    def create_diamond_refseq_match_id_indices(self):
        query = (
            "CREATE INDEX drmii ON diamond_refseq_match_id(match_id);"
        )
        self.server.adaptor.execute(query,)

    def create_refseq_hits_taxonomy(self):
        sql = (
            "CREATE TABLE IF NOT EXISTS refseq_hits_taxonomy(taxid INT, superkingdom INT, phylum INT, "
            "class INT, order_id INT, family INT, genus INT, specie INT, PRIMARY KEY(taxid));"
        )
        self.server.adaptor.execute(sql)

    def create_refseq_hits_taxonomy_indices(self):
        sql = (
            "CREATE INDEX rhti ON refseq_hits_taxonomy(taxid)"
        )
        self.server.adaptor.execute(sql)

    def create_taxonomy_mapping(self, hsh):
        sql = (
            "CREATE TABLE taxonomy_mapping(taxid INT, "
            "rank TEXT, value TEXT, PRIMARY KEY(taxid));"
        )
        self.server.adaptor.execute(sql)
        lst_values = []
        for key, (rank, value) in hsh.items():
            lst_values.append( (key, rank, value) )
        self.load_data_into_table("taxonomy_mapping", lst_values)

    def load_data_into_table(self, table, data):
        fmt_string = ", ".join("?" for i in range(len(data[0])))
        sql_string = f"INSERT into {table} VALUES ({fmt_string});"
        self.server.adaptor.executemany(sql_string, data)

    def load_refseq_hits_taxonomy(self, data):
        self.load_data_into_table("refseq_hits_taxonomy", data)

    def get_all_taxids(self):
        query = (
            "SELECT DISTINCT taxid from diamond_refseq_match_id;"
        )
        results = self.server.adaptor.execute_and_fetchall(query,)
        taxids = []
        for line in results:
            taxids.append(int(line[0]))
        return taxids

    # Utility class to make user code more readable
    class Taxon:
        hsh_taxo_key = {
                "superkingdom" : (1, 2),
                "phylum" :       (3, 4),
                "class" :        (5, 6),
                "order" :        (7, 8),
                "family" :       (9, 10),
                "genus" :        (11, 12),
                "species" :      (13, 14)
        }

        def __init__(self, line):
            self.raw_result = line

        def phylum(self):
            idx = DB.Taxon.hsh_taxo_key["phylum"][0]
            return self.raw_result[idx]

        def taxid(self):
            return self.raw_result[0]

        def update_hash(self, curr_hash):
            for rank, (name_idx, taxid_idx) in DB.Taxon.hsh_taxo_key.items():
                if self.raw_result[taxid_idx] in curr_hash:
                    continue
                curr_hash[self.raw_result[taxid_idx]] = (rank, self.raw_result[name_idx])

        def get_all_taxids(self):
            all_taxids = [self.taxid()]
            for taxon_name, (idx_name, idx_taxid) in DB.Taxon.hsh_taxo_key.items():
                all_taxids.append(self.raw_result[idx_taxid])
            return all_taxids

    def get_linear_taxonomy(self, args, taxids):
        conn_refseq = sqlite3.connect(args["databases_dir"] + "/ncbi-taxonomy/linear_taxonomy.db")
        cursor = conn_refseq.cursor()

        query_string = ",".join([str(i) for i in taxids])
        query = (
            "SELECT tax_id, `superkingdom`, superkingdom_taxid, "
            " `phylum`, phylum_taxid, `class`, class_taxid, "
            " `order`, order_taxid, `family`, family_taxid, "
            " `genus`, genus_taxid, `species`, species_taxid "
            f"FROM ncbi_taxonomy WHERE tax_id IN ({query_string});"
        )
        results = cursor.execute(query, ).fetchall()
        return (DB.Taxon(line) for line in results)

    def get_accession_to_taxid(self, accession, params):
        # reuse the connection to the database if already open
        if self.conn_ncbi_taxonomy == None:
            dtb_path = params["databases_dir"] + "/ncbi-taxonomy/prot_accession2taxid.db"
            self.conn_ncbi_taxonomy = sqlite3.connect(dtb_path)

        cursor = self.conn_ncbi_taxonomy.cursor()
        query_string = ",".join(["\"" + a + "\"" for a in accession])
        query = (
            f"SELECT accession, taxid FROM accession2taxid "
            f"WHERE accession IN ({query_string});" 
        )
        results = cursor.execute(query,).fetchall()
        hsh_results = {}
        for line in results:
            hsh_results[line[0]] = int(line[1])
        return hsh_results

    # NOTE: need to check which indices are necessary to add to this table
    def load_cog_hits(self, data):
        sql = (
            "CREATE TABLE cog_hits (hsh INTEGER, cog_id INT, evalue FLOAT);"
            # " FOREIGN KEY(cog_id) REFERENCES cog_names(cog_id)); "
        )
        self.server.adaptor.execute(sql,)

        sql = "CREATE INDEX chi_hsh ON cog_hits(hsh);"
        self.load_data_into_table("cog_hits", data)

    def load_cog_fun_data(self, data):
        sql = (
            "CREATE TABLE cog_functions (function TEXT, description TEXT);"
        )
        self.server.adaptor.execute(sql)
        self.load_data_into_table("cog_functions", data)


    def load_cdd_to_cog(self, data):
        sql = (
            "CREATE TABLE cdd_to_cog ( "
            "cdd INT, cog_id INT);"
        )
        self.server.adaptor.execute(sql)
        self.load_data_into_table("cdd_to_cog", data)


    def get_cdd_to_cog(self):
        query = (
            "SELECT cdd, cog_id "
            "FROM cdd_to_cog;"
        )
        results = self.server.adaptor.execute_and_fetchall(query)
        return DB.to_hsh(results, [1])

    
    def load_cog_ref_data(self, data):
        sql = (
            "CREATE TABLE cog_names (cog_id INTEGER, function TEXT, description TEXT, "
            "PRIMARY KEY(cog_id));"
        )
        self.server.adaptor.execute(sql)
        self.load_data_into_table("cog_names", data)


    def load_ko_pathway(self, data):
        sql = (
            "CREATE TABLE ko_pathway_def ("
            "pathway_id INTEGER, desc TEXT, "
            "PRIMARY KEY(pathway_id));"
        )
        self.server.adaptor.execute(sql)
        self.load_data_into_table("ko_pathway_def", data)
        sql = "CREATE INDEX kpd_i ON ko_pathway_def(pathway_id);"
        self.server.adaptor.execute(sql)


    def load_ko_module_classes(self, data):
        sql = (
            "CREATE TABLE ko_class( "
            "class_id INTEGER, descr TEXT, PRIMARY KEY(class_id));"
        )
        self.server.adaptor.execute(sql)
        self.load_data_into_table("ko_class", data)
        sql = "CREATE INDEX kc_i ON ko_class(class_id);"
        self.server.adaptor.execute(sql)


    def load_ko_module(self, data):
        sql = (
            "CREATE TABLE ko_module_def ("
            "module_id INTEGER, desc TEXT, definition TEXT, class INT, subclass INT, "
            "PRIMARY KEY(module_id), FOREIGN KEY(class) REFERENCES ko_class(class_id), "
            "FOREIGN KEY(subclass) REFERENCES ko_class(class_id));"
        )
        self.server.adaptor.execute(sql)
        self.load_data_into_table("ko_module_def", data)
        sql = "CREATE INDEX kmd_i ON ko_module_def(module_id);"
        self.server.adaptor.execute(sql)

    def load_ko_to_pathway(self, data):
        sql = (
            "CREATE TABLE ko_to_pathway ("
            "ko_id INT, pathway_id INT, "
            "PRIMARY KEY(ko_id, pathway_id), "
            "FOREIGN KEY(ko_id) REFERENCES ko_def(ko_id), "
            "FOREIGN KEY(pathway_id) REFERENCES ko_pathway_def(pathway_id));"
        )
        self.server.adaptor.execute(sql)
        self.load_data_into_table("ko_to_pathway", data)
        sql = "CREATE INDEX ktpk_i ON ko_to_pathway(ko_id);"
        self.server.adaptor.execute(sql)
        sql = "CREATE INDEX ktpp_i ON ko_to_pathway(pathway_id);"
        self.server.adaptor.execute(sql)

    def load_ko_to_module(self, data):
        sql = (
            "CREATE TABLE ko_to_module ("
            "ko_id INT, module_id INT, "
            "PRIMARY KEY(ko_id, module_id), "
            "FOREIGN KEY(ko_id) REFERENCES ko_def(ko_id), "
            "FOREIGN KEY(module_id) REFERENCES ko_module_def(module_id));"
        )
        self.server.adaptor.execute(sql)
        self.load_data_into_table("ko_to_module", data)
        sql = "CREATE INDEX ktmk_i ON ko_to_module(ko_id);"
        self.server.adaptor.execute(sql)
        sql = "CREATE INDEX ktmm_i ON ko_to_module(module_id);"
        self.server.adaptor.execute(sql)

    # Note: EC to add separately?
    def load_ko_def(self, data):
        sql = (
            "CREATE TABLE ko_def ( "
            "ko_id, descr TEXT, PRIMARY KEY(ko_id));"
        )
        self.server.adaptor.execute(sql)
        self.load_data_into_table("ko_def", data)


    def load_ko_hits(self, data):
        sql = (
            "CREATE TABLE ko_hits ("
            "hsh INTEGER, ko_id INT, threshold FLOAT, score FLOAT, evalue FLOAT, "
            "PRIMARY KEY(hsh, ko_id),"
            "FOREIGN KEY(ko_id) REFERENCES ko_def(ko_id));"
        )
        self.server.adaptor.execute(sql)
        self.load_data_into_table("ko_hits", data)
        sql = (
            "CREATE INDEX khi_i ON ko_hits(ko_id);"
        )
        self.server.adaptor.execute(sql)

    def get_all_modules_definition(self):
        query = "SELECT module_id, definition FROM ko_module_def;"
        results = self.server.adaptor.execute_and_fetchall(query)
        return [(line[0], line[1]) for line in results]

    def get_ko_total_count(self, ko_ids):
        entries = ",".join("?" for i in ko_ids)
        query = (
            "SELECT ko_hits.ko_id, seq.bioentry_id "
            "FROM ko_hits AS ko_hits "
            "INNER JOIN sequence_hash_dictionnary AS hsh ON hsh.hsh = ko_hits.hsh "
            "INNER JOIN seqfeature AS seq ON seq.seqfeature_id = hsh.seqid "
            f"WHERE ko_hits.ko_id IN ({entries})"
            "GROUP BY ko_hits.ko_id, seq.bioentry_id;"
        )
        results = self.server.adaptor.execute_and_fetchall(query, ko_ids)
        df = DB.to_pandas_frame(results, ["KO", "bioentry"])
        return df.groupby("KO").count()


    def to_hsh(results, val_range, has_multiple=False):
        hsh = {}
        for line in results:
            if len(val_range)==1:
                data = line[val_range[0]]
            else:
                data = [line[i] for i in val_range]
            if not has_multiple:
                hsh[line[0]] = data
            else:
                curr = hsh.get(line[0], [])
                curr.append(data)
                hsh[line[0]] = curr
        return hsh


    def get_ko_desc(self, ko_ids):
        entries = ",".join("?" for i in ko_ids)
        query = (
            "SELECT ko.ko_id, ko.descr "
            "FROM ko_def as ko "
            f"WHERE ko.ko_id IN ({entries});"
        )
        results = self.server.adaptor.execute_and_fetchall(query, ko_ids)
        return DB.to_hsh(results, [1])


    def get_ko_pathways(self, ko_ids):
        entries = ",".join("?" for i in ko_ids)
        query = (
            "SELECT ktp.ko_id, path.pathway_id, path.desc "
            "FROM ko_to_pathway AS ktp "
            "INNER JOIN ko_pathway_def AS path ON path.pathway_id = ktp.pathway_id "
            f"WHERE ktp.ko_id IN ({entries});"
        )
        results = self.server.adaptor.execute_and_fetchall(query, ko_ids)
        return DB.to_hsh(results, [1,2], has_multiple=True)


    def get_ko_modules(self, ko_ids):
        entries = ",".join("?" for i in ko_ids)
        query = (
            "SELECT ktm.ko_id, mod.module_id, mod.desc "
            "FROM ko_to_module AS ktm "
            "INNER JOIN ko_module_def AS mod ON mod.module_id = ktm.module_id "
            f"WHERE ktm.ko_id IN ({entries});"
        )
        results = self.server.adaptor.execute_and_fetchall(query, ko_ids)
        return DB.to_hsh(results, [1,2], has_multiple=True)


    def get_ko_count(self, bioentries):
        entries = ",".join("?" for i in bioentries)
        query = (
            "SELECT feature.bioentry_id, hit.ko_id, count(*) "
            "FROM seqfeature AS feature "
            "INNER JOIN sequence_hash_dictionnary AS hsh ON feature.seqfeature_id = hsh.seqid "
            "INNER JOIN ko_hits AS hit ON hit.hsh = hsh.hsh "
            f"WHERE feature.bioentry_id IN ({entries})"
            "GROUP BY feature.bioentry_id, hit.ko_id;"
        )
        results = self.server.adaptor.execute_and_fetchall(query, bioentries)
        return DB.to_pandas_frame(results, ["bioentry", "KO", "count"])


    def get_seqids_for_ko(self, ko_ids):
        entries = ",".join("?" for i in bioentries)
        query = (
            "SELECT seqid"
            "FROM ko_hits "
            "INNER JOIN sequence_hash_dictionnary"
            f"WHERE ko_id IN ({entries});"
        )
        return [line[0] for line in results].unique()


    def get_hsh_locus_to_seqfeature_id(self):
        query = (
            "SELECT t2.value, t1.seqfeature_id "
            "FROM seqfeature as t1 "
            "INNER JOIN seqfeature_qualifier_value AS t2 ON t1.seqfeature_id=t2.seqfeature_id "
            "INNER JOIN term as t3 on t3.term_id=t2.term_id AND t3.name=\"locus_tag\""
        )
        results = self.server.adaptor.execute_and_fetchall(query,)
        hsh_results = {}
        for line in results:
            hsh_results[line[0]] = int(line[1])
        return hsh_results

    def get_locus_to_og(self):
        query = (
            f"SELECT locus.value, og.value"
            f"FROM seqfeature_qualifier_value as og "
            f"INNER JOIN term AS term_og ON term_og.term_id=og.term_id AND term_og.name=\"orthogroup\" "
            f"INNER JOIN seqfeature_qualifier_value as locus ON locus.seqfeature_id=og.seqfeature_id "
            f"INNER JOIN term AS term_locus ON term_locus.term_id=locus.term_id "
            f"  AND term_locus.name=\"locus_tag\";"
        )
        query_results = self.server.adaptor.execute_and_fetchall(query,)
        locus_to_og = {}
        for locus, og in query_results:
            locus_to_og[locus] = int(og)
        return locus_to_og

    def create_new_og_matrix(self):
        query = (
            f"CREATE TABLE orthology_identity ( "
            f"orthogroup INT, "
            f"id_1 INT, id_2 INT, identity FLOAT(5), "
            f"PRIMARY KEY(orthogroup, id_1, id_2), "
            f"FOREIGN KEY(id_1) REFERENCES seqfeature(seqfeature_id)"
            f"FOREIGN KEY(id_2) REFERENCES seqfeature(seqfeature_id))"
        )
        self.server.adaptor.execute(query, )

    # To be removed
    # For each taxon, get the number of protein present in each of the existing
    # ortholog group
    # Returns a table of maps, organized as follows:
    #  group 0: hsh taxid -> number of proteins of taxid in group 0
    #  group 1: hsh taxid -> number of proteins of taxid in group 1
    #  ... 
    #  group N: hsh taxid -> number of proteins of taxid in group N
    def get_orthogroup_count_table(self):
        query = (
            "SELECT value, taxon_id, count(*) "
            "FROM seqfeature_qualifier_value as ortho_table "
            "INNER JOIN term ON ortho_table.term_id = term.term_id AND term.name=\"orthogroup\" "
            "INNER JOIN seqfeature AS feature ON feature.seqfeature_id=ortho_table.seqfeature_id "
            "INNER JOIN bioentry AS entry ON entry.bioentry_id=feature.bioentry_id "
            "GROUP BY value, taxon_id; "
        )
        results = self.server.adaptor.execute_and_fetchall(query,)
        size = max([int(group) for group, taxid, cnt in results])
        arr_group = [None]*(size+1)
        for line in results:
            group, taxid, cnt = int(line[0]), line[1], int(line[2])
            hsh_taxid_to_count = arr_group[group]
            if hsh_taxid_to_count==None:
                hsh_taxid_to_count = {}
                arr_group[group] = hsh_taxid_to_count
            hsh_taxid_to_count[taxid] = cnt
        return arr_group

    def load_og_matrix(self, matrix):
        for orthogroup_id, id_1, id_2, identity in matrix:
            query = (
                f"INSERT INTO orthology_identity "
                f"VALUES ({orthogroup_id}, {id_1}, {id_2}, {identity});"
            )
            self.server.adaptor.execute(query,)

    def create_BBH_phylogeny_table(self, data):
        sql = (
            "CREATE TABLE BBH_phylogeny (orthogroup_id INTEGER, tree TEXT, "
            " PRIMARY KEY(orthogroup_id));"
        )
        self.server.adaptor.execute(sql,)
        self.load_data_into_table("BBH_phylogeny", data)

    def create_gene_phylogeny_table(self, data):
        sql = (
            "CREATE TABLE gene_phylogeny (orthogroup_id INTEGER, tree TEXT, "
            "PRIMARY KEY(orthogroup_id)); "
        )
        self.server.adaptor.execute(sql,)
        self.load_data_into_table("gene_phylogeny", data)

    def create_og_matrix_indices(self):
        sql_1 = "CREATE INDEX oio ON orthology_identity(orthogroup);"
        self.server.adaptor.execute(sql_1,)

    def load_og_hits(self, lst_to_load):
        query = (
            "CREATE TABLE og_hits ( "
            "seqid INTEGER, orthogroup INTEGER, "
            "PRIMARY KEY(seqid), "
            "FOREIGN KEY(seqid) REFERENCES seqfeature(seqfeature_id));"
        )
        self.server.adaptor.execute(query,)
        query = "INSERT INTO og_hits VALUES(?, ?);"
        self.server.adaptor.executemany(query, lst_to_load)
        query = "CREATE INDEX og_og_idx ON og_hits(seqid);"
        self.server.adaptor.execute(query)
        query = "CREATE INDEX og_seqid_idx ON og_hits(orthogroup);"
        self.server.adaptor.execute(query)

    def taxon_ids(self):
        query = "SELECT taxon_id FROM taxon";
        result = self.server.adaptor.execute_and_fetchall(query,)
        arr_taxon_ids = []
        for line in result:
            arr_taxon_ids.append(int(line[0]))
        return arr_taxon_ids

    def set_status_in_config_table(self, status_name, status_val):
        sql = f"update biodb_config set status={status_val} where name={quote(status_name)};"
        self.server.adaptor.execute(sql,)

    def create_biosql_database(self, args):
        self.server.new_database(self.db_name)

    # NOTE:
    # As the primary key in sequence_hash_dictionnary is a tuple of hsh and seqid, we would 
    # need to have an entry in diamond_refseq for every possible combination of (hsh, seqid) to 
    # be able to use FOREIGN KEY. To avoid that and spare some disk space, foreign keys were not
    # used..
    def create_refseq_hits_table(self):
        sql = (
            f"CREATE TABLE IF NOT EXISTS diamond_refseq(hit_count INT, seq_hash INTEGER, "
            f"sseqid INT, pident FLOAT, length INT, mismatch INT, "
            f"gapopen INT, qstart INT, qend INT, sstart INT, "
            f"send INT, evalue FLOAT, bitscore FLOAT, "
            "FOREIGN KEY(sseqid) REFERENCES diamond_refseq_match_id(match_id), "
            f"PRIMARY KEY (hit_count, seq_hash));"
        )
        self.server.adaptor.execute(sql,)

    def load_refseq_hits(self, data):
        self.load_data_into_table("diamond_refseq", data)

    def create_refseq_hits_indices(self):
        sql = (
            "CREATE INDEX dri ON diamond_refseq(seq_hash);"
        )
        self.server.adaptor.execute(sql,)

    def n_accession_to_n_tRNA(self):
        query = (
            "SELECT accession, count(*) FROM seqfeature AS seq "
            " INNER JOIN term AS t ON t.term_id = seq.type_term_id "
            " INNER JOIN bioentry AS entry ON seq.bioentry_id=entry.bioentry_id "
            " AND t.name=\"tRNA\" GROUP BY accession; "
        )
        n_tRNA_results = self.server.adaptor.execute_and_fetchall(query)
        results = {}
        for line in n_tRNA_results:
            results[line[0]] = int(line[1])
        return results
    
    def get_genomes_sequences(self):
        query = (
            "SELECT bioentry_id, seq FROM biosequence;"
        )
        results = self.server.adaptor.execute_and_fetchall(query)
        hsh_results = {}
        for str_id, seq in results:
            hsh_results[int(str_id)] = seq
        return hsh_results

    def get_term_id(self, term, create_if_absent=False, ontology="Annotation Tags"):
        query = f"SELECT term_id FROM term WHERE name=\"{term}\";"
        result = self.server.adaptor.execute_and_fetchall(query)
        gc_term_id = None
        if len(result) == 0 and create_if_absent:
            sql1 = (
                "SELECT ontology_id FROM ontology "
                f"WHERE name=\"{ontology}\";"
            )
            ontology_id = self.server.adaptor.execute_and_fetchall(sql1)[0][0]
            sql = (
                f"INSERT INTO term (name, ontology_id) VALUES (\"{term}\", {ontology_id});"
            )
            self.server.adaptor.execute(sql)
            gc_term_id = self.server.adaptor.cursor.lastrowid
        else:
            gc_term_id = result[0][0]
        return gc_term_id

    def get_genomes_description(self, entries=None, indexing="bioentry", exclude_plasmids=False,
            indexing_type="str"):
        if indexing != "bioentry" and indexing != "accession":
            raise RuntimeError(f"{bioentry} indexing is not supported")

        if indexing_type != "int" and indexing_type != "str":
            raise RuntimeError(f"{indexing_type} not supported, must be int or str")

        entry_filter = ""
        if entries != None:
            string = ",".join(["?"] * len(entries))
            entry_filter = f"WHERE bioentry_id IN ({string}) "

        selection = ("accession", "bioentry_id")[indexing=="bioentry"]
        query = (
            f"SELECT {selection}, description "
            " FROM bioentry "
            f"{entry_filter}"
            " ORDER BY description;"
        )
        if entries == None:
            results = self.server.adaptor.execute_and_fetchall(query)
        else:
            results = self.server.adaptor.execute_and_fetchall(query, entries)

        hsh_results = {}
        for line in results:
            if exclude_plasmids and ("plasmid" in line[1] or "phage" in line[1]):
                continue
            idx = line[0] if indexing_type == "int" else str(line[0])
            hsh_results[idx] = line[1].strip()
        return hsh_results

    # Used for example in the 'home' page of Chlamdb
    # Returns a table with a row for each genome, with the following entries
    # - accession
    # - GC
    # - number of proteins
    # - number of contigs
    # - size
    # - percent non-coding
    # - description
    def get_genomes_infos(self, filter_out_plasmids=True):
        query = (
            "SELECT accession, COUNT(*) from seqfeature AS seq"
            " INNER JOIN term AS cds ON cds.term_id = seq.type_term_id AND cds.name=\"CDS\" "
            " INNER JOIN bioentry AS entry ON entry.bioentry_id = seq.bioentry_id " 
            " GROUP BY accession ORDER BY accession;"
        )
        n_prot_results = self.server.adaptor.execute_and_fetchall(query)

        query = (
            "SELECT accession, gc.value, l.value, contigs.value, cd.value, description "
            "FROM bioentry AS entry "
            " INNER JOIN bioentry_qualifier_value AS gc ON entry.bioentry_id = gc.bioentry_id "
            "  INNER JOIN term as gc_term ON gc_term.term_id = gc.term_id AND gc_term.name = \"gc\" "
            " INNER JOIN bioentry_qualifier_value AS l ON l.bioentry_id=entry.bioentry_id "
            "  INNER JOIN term as l_term ON l_term.term_id = l.term_id AND l_term.name = \"length\" "
            " INNER JOIN bioentry_qualifier_value AS contigs ON contigs.bioentry_id=entry.bioentry_id "
            "  INNER JOIN term as contigs_term ON contigs_term.term_id = contigs.term_id " 
            "     AND contigs_term.name = \"n_contigs\""
            " INNER JOIN bioentry_qualifier_value AS cd ON cd.bioentry_id=entry.bioentry_id "
            "  INNER JOIN term as cd_term ON cd_term.term_id = cd.term_id " 
            "     AND cd_term.name = \"coding_density\" "
            " GROUP BY accession ORDER BY accession;"
        )
        all_other_results = self.server.adaptor.execute_and_fetchall(query)

        # Note: both queries return the values ordered by accession
        # TODO: this should not happen, but it would be a good idea to add
        # an assert to make sure the entries correspond
        assembled_results = []
        for index, value in enumerate(n_prot_results):
            other_entry = all_other_results[index]
            description = other_entry[5]
            if filter_out_plasmids and "plasmid" in description:
                continue

            n_prot = int(value[1])
            gc, length = round(float(other_entry[1]),2), int(other_entry[2])
            n_contigs, non_coding_density = int(other_entry[3]), round(100.0 - float(other_entry[4]), 2)
            entry = (value[0], gc, n_prot, n_contigs, length, non_coding_density, description)
            assembled_results.append(entry)
        return assembled_results

    def load_genomes_gc(self, gcs):
        term_id = self.get_term_id("gc", create_if_absent=True)
        sql = "INSERT INTO bioentry_qualifier_value (bioentry_id, term_id, value) VALUES (?, ?, ?);";
        for bioentry_id, gc in gcs:
            self.server.adaptor.execute(sql , (bioentry_id, term_id, gc))

    def load_genomes_lengths(self, lengths):
        term_id = self.get_term_id("length", create_if_absent=True)
        sql = "INSERT INTO bioentry_qualifier_value (bioentry_id, term_id, value) VALUES (?, ?, ?);";
        for bioentry_id, length in lengths:
            self.server.adaptor.execute(sql , (bioentry_id, term_id, length))

    def load_genomes_n_contigs(self, n_contigs):
        term_id = self.get_term_id("n_contigs", create_if_absent=True)
        sql = "INSERT INTO bioentry_qualifier_value (bioentry_id, term_id, value) VALUES (?, ?, ?);";
        for bioentry_id, n_contig in n_contigs:
            self.server.adaptor.execute(sql , (bioentry_id, term_id, n_contig))

    def load_genomes_coding_density(self, coding_density):
        term_id = self.get_term_id("coding_density", create_if_absent=True)
        sql = "INSERT INTO bioentry_qualifier_value (bioentry_id, term_id, value) VALUES (?, ?, ?);";
        for bioentry_id, density in coding_density:
            self.server.adaptor.execute(sql , (bioentry_id, term_id, density))

    def get_coding_region_total_length(self):
        query = (
            "SELECT gen.bioentry_id, SUM(end_pos-start_pos) "
            " FROM location AS loc "
            " INNER JOIN seqfeature AS seq ON loc.seqfeature_id=seq.seqfeature_id "
            " INNER JOIN bioentry AS gen ON gen.bioentry_id=seq.bioentry_id "
            " INNER JOIN term AS t ON seq.type_term_id=t.term_id AND t.name=\"gene\" "
            " GROUP BY gen.bioentry_id; "
        )
        results = self.server.adaptor.execute_and_fetchall(query)
        hsh_entry_to_sum = {}
        for str_id, str_sum in results:
            hsh_entry_to_sum[int(str_id)] = int(str_sum)
        return hsh_entry_to_sum

    def get_accession_to_entry(self):
        query = (
            "SELECT accession, bioentry_id FROM bioentry;"
        )
        results = self.server.adaptor.execute_and_fetchall(query)
        hsh_results = {}
        for line in results:
            hsh_results[line[0]] = line[1]
        return hsh_results

    # Return the number of proteins, indexed either by
    # - accession
    # - bioentry_id
    def n_CDS(self, indexing="bioentry", taxons=None, bioentries=None):

        if taxons!=None and bioentries!=None:
            raise RuntimeError("Specifying both taxons and bioentries not supported")

        # This needs some testing. Also needs query sanitizing
        filtering = ""
        if taxons != None:
            j = ",".join([str(taxon) for taxon in taxons])
            filtering = (
                " INNER JOIN taxon AS t ON t.taxon_id=entry.taxon_id "
                f" AND ncbi_taxon_id IN ({j});"
            )
        elif bioentries != None:
            j = ",".join([str(entry) for entry in bioentries])
            filtering = f"WHERE entry.bioentry_id IN {j}"

        index = None
        if indexing == "bioentry":
            index = "entry.bioentry_id"
        elif indexing == "accession":
            index = "accession"
        elif indexing == "ncbi":
            index = "ncbi_taxon_id"
        else:
            raise RuntimeError(f"Unsupported indexing method: {indexing}")

        query = (
            f"SELECT {index}, COUNT(*) FROM seqfeature AS prot "
            " INNER JOIN term AS t ON t.term_id = prot.type_term_id AND t.name=\"CDS\" "
            " INNER JOIN bioentry AS entry ON entry.bioentry_id = prot.bioentry_id "
            f"{filtering}"
            f"GROUP BY {index};"
        )
        results = self.server.adaptor.execute_and_fetchall(query)
        hsh_results = {}
        for index, str_count in results:
            hsh_results[index] = int(str_count)
        return hsh_results

    def get_cog_code_description(self):
        query = (
            "SELECT function, description FROM cog_functions;"
        )
        results = self.server.adaptor.execute_and_fetchall(query)
        hsh_results = {}
        for line in results:
            hsh_results[line[0].strip()] = line[1].strip()
        return hsh_results

    # Note: need to check whether this is a more correct way to proceed
    # as my Rhabdo genomes currently all have the same taxon id, the following
    # code wouldn't work
    # sql_n_genomes = 'select count(*) from (select distinct taxon_id from bioentry t1 ' \
    # ' inner join biodatabase t2 on t1.biodatabase_id=t2.biodatabase_id where t2.name="%s") A;' % biodb
    def get_n_genomes(self):
        query = "SELECT COUNT(*) FROM bioentry;"
        result = self.server.adaptor.execute_and_fetchall(query)
        return result[0][0]

    def load_reference_phylogeny(self, tree):
        sql = "CREATE TABLE IF NOT EXISTS reference_phylogeny (tree TEXT);"
        self.server.adaptor.execute(sql,)

        sql = f"INSERT INTO reference_phylogeny VALUES (\"{tree}\");"
        self.server.adaptor.execute(sql,)

    def get_reference_phylogeny(self):
        sql = "SELECT tree FROM reference_phylogeny;"
        values = self.server.adaptor.execute_and_fetchall(sql)
        return values[0][0]

    def get_bioentry_id_for_record(self, record):
        locus_tag = record["gene"]["locus_tag"]
        sql = (
            "SELECT bioentry_id FROM seqfeature_qualifier_value AS tag"
            " INNER JOIN seqfeature AS seq ON seq.seqfeature_id = tag.seqfeature_id "
            f" where value = \"{locus_tag}\";" 
        )
        results = self.server.adaptor.execute_and_fetchall(sql)
        return results[0][0]

    def load_filenames(self, data):
        sql = (
            "CREATE TABLE filenames (bioentry_id INTEGER, filename TEXT, "
            " FOREIGN KEY(bioentry_id) REFERENCES bioentry(bioentry_id));"
        )
        self.server.adaptor.execute(sql)
        self.load_data_into_table("filenames", data)

    def load_checkm_results(self, data):
        sql = (
            "CREATE TABLE checkm (bioentry_id INTEGER, completeness FLOAT, "
            " contamination FLOAT, FOREIGN KEY(bioentry_id) REFERENCES bioentry(bioentry_id));"
        )
        self.server.adaptor.execute(sql)
        self.load_data_into_table("checkm", data)

    def get_checkm_results(self):
        query = (
            "SELECT accession, completeness, contamination "
            "FROM checkm as chk "
            "INNER JOIN bioentry AS entry ON entry.bioentry_id = chk.bioentry_id;"
        )
        results = self.server.adaptor.execute_and_fetchall(query)
        hsh_accession_to_checkm = {}
        for line in results:
            hsh_values = {"completeness": float(line[1]), "contamination": float(line[2])}
            hsh_accession_to_checkm[line[0]] = hsh_values
        return hsh_accession_to_checkm

    # Optional arguments: this function will return all cog counts 
    # grouped by bioentry and cog function, except one of the following argument is passed
    #  - bioentry_ids : restrict the query to a subset of bioentries
    #  - taxon_ids : restrict the query to the bioentries whose taxid matches those in this list
    def get_cog_count_for_genomes(self, bioentry_ids=None, taxon_ids=None):
        additional_query = ""

        if bioentry_ids != None and taxon_ids != None:
            raise RuntimeError("For the moment, cannot have both bioentry_ids and taxon_ids")

        if bioentry_ids != None:
            # to test
            condition = ", ".join([str(bid) for bid in bioentry_ids])
            additional_query = f" WHERE entry.bioentry_id IN ({condition}) "

        if taxon_ids != None:
            # to test
            condition = ", ".join([str(tid) for tid in taxon_ids])
            additional_query = (
                " INNER JOIN taxon AS t ON t.taxon_id = entry.taxon_id "
                f" WHERE t.ncbi_taxon_id IN ({condition})"
            )

        query = (
            "SELECT entry.bioentry_id, cog.function, MIN(hits.evalue)"
            " FROM cog_hits AS hits "
            " INNER JOIN sequence_hash_dictionnary AS seq ON seq.hsh = hits.hsh "
            " INNER JOIN cog_names AS cog ON cog.cog_id = hits.cog_id "
            " INNER JOIN seqfeature AS fet ON fet.seqfeature_id = seq.seqid "
            " INNER JOIN bioentry AS entry ON fet.bioentry_id = entry.bioentry_id "
            f"{additional_query}"
            " GROUP BY fet.seqfeature_id;"
        )

        results = self.server.adaptor.execute_and_fetchall(query)
        hsh_results = {}
        for line in results:
            bioentry = line[0]
            func = line[1]
            if bioentry in hsh_results:
                cnt = hsh_results[bioentry].get(func, 0)
                hsh_results[bioentry][func] = cnt+1
            else:
                hsh_results[bioentry] = {func: 1}
        return hsh_results

    def get_cog_summaries(self, cog_ids, only_cog_desc=False):
        ids = ",".join(["?"] * len(cog_ids))
        query = (
            "SELECT cog_id, function, description "
            "FROM cog_names "
            f"WHERE cog_id IN ({ids});"
        )
        results = self.server.adaptor.execute_and_fetchall(query, cog_ids)
        if only_cog_desc:
            hsh_results = {}
            for line in results:
                hsh_results[line[0]] = (line[1], line[2])
            return hsh_results

        funcs = "SELECT function, description FROM cog_functions;"
        functions = self.server.adaptor.execute_and_fetchall(funcs)
        hsh_func_to_description = {}
        for function, description in functions:
            hsh_func_to_description[function] = description

        hsh_results = { }
        for cog_id, function, cog_description in results:
            hsh_results[cog_id] = []
            for i in range(0, len(function)):
                func = function[i]
                func_descr = hsh_func_to_description[func]
                hsh_results[cog_id].append((func, func_descr, cog_description))
        return hsh_results

    # Note: ordering by seqid makes it faster to assemble informations
    # from several queries if the index is the same.
    def get_gene_loc(self, seqids, as_hash=True):
        seqids_query = ",".join(["?"] * len(seqids))

        query = (
            "SELECT seqfeature_id, strand, start_pos, end_pos "
            f"FROM location WHERE seqfeature_id IN ({seqids_query});"
        )

        results = self.server.adaptor.execute_and_fetchall(query, seqids)
        if not as_hash:
            return results

        hsh_results = {}
        for line in results:
            seqid = line[0]
            strand = line[1]
            start = line[2]
            end = line[3]
            hsh_results[seqid] = [strand, start, end]
        return hsh_results
    
    # Returns the seqid, locus tag, protein id, product and gene for a given
    # list of seqids, ordered by seqids
    #
    # WARNING: fam_cog relies on the result being ordered by seqfeature_id
    def get_proteins_info(self, seqids):
        seqids_query = ",".join(["?"] * len(seqids))
        term_names = ["locus_tag", "protein_id", "gene", "product"]
        term_names_query = ",".join([f"\"{name}\"" for name in term_names])
        query = (
            "SELECT seqfeature_id, name, value "
            "FROM seqfeature_qualifier_value AS v "
            "INNER JOIN term AS t ON t.term_id = v.term_id "
            f"WHERE seqfeature_id IN ({seqids_query}) AND name IN ({term_names_query}) "
            "ORDER BY seqfeature_id ASC;"
        )
        results = self.server.adaptor.execute_and_fetchall(query, seqids)

        hsh_results = {}
        for line in results:
            seqid = line[0]
            term_name = line[1]
            value = line[2]
            if seqid not in hsh_results:
                hsh_results[seqid] = [None]*len(term_names)
            hsh_results[seqid][term_names.index(term_name)] = value
        return hsh_results

    def get_organism(self, ids, as_hash=True, id_type="seqid"):
        seqids_query = ",".join(["?"] * len(ids))

        if id_type == "seqid":
            query = (
                "SELECT feature.seqfeature_id, organism.value "
                "FROM seqfeature AS feature "
                "INNER JOIN bioentry AS entry ON feature.bioentry_id = entry.bioentry_id "
                "INNER JOIN bioentry_qualifier_value AS organism "
                "    ON organism.bioentry_id = entry.bioentry_id "
                "INNER JOIN term AS organism_term ON organism.term_id = organism_term.term_id "
                " AND organism_term.name = \"organism\" "
                f"WHERE feature.seqfeature_id IN ({seqids_query});"
            )
        elif id_type == "bioentry":
            query = (
                "SELECT organism.bioentry_id, organism.value "
                "FROM bioentry_qualifier_value AS organism "
                "INNER JOIN term AS organism_term ON organism.term_id = organism_term.term_id "
                " AND organism_term.name = \"organism\" "
                f"WHERE organism.bioentry_id IN ({seqids_query});"
            )
        else:
            raise RuntimeError("Id_type should be either seqid or bioentry")

        results = self.server.adaptor.execute_and_fetchall(query, ids)
        if not as_hash:
            return results 

        hsh_results = {}
        for line in results:
            seqid = line[0]
            organism = line[1]
            hsh_results[seqid] = organism
        return hsh_results

    # Note:
    # First extracting the data in memory and creating the dataframe
    # from it is much faster than iterating over results and adding
    # elements separately
    def to_pandas_frame(db_results, columns):
        data = []
        for line in db_results:
            data.append(line)
        return pd.DataFrame(data, columns=columns)

    # For each genome, return the number of gene that were assigned
    # to each orthogroup passed in argument
    def get_og_count(self, lookup_term, search_on="orthogroup", indexing="bioentry"):
        if indexing != "bioentry" and indexing != "taxid":
            raise RuntimeError("Only bioentry and taxid indexing are supported")
        if search_on!="orthogroup" and search_on!="bioentry":
            raise RuntimeError("Can only search on orthogroup or bioentry")

        where_clause = "orthogroup"
        if search_on=="bioentry":
            where_clause = "entry.bioentry_id"
        
        entries = ",".join("?" for i in lookup_term)
        query = (
            "SELECT feature.bioentry_id, orthogroup, count(*) "
            "FROM bioentry AS entry "
            "INNER JOIN seqfeature AS feature ON entry.bioentry_id = feature.bioentry_id " 
            "INNER JOIN og_hits AS og ON og.seqid = feature.seqfeature_id "
            f"WHERE {where_clause} IN ({entries}) "
            "GROUP BY feature.bioentry_id, orthogroup;"
        )
        results = self.server.adaptor.execute_and_fetchall(query, lookup_term)
        return DB.to_pandas_frame(results, ["bioentry", "orthogroup", "count"])

    # NOTE: should be removed (and use get_cog_hits)
    # Note: only takes the best hits into account
    # Note II: code written on a Friday at 7pm... need to check it thoroughly
    #
    # For now, the query identifies the best cog hit for all seqid
    # where cogs from cog_list are located. The code below the query then only
    # count the hits that correspond to entries in cog_list.
    # This may need to be rewritten, especially if we end up choosing to not
    # include all cog_hits.
    def get_cog_counts(self, cog_list):
        entries = ",".join(["?"] * len(cog_list))

        query = (
            "SELECT feature.bioentry_id, all_hits.cog_id, MIN(all_hits.evalue) "
            " FROM sequence_hash_dictionnary AS hsh "
            " INNER JOIN seqfeature AS feature ON hsh.seqid = feature.seqfeature_id "
            " INNER JOIN cog_hits AS hit ON hsh.hsh = hit.hsh "
            " INNER JOIN sequence_hash_dictionnary AS all_hsh ON all_hsh.seqid = hsh.seqid "
            " INNER JOIN cog_hits AS all_hits ON all_hsh.hsh = all_hits.hsh "
            f"WHERE hit.cog_id IN ({entries}) "
            " GROUP BY all_hsh.seqid; "
        )
        results = self.server.adaptor.execute_and_fetchall(query, cog_list)
        hsh_results = {}
        for line in results:
            bioentry, cog_id = line[0:2]
            if not cog_id in cog_list:
                continue

            if bioentry in hsh_results:
                cnt = hsh_results[bioentry].get(cog_id, 0)
                hsh_results[bioentry][cog_id] = cnt+1
            else:
                hsh_results[bioentry] = {cog_id : 1}
        return hsh_results

    # NOTE: should be modified as selecting the cog with the best value is not
    # necessary anymore
    # For now, only returns the best hits, may be interesting to modify it later
    # to return the count for all hits or to use a threshold.
    def get_cog_counts_per_category(self, bioentries):
        entries = ",".join(["?"] * len(bioentries))
        query = (
            "SELECT bioentry_id, cog.function, MIN(hit.evalue) "
            "FROM seqfeature "
            "INNER JOIN sequence_hash_dictionnary AS hsh on hsh.seqid = seqfeature_id "
            "INNER JOIN cog_hits AS hit ON hit.hsh = hsh.hsh "
            "INNER JOIN cog_names AS cog ON cog.cog_id = hit.cog_id "
            f"WHERE bioentry_id IN ({entries}) "
            "GROUP BY bioentry_id, seqfeature_id;"
        )
        results = self.server.adaptor.execute_and_fetchall(query, bioentries)
        hsh_results = {}
        for line in results:
            entry_id, func, foo = line[0:3]
            if entry_id not in hsh_results:
                hsh_results[entry_id] = {func : 1}
            else:
                cnt = hsh_results[entry_id].get(func, 0)
                hsh_results[entry_id][func] = cnt+1
        return hsh_results

    def get_og(self, seqids, order=True):
        entries = ",".join("?" for i in seqids)
        ordering = ""
        if order:
            ordering = "ORDER BY seqid ASC"
        query = (
            "SELECT seqid, orthogroup "
            f"FROM og_hits WHERE seqid IN ({entries})"
            f"{ordering};"
        )
        results = self.server.adaptor.execute_and_fetchall(query, seqids)
        hsh_results = {}
        for line in results:
            hsh_results[line[0]] = line[1]
        return hsh_results

    # Get all cog hits for a given list of bioentries
    # The results are either indexed by the bioentry or by the seqid
    def get_cog_hits(self, bioentries, index_by_seqid=False, as_count=False):
        if index_by_seqid and as_count:
            raise RuntimeError("Counting cog per seqid is useless (only best hits are kept)")

        entries = ",".join("?" for i in bioentries)

        index = "entry.bioentry_id"
        if index_by_seqid:
            index = "seqid.seqfeature_id"

        prefix, suffix = "", ""
        if as_count:
            prefix = ", COUNT(*)"
            suffix = f"GROUP BY entry.bioentry_id, cogs.cog_id "

        query = (
            f"SELECT {index}, cogs.cog_id {prefix}"
            "FROM bioentry AS entry "
            "INNER JOIN seqfeature AS seqid ON seqid.bioentry_id = entry.bioentry_id "
            "INNER JOIN sequence_hash_dictionnary AS hsh ON seqid.seqfeature_id = hsh.seqid "
            "INNER JOIN cog_hits AS cogs ON cogs.hsh = hsh.hsh "
            f"WHERE entry.bioentry_id IN ({entries}) " 
            f"{suffix};"
        )
        results = self.server.adaptor.execute_and_fetchall(query, bioentries)
        if as_count:
            return DB.to_pandas_frame(results, ["bioentry", "cog", "count"])

        hsh_results = {}
        # legacy code, to be replaced in profit of pandas dataframe
        for line in results:
            entry = line[0]
            cog = line[1]
            if entry in hsh_results:
                hsh_results[entry].append(cog)
            else:
                hsh_results[entry] = [cog]
        return hsh_results

    # May need to modify this so it can also index results by bioentry
    def get_all_seqfeature_for_cog(self, cogs):
        query_str = ",".join("?" for i in cogs)
        query = (
            "SELECT hsh.seqid, hit.cog_id "
            "FROM cog_hits AS hit "
            "INNER JOIN sequence_hash_dictionnary AS hsh ON hsh.hsh = hit.hsh "
            f"WHERE hit.cog_id IN ({query_str}); "
        )
        results = self.server.adaptor.execute_and_fetchall(query, cogs)

        hsh_results = {}
        for line in results:
            seqid, cog_id = line
            if cog_id not in cogs:
                continue
            assert seqid not in hsh_results
            hsh_results[seqid] = cog_id
        return hsh_results


    def get_filenames_to_bioentry(self):
        sql = (
            "SELECT * FROM filenames;"
        )
        hsh_filenames_to_bioentry = {}
        results = self.server.adaptor.execute_and_fetchall(sql)
        for line in results:
            hsh_filenames_to_bioentry[line[1].replace(".gbk", "")] = line[0]
        return hsh_filenames_to_bioentry


    def load_chlamdb_config_tables(self, entries):
        sql = (
            "CREATE TABLE biodb_config"
            "(name varchar(200), type varchar(200), status BOOLEAN);"
        )
        self.server.adaptor.execute(sql)
        self.load_data_into_table("biodb_config", entries)


    # wrapper methods
    def commit(self):
        self.server.commit()

    def load_gbk_wrapper(self, records):
        self.server[self.db_name].load(records)

    # Maybe return different instance of a subclass depending on the type
    # of database? Would allow to avoid code duplication if several database
    # types are to be included.
    def load_db(db_file, params):
        sqlpsw = params.get("chlamdb.db_psswd", "")
        db_type = params["chlamdb.db_type"]
        db_name = params["chlamdb.db_name"]

        if db_type != "sqlite":
            server = BioSeqDatabase.open_database(driver="MySQLdb", 
                                                  user="root",
                                                  passwd = sqlpsw, 
                                                  host = "127.0.0.1", 
                                                  db=db_file, 
                                                  charset='utf8',
                                                  use_unicode=True)
        else:
            server = BioSeqDatabase.open_database(driver="sqlite3", 
                                                  user="root",
                                                  passwd = sqlpsw, 
                                                  host = "127.0.0.1",
                                                  db=db_file)
        return DB(server, db_name)

    def load_db_from_name(db_name, db_type = "sqlite"):
        params = {"chlamdb.db_type" : db_type, "chlamdb.db_name" : db_name}
        return DB.load_db(params)
