import os
import sys

from BioSQL import BioSeqDatabase
from Bio import SeqIO
from Bio.SeqUtils import GC

import sqlite3


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

    def add_orthogroups_to_seq(self, hsh_locus_to_group, hsh_locus_to_feature_id, orthogroup_term_id):
        rank = 1
        for locus, group_id in hsh_locus_to_group.items():
            feature_id = hsh_locus_to_feature_id[locus]
            insert = (
                f"INSERT INTO seqfeature_qualifier_value "
                f"VALUES ({feature_id}, {orthogroup_term_id}, {rank}, {quote(group_id)})"
            )
            self.server.adaptor.execute(insert)

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

    # Returns a hash that maps the accesion to the orthogroup it was linked to
    # Note: this function already filters the top n hits
    def get_non_PVC_refseq_matches(self, params):
        pvc = params["refseq_diamond_BBH_phylogeny_phylum_filter"]
        max_hit_count = params["refseq_diamond_BBH_phylogeny_top_n_hits"]
        pvc_fmt_string = ",".join([f"\"{phylum}\"" for phylum in pvc])
        query = (
            "SELECT feature.seqfeature_id, refseq_hit_id.accession, feature.value, refseq_hit.hit_count "
            "FROM diamond_refseq AS refseq_hit "
            "INNER JOIN diamond_refseq_match_id AS refseq_hit_id "
            "   ON refseq_hit.sseqid=refseq_hit_id.match_id "
            " INNER JOIN sequence_hash_dictionnary AS hsh_dict ON hsh_dict.hsh = refseq_hit.seq_hash"
            " INNER JOIN seqfeature_qualifier_value as feature "
            "   ON feature.seqfeature_id=hsh_dict.seqid "
            "INNER JOIN term as t ON t.term_id = feature.term_id AND t.name=\"orthogroup\" "
            "INNER JOIN refseq_hits_taxonomy AS taxo ON taxo.taxid=refseq_hit_id.taxid "
            "INNER JOIN taxonomy_mapping taxo_name "
            f"   ON taxo.phylum=taxo_name.taxid AND taxo_name.value NOT IN ({pvc_fmt_string}) "
            "ORDER BY feature.seqfeature_id ASC, hit_count ASC;"
        )
        results = self.server.adaptor.execute_and_fetchall(query, )

        # The following code relies on the results being grouped by qseqid
        # and ordered by hit count by the SQL query
        hsh_results = {}
        curr_tab = []
        curr_qseq_id = None
        curr_count = 0

        # The accession not in curr_tab are not really efficient, particularly
        # since they compare strings and perform in O(n). May need to improve
        # on it lots of sequences are added for each orthogroup.
        for result in results:
            qseqid, accession, orthogroup = result[0], result[1], result[2]
            if qseqid != curr_qseq_id:
                curr_tab = hsh_results.get(orthogroup, [])
                hsh_results[orthogroup] = curr_tab
                if accession not in curr_tab:
                    curr_tab.append(accession)
                curr_count = 1
                curr_qseq_id = qseqid
            elif curr_count < max_hit_count:
                if accession not in curr_tab:
                    curr_tab.append(accession)
                curr_count += 1
        return hsh_results

    def get_all_sequences_for_orthogroup(self, orthogroup):
        query = (
            "SELECT locus.value, seq.value "
            "FROM seqfeature_qualifier_value AS ortho "
            "INNER JOIN term as ortho_term "
            "   ON ortho_term.term_id=ortho.term_id AND ortho_term.name=\"orthogroup\" "
            "INNER JOIN seqfeature_qualifier_value AS seq ON seq.seqfeature_id=ortho.seqfeature_id "
            " INNER JOIN term AS seq_term "
            "   ON seq_term.term_id=seq.term_id AND seq_term.name=\"translation\""
            "INNER JOIN seqfeature_qualifier_value locus ON locus.seqfeature_id=ortho.seqfeature_id "
            "INNER JOIN term as locus_term "
            "   ON locus_term.term_id=locus.term_id AND locus_term.name=\"locus_tag\" "
            f"WHERE ortho.value={orthogroup};"
        )
        results = self.server.adaptor.execute_and_fetchall(query,)
        tab = []
        for result in results:
            tab.append([result[0], result[1]])
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
        fmt_string = ", ".join(["?"] * len(data[0]))
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

    hsh_taxo_key = {
            "superkingdom" : (1, 2),
            "phylum" :       (3, 4),
            "class" :        (5, 6),
            "order" :        (7, 8),
            "family" :       (9, 10),
            "genus" :        (11, 12),
            "species" :      (13, 14)
    }
    def parse_taxo(self, results):
        # taxid
        lst = [ results[0] ]
        for (idx_name, idx_taxid) in self.hsh_taxo_key.values():
            lst.append( results[idx_name] )
            lst.append( results[idx_taxid] )
        return lst

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
        lst_results = []
        for line in results:
            lst_results.append(self.parse_taxo(line))
        return lst_results

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

    def get_accession_to_prot(self, accession, params):
        if self.conn_refseq == None:
            self.conn_refseq = sqlite3.connect(params["databases_dir"] + "/refseq/merged_refseq.db")
        cursor = self.conn_refseq.cursor()

        fmt_string = ",".join(["\"" + a + "\"" for a in accession])
        query = (
            f"SELECT accession, description, sequence_length FROM refseq "
            f"WHERE accession IN ({fmt_string});"
        )
        results = cursor.execute(query,).fetchall()
        hsh_results = {}
        for line in results:
            hsh_results[line[0]] = (line[1], int(line[2]))
        return hsh_results

    # NOTE: need to check which indices are necessary to add to this table
    def load_cog_hits(self, data):
        sql = (
            "CREATE TABLE cog_hits (hsh INTEGER, "
            " cog_id INT, query_start INT, query_end INT, hit_start INT, "
            " hit_end INT, query_coverage FLOAT, hit_coverage FLOAT, identity FLOAT, "
            " evalue FLOAT, bitscore FLOAT, "
            " FOREIGN KEY(cog_id) REFERENCES cog_names(cog_id)); "
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
    
    def load_cog_ref_data(self, data):
        sql = (
            "CREATE TABLE cog_names (cog_id INTEGER, function TEXT, description TEXT, "
            "PRIMARY KEY(cog_id));"
        )
        self.server.adaptor.execute(sql)
        self.load_data_into_table("cog_names", data)

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
            f"FOREIGN KEY(orthogroup) REFERENCES orthology_orthogroup(orthogroup_id)"
            f"FOREIGN KEY(id_1) REFERENCES seqfeature(seqfeature_id)"
            f"FOREIGN KEY(id_2) REFERENCES seqfeature(seqfeature_id))"
        )
        self.server.adaptor.execute(query, )

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

    def load_og_averages(self, averages):
        for og, average in averages:
            sql = (
                f"UPDATE orthology_orthogroup "
                f"SET average_identity = {average} "
                f"WHERE orthogroup_id = {og};"
            )
            self.server.adaptor.execute(sql,)

    def create_orthology_table(self, arr_cnt_tables):
        sql = (
            f"CREATE TABLE IF NOT EXISTS orthology_orthogroup( "
            f"orthogroup_id INTEGER, orthogroup_size INT, n_genomes INT, average_identity FLOAT(5), "
            f"PRIMARY KEY(orthogroup_id)) "
        )
        self.server.adaptor.execute(sql,)

        for group, hsh_cnt_tables in enumerate(arr_cnt_tables):
            size = 0
            n_genomes = 0
            for taxon_id, cnt in hsh_cnt_tables.items():
                size += cnt
                if cnt>0:
                    n_genomes += 1
            sql = (
                f"INSERT INTO orthology_orthogroup VALUES("
                f"{group}, {size}, {n_genomes}, NULL"
                f")"
            )
            self.server.adaptor.execute(sql,)

    def create_BBH_phylogeny_table(self, data):
        sql = (
            "CREATE TABLE BBH_phylogeny (orthogroup_id INTEGER, tree TEXT, "
            " PRIMARY KEY(orthogroup_id), "
            " FOREIGN KEY(orthogroup_id) REFERENCES orthology_orthogroup(orthogroup_id));"
        )
        self.server.adaptor.execute(sql,)
        self.load_data_into_table("BBH_phylogeny", data)

    def create_gene_phylogeny_table(self, data):
        sql = (
            "CREATE TABLE gene_phylogeny (orthogroup_id INTEGER, tree TEXT, "
            "PRIMARY KEY(orthogroup_id), "
            " FOREIGN KEY(orthogroup_id) REFERENCES orthology_orthogroup(orthogroup_id));"
        )
        self.server.adaptor.execute(sql,)
        self.load_data_into_table("gene_phylogeny", data)


    # NOTE: may be more efficient to create indices on a combination of keys, depending
    # on how the table is used.
    def create_og_matrix_indices(self):
        sql_1 = "CREATE INDEX oio ON orthology_identity(orthogroup);"
        sql_2 = "CREATE INDEX oiid_1 ON orthology_identity(id_1);"
        sql_3 = "CREATE INDEX oiid_2 ON orthology_identity(id_2);"
        self.server.adaptor.execute(sql_1,)
        self.server.adaptor.execute(sql_2,)
        self.server.adaptor.execute(sql_3,)

    def load_orthology_table(self, arr_count_tables):
        create_table_sql = (
            f"CREATE TABLE comparative_tables_orthology (orthogroup INT, "
            f" taxid INT, count INT"
            f", PRIMARY KEY(orthogroup, taxid)"
            f", FOREIGN KEY(orthogroup) REFERENCES orthology_orthogroup(orthogroup_id)"
            f", FOREIGN KEY(taxid) REFERENCES taxon(taxon_id)"
            f")"
        )
        self.server.adaptor.execute(create_table_sql)

        for group, hsh_taxid_to_count in enumerate(arr_count_tables):
            for taxid, count in hsh_taxid_to_count.items():
                sql = (
                    f"INSERT INTO comparative_tables_orthology "
                    f"VALUES ({group}, {taxid}, {count})" 
                )
                self.server.adaptor.execute(sql,)

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
            f"PRIMARY KEY (hit_count, seq_hash));"
        )
        self.server.adaptor.execute(sql,)

    def hsh_to_prot_length(self):
        query = (
            "SELECT DISTINCT hsh, length(value) "
            "FROM seqfeature_qualifier_value AS value "
            "INNER JOIN term ON term.term_id = value.term_id AND term.name=\"translation\" "
            "INNER JOIN sequence_hash_dictionnary AS hashes ON hashes.seqid = value.seqfeature_id;"
        )
        results = self.server.adaptor.execute_and_fetchall(query,)
        hsh = {}
        for line in results:
            hsh_id = line[0]
            length = line[1]
            hsh[hsh_id] = length
        return hsh

    def load_refseq_hits(self, data):
        self.load_data_into_table("diamond_refseq", data)

    def setup_orthology_table(self):
        sql1 = 'SELECT ontology_id FROM ontology WHERE name="SeqFeature Keys"'
        ontology_id = self.server.adaptor.execute_and_fetchall(sql1)[0][0]
        sql2 = f"INSERT INTO term (name, ontology_id) VALUES (\"orthogroup\", {ontology_id});"
        self.server.adaptor.execute(sql2)
        return self.server.adaptor.cursor.lastrowid

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

    def get_genomes_description(self, indexing="bioentry", exclude_plasmids=False, indexing_type="str"):
        if indexing != "bioentry" and indexing != "accession":
            raise RuntimeError(f"{bioentry} indexing is not supported")

        if indexing_type != "int" and indexing_type != "str":
            raise RuntimeError(f"{indexing_type} not supported, must be int or str")

        selection = ("accession", "bioentry_id")[indexing=="bioentry"]
        query = (
            f"SELECT {selection}, description "
            " FROM bioentry "
            " ORDER BY description;"
        )
        results = self.server.adaptor.execute_and_fetchall(query)
        hsh_results = {}
        for line in results:
            if exclude_plasmids and ("plasmid" in line[1] or "phage" in line[1]):
                continue
            idx = line[0] if indexing_type == "int" else str(line[0])
            hsh_results[idx] = line[1]
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
            gc, length = float(other_entry[1]), int(other_entry[2])
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

        # This needs some testing
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
            index = "bioentry_id"
        elif indexing == "accession":
            index = "accession"
        elif indexing == "ncbi":
            index = "ncbi_taxon_id"
        else:
            raise RuntimeError(f"Unsupported indexing method: {indexing}")

        query = (
            f"SELECT {indexing}, COUNT(*) FROM seqfeature AS prot "
            " INNER JOIN term AS t ON t.term_id = prot.type_term_id "
            " INNER JOIN bioentry AS entry ON entry.bioentry_id = prot.bioentry_id "
            " AND t.name=\"CDS\" GROUP BY accession;"
            f"{filtering}"
        )
        results = self.server.adaptor.execute_and_fetchall(query)
        hsh_results = {}
        for accession, str_count in results:
            hsh_results[accession] = int(str_count)
        return hsh_results

    def get_cog_code_description(self):
        query = (
            "SELECT function, description FROM cog_functions;"
        )
        results = self.server.adaptor.execute_and_fetchall(query)
        hsh_results = {}
        for line in results:
            hsh_results[line[0]] = line[1]
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
    def get_COG_counts(self, bioentry_ids=None, taxon_ids=None):
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
            "SELECT entry.bioentry_id, cog.function, COUNT(*)"
            " FROM cog_hits AS hits "
            " INNER JOIN sequence_hash_dictionnary AS seq ON seq.hsh = hits.hsh "
            " INNER JOIN cog_names AS cog ON cog.cog_id = hits.cog_id "
            " INNER JOIN seqfeature AS fet ON fet.seqfeature_id = seq.seqid "
            " INNER JOIN bioentry AS entry ON fet.bioentry_id = entry.bioentry_id "
            f"{additional_query}"
            " GROUP BY entry.bioentry_id, cog.function;"
        )

        results = self.server.adaptor.execute_and_fetchall(query)
        hsh_results = {}
        for line in results:
            bioentry = line[0]
            func = line[1]
            count = line[2]
            if bioentry in hsh_results:
                hsh_results[bioentry][func] = count
            else:
                hsh_results[bioentry] = {func: count}
        return hsh_results

    def get_n_prot_without_cog(self):
        query = (
            "SELECT entry.bioentry_id, COUNT(*) "
            " FROM bioentry AS entry "
            " INNER JOIN seqfeature AS feature ON entry.bioentry_id = feature.bioentry_id "
            " INNER JOIN sequence_hash_dictionnary AS dict ON dict.seqid = feature.seqfeature_id"
            " WHERE NOT EXISTS ( "
            "   SELECT 1 FROM cog_hits WHERE hsh = dict.hsh  "
            " ) "
            " GROUP BY entry.bioentry_id; "
        )

        results = self.server.adaptor.execute_and_fetchall(query)
        hsh_results = {}
        for line in results:
            hsh_results[line[0]] = results[1]
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
            return results

        funcs = "SELECT function, description FROM cog_functions;"
        functions = self.server.adaptor.execute_and_fetchall(funcs)
        hsh_func_to_description = {}
        for function, description in functions:
            hsh_func_to_description[function] = description

        all_results = []
        for cog_id, function, cog_description in results:
            for i in range(0, len(function)):
                func = function[i]
                func_descr = hsh_func_to_description[func]
                all_results.append( (cog_id, func, func_descr, cog_description) )
        return all_results

    # Note: ordering by seqid makes it faster to assemble informations
    # from several queries if the index is the same.
    def get_gene_loc(self, seqids, order=True):
        seqids_query = ",".join(["?"] * len(seqids))
        ordering = ""
        if order:
            ordering = "ORDER BY seqfeature_id ASC"

        query = (
            "SELECT seqfeature_id, strand, start_pos, end_pos "
            f"FROM location WHERE seqfeature_id IN ({seqids_query})"
            f"{ordering};"
        )
        return self.server.adaptor.execute_and_fetchall(query, seqids)
    
    # Returns the locus tag, protein id of the protein from seqids
    def get_proteins_info(self, seqids, order=True):
        seqids_query = ",".join(["?"] * len(seqids), order=True)

        ordering = ""
        if order:
            ordering = "ORDER BY seqfeature_id ASC"

        query = (
            "SELECT seqfeature_id, locus_tag.value, prot_id.value, product.value, gene.value "
            "FROM seqfeature_qualifier_value AS locus_tag "
            "INNER JOIN seqfeature_qualifier_value AS prot_id "
            "  ON locus_tag.seqfeature_id = prot_id.seqfeature_id "
            " INNER JOIN term AS prot_id_term ON prot_id_term.term_id = prot_id.term "
            "  AND prot_id_term.name = \"protein_id\" "
            "INNER JOIN seqfeature_qualifier_value AS product " 
            "  ON product.seqfeature_id=locus_tag.seqfeature_id "
            " INNER JOIN term AS product_term ON product_term.term_id = product.term_id "
            "  AND product_term.name = \"product\" "
            "INNER JOIN seqfeature_qualifier_value AS gene "
            "  ON gene.seqfeature_id=locus_tag.seqfeature_id "
            " INNER JOIN term AS gene_term ON gene_term.term_id=gene.term_id "
            "  AND gene_term.name = \"gene\" "
            f"WHERE locus_tag.seqfeature_id IN ({seqids_query})"
            f"{ordering};"
        )
        return self.server.adaptor.execute_and_fetchall(query, seqids)

    def get_organism(self, seqids, order=True):
        seqids_query = ",".join(["?"] * len(seqids))
        ordering = ""
        if order:
            ordering = "ORDER BY seqfeature_id ASC"

        query = (
            "SELECT feature.seqfeature_id, organism.value "
            "FROM seqfeature AS feature "
            "INNER JOIN bioentry AS entry ON feature.bioentry_id = entry.bioentry_id "
            "INNER JOIN bioentry_qualifier_value AS organism ON organism.bientry_id = entry.bioentry_id "
            "INNER JOIN term AS organism_term ON organism.term_id = organism_term.term_id "
            " AND organism_term.name = \"organism\" "
            f"WHERE feature.seqfeature_id IN ({seqids_query})"
            f"{ordering};"
        )

        return self.server.adaptor.execute_and_fetchall(query, seqids)

    def get_og(self, seqids, order=True):
        entries = ",".join(["?"] * len(seqids))
        ordering = ""
        if order:
            ordering = "ORDER BY seqfeature_id ASC"
        query = (
            "SELECT seqfeature_id, value "
            "FROM seqfeature_qualifier_value "
            "INNER JOIN term ON term_id AND name =\"orthogroup\" "
            f"WHERE seqfeature_id IN ({entries}) "
            f"{ordering};"
        )
        results = self.server.adaptor.execute_and_fetchall(query, seqids)
        hsh_results = {}
        for line in results:
            hsh_results[line[0]] = line[1]
        return hsh_results

    def get_cog_hits(self, bioentries):
        entries = ",".join([str(entry) for entry in bioentries])
        query = (
            "SELECT entry.bioentry_id, cogs.cog_id "
            "FROM bioentry AS entry "
            "INNER JOIN seqfeature AS seqid ON seqid.bioentry_id = entry.bioentry_id "
            "INNER JOIN sequence_hash_dictionnary AS hsh ON seqid.seqfeature_id = hsh.seqid "
            "INNER JOIN cog_hits AS cogs ON cogs.hsh = hsh.hsh "
            f"WHERE entry.bioentry_id IN ({entries}) " 
            "ORDER BY entry.bioentry_id ASC;"
        )
        results = self.server.adaptor.execute_and_fetchall(query)
        hsh_results = {}
        for line in results:
            entry = line[0]
            cog = line[1]
            if entry in hsh_results:
                hsh_results[entry].append(cog)
            else:
                hsh_results[entry] = [cog]
        return hsh_results

    def get_all_seqfeature_for_cog(self, cog):
        query = (
            "SELECT hsh.seqid "
            "FROM cog_hits AS hit "
            "INNER JOIN sequence_hash_dictionnary AS hsh ON hsh.hsh = hit.hsh "
            f"WHERE hit.cog_id = ?;"
        )
        results = self.server.adaptor.execute_and_fetchall(query, cog)
        return results

    def get_filenames_to_bioentry(self):
        sql = (
            "SELECT * FROM filenames;"
        )
        hsh_filenames_to_bioentry = {}
        results = self.server.adaptor.execute_and_fetchall(sql)
        for line in results:
            hsh_filenames_to_bioentry[line[1].replace(".gbk", "")] = line[0]
        return hsh_filenames_to_bioentry

    # wrapper methods
    def commit(self):
        self.server.commit()

    def load_gbk_wrapper(self, records):
        self.server[self.db_name].load(records)

    # Maybe return different instance of a subclass depending on the type
    # of database? Would allow to avoid code duplication if several database
    # types are to be included.
    def load_db(params):
        sqlpsw = os.environ['SQLPSW']
        db_type = params["chlamdb.db_type"]
        db_name = params["chlamdb.db_name"]

        if db_type != "sqlite":
            server = BioSeqDatabase.open_database(driver="MySQLdb", 
                                                  user="root",
                                                  passwd = sqlpsw, 
                                                  host = "127.0.0.1", 
                                                  db=db_name, 
                                                  charset='utf8',
                                                  use_unicode=True)
        else:
            server = BioSeqDatabase.open_database(driver="sqlite3", 
                                                  user="root",
                                                  passwd = sqlpsw, 
                                                  host = "127.0.0.1",
                                                  db=f"{db_name}")
        return DB(server, db_name)

    def load_db_from_name(db_name, db_type = "sqlite"):
        params = {"chlamdb.db_type" : db_type, "chlamdb.db_name" : db_name}
        return DB.load_db(params)
