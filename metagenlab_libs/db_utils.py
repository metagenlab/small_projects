import os
import sys

from BioSQL import BioSeqDatabase
from Bio import SeqIO

from Bio.Seq import Seq
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
        # this will need to be changed in case a MySQL database is used
        self.placeholder = "?"

    # the next two methods are necessary for DB objects to be used
    # in 'with' blocks.
    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.server.close()

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


    def get_taxid_from_seqid(self, seqids):
        query = ",".join("?" for _ in seqids)
        sql = (
            "SELECT fet.seqfeature_id, entry.taxon_id "
            "FROM seqfeature AS fet "
            "INNER JOIN bioentry AS entry ON entry.bioentry_id=fet.bioentry_id "
            f"WHERE fet.seqfeature_id IN ({query});"
        )
        results = self.server.adaptor.execute_and_fetchall(sql, seqids)
        return {seqid: taxon_id for seqid, taxon_id in results}


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
        hsh_results = {}
        for line in results:
            hsh_results[line[0]] = line[1]
        return hsh_results

    
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
            "module_id INTEGER, desc TEXT, definition TEXT, is_signature_module BOOL, class INT, subclass INT, "
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
            "ko_id INTEGER, descr TEXT, PRIMARY KEY(ko_id));"
        )
        self.server.adaptor.execute(sql)
        self.load_data_into_table("ko_def", data)
        sql = "CREATE INDEX kdko_i ON ko_def(ko_id);"
        self.server.adaptor.execute(sql)


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


    def get_all_modules_definition(self, allow_signature=False):
        where = ""
        if not allow_signature:
            where = "WHERE is_signature_module = 0"

        query = (
            "SELECT module_id, definition FROM ko_module_def "
            f"{where};"
        )
        results = self.server.adaptor.execute_and_fetchall(query)
        return [(line[0], line[1]) for line in results]


    def get_module_categories(self, module_ids=None):
        if module_ids!=None:
            selection_str = ",".join(str(mod_id) for mod_id in module_ids)
            selection = f"AND ko_module_def.module_id IN ({selection_str})"
        else:
            selection = ""

        query = (
            "SELECT class_id, descr "
            "FROM ko_module_def "
            "INNER JOIN ko_class ON class_id = class "
            f"WHERE is_signature_module = 0 {selection}"
            "GROUP BY (descr);"
        )
        results = self.server.adaptor.execute_and_fetchall(query)
        return [(line[0], line[1]) for line in results]


    def get_module_sub_categories(self, module_ids=None):
        if module_ids!=None:
            selection_str = ",".join(str(mod_id) for mod_id in module_ids)
            selection = f"AND ko_module_def.module_id IN ({selection_str})"
        else:
            selection = ""

        query = (
            "SELECT class_id, descr "
            "FROM ko_module_def "
            "INNER JOIN ko_class ON class_id = subclass "
            f"WHERE is_signature_module = 0 {selection}"
            "GROUP BY (descr);"
        )
        results = self.server.adaptor.execute_and_fetchall(query)
        return [(line[0], line[1]) for line in results]


    # NOTE: when new dtb, to modify so that the input 
    # is sanitized.
    def get_ko_desc(self, ko_ids):
        entries = ",".join(f"{ko_id}" for ko_id in ko_ids)
        query = (
            "SELECT ko.ko_id, ko.descr "
            "FROM ko_def as ko "
            f"WHERE ko.ko_id IN ({entries});"
        )
        results = self.server.adaptor.execute_and_fetchall(query)
        hsh_results = {}
        for line in results:
            hsh_results[line[0]] = line[1]
        return hsh_results


    def get_ko_pathways(self, ko_ids):
        entries = ",".join("?" for i in ko_ids)
        query = (
            "SELECT ktp.ko_id, path.pathway_id, path.desc "
            "FROM ko_to_pathway AS ktp "
            "INNER JOIN ko_pathway_def AS path ON path.pathway_id = ktp.pathway_id "
            f"WHERE ktp.ko_id IN ({entries});"
        )
        results = self.server.adaptor.execute_and_fetchall(query, ko_ids)
        hsh_results = {}
        for line in results:
            ko_id = line[0]
            data = hsh_results.get(ko_id, [])
            data.append((line[1], line[2]))
            hsh_results[ko_id] = data
        return hsh_results


    def get_module_kos(self, module_id):
        query = (
            "SELECT ktm.ko_id "
            "FROM ko_to_module AS ktm "
            f"WHERE ktm.module_id = {module_id};"
        )
        results = self.server.adaptor.execute_and_fetchall(query)
        lst_results = []
        for line in results:
            lst_results.append(line[0])
        return lst_results


    # compact: do not return the module description if true
    def get_ko_modules(self, ko_ids, as_pandas=False, compact=False):
        entries = ",".join("?" for i in ko_ids)
        if compact:
            supp_query = ""
            ids = ["ko_id", "module_id"]
        else:
            supp_query = ", mod.desc"
            ids = ["ko_id", "module_id", "desc"]

        query = (
            f"SELECT ktm.ko_id, mod.module_id {supp_query} "
            "FROM ko_to_module AS ktm "
            "INNER JOIN ko_module_def AS mod ON mod.module_id = ktm.module_id "
            f"WHERE ktm.ko_id IN ({entries});"
        )
        results = self.server.adaptor.execute_and_fetchall(query, ko_ids)
        hsh_results = {}
        if as_pandas:
            return DB.to_pandas_frame(results, ids)
        for line in results:
            ko_id = line[0]
            data = hsh_results.get(ko_id, [])

            if compact:
                data.append(line[1])
            else:
                data.append((line[1], line[2]))
            hsh_results[ko_id] = data
        return hsh_results


    def get_ko_count_cat(self, category=None, taxon_ids=None, category_name=None, index=True):
        if category!=None and category_name!=None:
            raise RuntimeError("Selection on both category and category name not supported")
        if category==None and category_name==None:
            raise RuntimeError("Need at least category or category name")
        sel = ""

        args = []
        if taxon_ids != None:
            sel_str = ",".join("?" for _ in taxon_ids)
            sel = f" AND entry.taxon_id IN ({sel_str})"
            args = taxon_ids

        if category != None:
            where = f"WHERE module.subclass = ? AND is_signature_module = 0 {sel}"
            args = [category] + args
        if category_name != None:
            where = (
                "INNER JOIN ko_class AS class ON module.subclass = class.class_id "
                f"WHERE class.descr = ? AND is_signature_module = 0 {sel}"
            )
            args = [category_name] + args

        query = (
            "SELECT entry.taxon_id, module.module_id, ktm.ko_id, COUNT(*) "
            "FROM ko_module_def AS module "
            "INNER JOIN ko_to_module AS ktm ON module.module_id = ktm.module_id "
            "INNER JOIN ko_hits AS hit ON hit.ko_id = ktm.ko_id "
            "INNER JOIN sequence_hash_dictionnary AS hsh ON hit.hsh = hsh.hsh "
            "INNER JOIN seqfeature AS fet ON fet.seqfeature_id = hsh.seqid "
            "INNER JOIN bioentry AS entry ON fet.bioentry_id=entry.bioentry_id "
            f"{where}"
            "GROUP BY entry.taxon_id, module.module_id, ktm.ko_id;"
        )
        results = self.server.adaptor.execute_and_fetchall(query, args)
        columns = ["taxon_id", "module_id", "KO", "count"]
        df = DB.to_pandas_frame(results, columns)
        if index == False:
            return df
        return df.set_index(["taxon_id", "module_id", "KO"])


    def get_modules_info(self, modules_id, as_pandas=False):
        fmt = ",".join("?" for i in modules_id)
        query = (
            "SELECT module_id, desc, definition, cat.descr, subcat.descr "
            "FROM ko_module_def AS def "
            "INNER JOIN ko_class AS subcat ON subcat.class_id = def.subclass "
            "INNER JOIN ko_class AS cat ON cat.class_id = def.class "
            f"WHERE is_signature_module = 0 AND module_id IN ({fmt});"
        )
        results = self.server.adaptor.execute_and_fetchall(query, modules_id)

        if as_pandas:
            return DB.to_pandas_frame(results, ["module_id", "descr", "definition", "cat", "subcat"])
        return [(line[0], line[1], line[2], line[3], line[4]) for line in results]


    def get_ko_count_for_ko(self, ko_id):
        has_multiple = isinstance(ko_id, list)
        if has_multiple:
            search_str = ",".join(str(ko) for ko in ko_id)
            selection_query = f"IN ({search_str})"
        else:
            selection_query = f" = {ko_id}"

        query = (
            "SELECT fet.bioentry_id, hits.ko_id, COUNT(*) "
            "FROM ko_hits as hits "
            "INNER JOIN sequence_hash_dictionnary AS hsh ON hsh.hsh = hits.hsh "
            "INNER JOIN seqfeature AS fet ON fet.seqfeature_id = hsh.seqid "
            f"WHERE hits.ko_id {selection_query} "
            "GROUP BY fet.bioentry_id, hits.ko_id;"
        )
        results = self.server.adaptor.execute_and_fetchall(query)

        if has_multiple:
            return DB.to_pandas_frame(results, ["bioentry", "ko_id", "count"])

        hsh_results = {}
        for line in results:
            bioentry, ko_id, cnt = line[0], line[1], line[2]
            hsh_results[bioentry] = cnt
        return hsh_results


    # NOTE:
    # All those get_*_hits could be improved in several ways:
    #  - probably code refactoring: since the code is quite similar
    #    it may be possible to factor out the redundant code
    #  - when querying for seqids, two many joins are being performed
    #    it may be worth it to simplify this
    def get_ko_hits(self, ids, search_on="taxid", keep_taxid=False):
        """
        Note: if search_on = ko, only the seqid are returned by default, if the keep_taxids
        is set, taxid are returned in an additional column.
        """

        prot_query = self.gen_placeholder_string(ids)

        header      = None
        select      = "SELECT seqid.seqfeature_id, hit.ko_id "
        group_by    = ""
        search_term = ""
        if search_on == "seqid":
            search_term = "seqid.seqfeature_id"
            header      = ["seqid", "ko"]
        elif search_on == "taxid":
            select      = "SELECT entry.taxon_id, hit.ko_id, COUNT(*) "
            search_term = "entry.taxon_id"
            group_by    = "GROUP BY entry.taxon_id, hit.ko_id "
            header      = ["taxid", "ko", "count"]
        elif search_on == "ko":
            search_term = "hit.ko_id"
            header      = ["seqid", "ko"]
        else:
            raise RuntimeError(f"Searching on {search_on} not supported")

        if keep_taxid and not search_on=="taxid":
            select += ", entry.taxon_id "
            header.append("taxid")
        elif keep_taxid:
            raise RuntimeError(("Are you mocking me? Taxid is already returned! "
                "Please remove this pesky keep_taxid flag or use another search term!"))

        query = (
            f"{select}"
            "FROM bioentry AS entry "
            "INNER JOIN seqfeature AS seqid ON seqid.bioentry_id=entry.bioentry_id "
            "INNER JOIN sequence_hash_dictionnary AS hsh ON hsh.seqid = seqid.seqfeature_id "
            "INNER JOIN ko_hits AS hit ON hit.hsh=hsh.hsh "
            f"WHERE {search_term} IN ({prot_query}) "
            f"{group_by};"
        )
        results = self.server.adaptor.execute_and_fetchall(query, ids)

        df = DB.to_pandas_frame(results, header)
        if df.empty:
            return df

        if search_on=="taxid":
            df = df.set_index(["taxid", "ko"]).unstack(level=0, fill_value=0)
            df.columns = [col for col in df["count"].columns.values]
        elif search_on=="seqid" or search_on=="ko":
            df = df.set_index(["seqid"])
        return df


    def get_ko_count(self, search_entries, keep_seqids=False, search_on="taxid", as_multi=True):
        if search_on=="ko_id":
            where_clause = "hit.ko_id"
        elif search_on=="taxid":
            where_clause = "entry.taxon_id"
        else:
            raise RuntimeError(f"Searching on {search_on} not supported, must be taxid or ko_id")

        if keep_seqids:
            keep_sel = " hsh.seqid, "
            keep_grp = " , hsh.seqid"
            ids = ["taxid", "KO", "seqid", "count"]
        else:
            keep_sel = ""
            keep_grp = ""
            ids = ["taxid", "KO", "count"]

        entries = self.gen_placeholder_string(search_entries)
        query = (
            f"SELECT entry.taxon_id, hit.ko_id, {keep_sel} count(*) "
            "FROM seqfeature AS feature "
            "INNER JOIN sequence_hash_dictionnary AS hsh ON feature.seqfeature_id = hsh.seqid "
            "INNER JOIN ko_hits AS hit ON hit.hsh = hsh.hsh "
            "INNER JOIN bioentry AS entry ON entry.bioentry_id=feature.bioentry_id "
            f"WHERE {where_clause} IN ({entries})"
            f"GROUP BY entry.taxon_id, hit.ko_id {keep_grp};"
        )
        results = self.server.adaptor.execute_and_fetchall(query, search_entries)
        if not as_multi:
            return DB.to_pandas_frame(results, ids)
        return DB.to_pandas_frame(results, ids).set_index(["taxid", "KO"])


    def get_seqids_for_ko(self, ko_ids, only_seqids=False):
        if only_seqids:
            selection = ""
        else:
            selection = ", hits.ko_id "

        entries = ",".join("?" for i in ko_ids)
        query = (
            f"SELECT hsh.seqid {selection}"
            "FROM ko_hits AS hits "
            "INNER JOIN sequence_hash_dictionnary as hsh ON hsh.hsh = hits.hsh "
            f"WHERE ko_id IN ({entries}) GROUP BY hsh.seqid;"
        )
        results = self.server.adaptor.execute_and_fetchall(query, ko_ids)

        if only_seqids:
            return [line[0] for line in results]
        else:
            hsh_results = {}
            for line in results:
                assert line[0] not in hsh_results
                hsh_results[line[0]] = line[1]
            return hsh_results


    def get_hsh_locus_to_seqfeature_id(self, only_CDS=False):
        filtering = ""
        if only_CDS:
            filtering = "INNER JOIN term as t4 ON t4.term_id=t1.type_term_id AND t4.name = \"CDS\" "

        query = (
            "SELECT t2.value, t1.seqfeature_id "
            "FROM seqfeature as t1 "
            f"{filtering}"
            "INNER JOIN seqfeature_qualifier_value AS t2 ON t1.seqfeature_id=t2.seqfeature_id "
            "INNER JOIN term as t3 on t3.term_id=t2.term_id AND t3.name=\"locus_tag\" "
        )
        results = self.server.adaptor.execute_and_fetchall(query,)
        hsh_results = {}
        for line in results:
            hsh_results[line[0]] = int(line[1])
        return hsh_results


    def get_og_phylogeny(self, og):
        query = (
            "SELECT tree FROM gene_phylogeny WHERE orthogroup_id=?;"
        )
        results = self.server.adaptor.execute_and_fetchall(query, [og])
        if len(results)!=1:
            raise RuntimeError(f"Could not find phylogeny for orthogroup {og}")
        return results[0][0]


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


    def get_config_table(self, ret_mandatory=False):
        sql = "SELECT * from biodb_config;"
        values = self.server.adaptor.execute_and_fetchall(sql)
        if ret_mandatory:
            return {val[0]: (val[1]=="mandatory", val[2]) for val in values}
        else:
            return {val[0]: val[2] for val in values}


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

    
    def get_locus_to_genomes(self, locus_lst):
        # quick and dirty, may need to regroup it with another function
        values = ",".join("?" for _ in locus_lst)
        query = (
            "SELECT locus_tag.value, entry.description "
            "FROM seqfeature_qualifier_value AS locus_tag "
            "INNER JOIN term AS locus_term ON locus_term.term_id=locus_tag.term_id "
            " AND locus_term.name=\"locus_tag\" "
            "INNER JOIN seqfeature AS feature ON feature.seqfeature_id=locus_tag.seqfeature_id "
            "INNER JOIN bioentry AS entry ON feature.bioentry_id=entry.bioentry_id "
            f"WHERE locus_tag.value IN ({values});"
        )
        results = self.server.adaptor.execute_and_fetchall(query, locus_lst)
        return {locus: description for locus, description in results}

    
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


    def get_genomes_description(self, lst_plasmids=False):
        """
        Returns the description of the genome as it has been read from the genbank
        files, indexed by taxon_id. The output also contains a flag has_plasmid
        indicating whether the genome contains a plasmid or not, if the lst_plasmid flag
        has been set.
        """

        has_plasmid_query = (
            "SELECT * "
            "FROM bioentry_qualifier_value AS has_plasmid "
            "INNER JOIN term AS pls_term ON pls_term.term_id=has_plasmid.term_id "
            " AND pls_term.name=\"plasmid\" " 
            "INNER JOIN bioentry AS plasmid ON has_plasmid.bioentry_id=plasmid.bioentry_id "
            "WHERE plasmid.taxon_id=entry.taxon_id"
        )
        query = (
            "SELECT entry.taxon_id, entry.description,  "
            f" CASE WHEN EXISTS ({has_plasmid_query}) THEN 1 ELSE 0 END "
            "FROM bioentry AS entry "
            "INNER JOIN bioentry_qualifier_value AS orga " 
            "INNER JOIN term AS orga_term ON orga.term_id=orga_term.term_id "
            " AND orga_term.name=\"organism\" "
            "GROUP BY taxon_id;" 
        )
        descr = self.server.adaptor.execute_and_fetchall(query)
        columns = ["taxon_id", "description"]
        if lst_plasmids:
            columns.append("has_plasmid")
        else:
            # remove the third row to keep pandas happy
            descr = ((taxon_id, entry_desc) for taxon_id, entry_desc, _ in descr)

        return DB.to_pandas_frame(descr, columns).set_index(["taxon_id"])


    def get_genomes_infos(self):
        """
        Note: for efficiency sake, it would be possible to order the different
         tables by taxon_id, walk through the table using zip in an iteator
         and pass the iterator to pandas.

         However, given the small size of the dataset, it's probably not worth
         the implementation effort.
        """

        query = (
            "SELECT entry.taxon_id, COUNT(*) "
            " FROM seqfeature AS seq "
            " INNER JOIN term AS cds ON cds.term_id = seq.type_term_id AND cds.name=\"CDS\" "
            " INNER JOIN bioentry AS entry ON entry.bioentry_id = seq.bioentry_id " 
            " GROUP BY entry.taxon_id;"
        )
        n_prot_results = self.server.adaptor.execute_and_fetchall(query)

        query = (
            "SELECT entry.taxon_id, COUNT(*) "
            "FROM bioentry AS entry "
            "GROUP BY entry.taxon_id;"
        )
        n_contigs = self.server.adaptor.execute_and_fetchall(query)

        cols = ["taxon_id", "completeness", "contamination", "gc", "length", "coding_density"]
        query = (
            "SELECT taxon_id, completeness, contamination, gc, length, coding_density "
            " from genome_summary;"
        )
        all_other_results = self.server.adaptor.execute_and_fetchall(query)

        df_n_prot = DB.to_pandas_frame(n_prot_results, ["taxon_id", "n_prot"]).set_index("taxon_id")
        df_n_contigs = DB.to_pandas_frame(n_contigs, ["taxon_id", "n_contigs"]).set_index("taxon_id")
        df_stats = DB.to_pandas_frame(all_other_results, cols).set_index("taxon_id")
        return df_n_prot.join(df_n_contigs).join(df_stats)


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
        query = "SELECT COUNT(*) FROM bioentry GROUP BY taxon_id;"
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
            "CREATE TABLE filenames (taxon_id INTEGER, filename TEXT);"
        )
        self.server.adaptor.execute(sql)
        self.load_data_into_table("filenames", data)


    def load_genomes_info(self, data):
        sql = (
            "CREATE TABLE genome_summary (taxon_id INTEGER, completeness FLOAT, "
            " contamination FLOAT, gc INTEGER, length INTEGER, coding_density FLOAT);"
        )
        self.server.adaptor.execute(sql)
        self.load_data_into_table("genome_summary", data)


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


    def get_cog_summaries(self, cog_ids, only_cog_desc=False, as_df=False):
        ids = ",".join(["?"] * len(cog_ids))
        query = (
            "SELECT cog_id, function, description "
            "FROM cog_names "
            f"WHERE cog_id IN ({ids});"
        )
        results = self.server.adaptor.execute_and_fetchall(query, cog_ids)
        if only_cog_desc and not as_df:
            hsh_results = {}
            for line in results:
                hsh_results[line[0]] = (line[1], line[2])
            return hsh_results
        elif only_cog_desc:
            return DB.to_pandas_frame(results, ["cog", "function", "description"])

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


    def get_CDS_from_locus_tag(self, locus_tag):
        # NOTE: I did not add a join to filter on locus_tag,
        # it may be worth it performance-wise to pre-filter the 
        # database entries if this query were to become an issue.

        query = (
            "SELECT seqfeature_id "
            "FROM seqfeature_qualifier_values AS locus_tag "
            "INNER JOIN seqfeature AS feature ON feature.seqfeature_id=locus_tag.seqfeature_id "
            "INNER JOIN term AS t ON feature.type_term_id=t.term_id AND t.name = \"CDS\" "
            "WHERE locus_tag.value=?;"
        )
        ret = self.server.execute_and_fetchall(query, locus_tag)
        if ret==None or len(ret)==0:
            return None
        return ret[0][0]


    def get_seqid(self, locus_tag):
        query = (
            "SELECT locus_tag.seqfeature_id "
            "FROM seqfeature_qualifier_value AS locus_tag "
            "INNER JOIN seqfeature AS cds "
            " ON cds.seqfeature_id=locus_tag.seqfeature_id "
            "INNER JOIN term AS cds_term ON cds_term.term_id=cds.type_term_id "
            " AND cds_term.name=\"CDS\""
            "WHERE locus_tag.value = ?;"
        )

        values =  self.server.adaptor.execute_and_fetchall(query, [locus_tag,])
        if len(values)==0:
            raise RuntimeError("No such entry")
        return values[0][0]

    
    def get_bioentry(self, from_val, val_type="seqid"):
        if val_type != "seqid":
            raise Exception("For now, only seqid indexing is supported")
        # may be extended in the future for other types of indexing

        sql = (
            "SELECT bioentry_id "
            "FROM seqfeature "
            "WHERE seqfeature_id = ?;"
        )
        values = self.server.adaptor.execute_and_fetchall(sql, [from_val])

        # assumes result has been found
        return values[0][0]


    def get_seqid_in_neighborhood(self, bioentry_id, start_loc, stop_loc):
        query = (
            "SELECT seqfeature_id "
            "FROM location AS loc " 
            "INNER JOIN seqfeature AS seq ON seq.seqfeature_id=loc.seqfeature_id "
            "  AND seq.bioentry_id= ? "
            "WHERE loc.start_pos > ? AND loc.start_pos < ?;"
        )
        results = self.server.execute_and_fetchall(query, [bioentry_id, start_pos, stop_loc])
        if results==None:
            return []
        return [line[0] for line in results]


    def get_DNA_sequence(self, bioentry_id, alphabet="dna"):
        query = (
            "SELECT seq FROM biosequence WHERE bioentry_id=? AND alphabet = ?;"
        )

        results = self.server.adaptor.execute_and_fetchall(query, [bioentry_id, alphabet])
        return Seq(results[0][0])


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
    def get_proteins_info(self, seqids, bioentries=None):
        seqids_query = ",".join(["?"] * len(seqids))
        term_names = ["locus_tag", "protein_id", "gene", "product"]
        term_names_query = ",".join([f"\"{name}\"" for name in term_names])

        if bioentries!=None:
            entries = ",".join(str(entry) for entry in bioentries)
            sel = (
                f"INNER JOIN seqfeature AS seq ON seq.seqfeature_id = v.seqfeature_id "
                f" AND seq.bioentry_id IN ({entries}) "
            )
        else:
            sel = ""
        query = (
            "SELECT v.seqfeature_id, t.name, v.value "
            "FROM seqfeature_qualifier_value AS v "
            "INNER JOIN term AS t ON t.term_id = v.term_id "
            f"{sel}"
            f"WHERE v.seqfeature_id IN ({seqids_query}) AND name IN ({term_names_query}) "
            "ORDER BY v.seqfeature_id ASC;"
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
    # elements separately.
    #
    # NOTE: may be interesting to use int8/16 whenever possible 
    # to spare memory.
    def to_pandas_frame(db_results, columns, types=None):
        return pd.DataFrame(db_results, columns=columns)


    def get_bioentries_in_taxon(self, bioentries=None):
        # NOTE: need to write the code for the base where bioentries is None
        # -> returns all the entries
        query_str = ",".join("?" for entry in bioentries)
        query = (
            "SELECT two.bioentry_id, one.taxon_id, one.bioentry_id "
            "FROM bioentry AS one "
            "INNER JOIN bioentry AS two ON one.taxon_id = two.taxon_id "
            f"WHERE one.bioentry_id IN ({query_str});"
        )
        results = self.server.adaptor.execute_and_fetchall(query, bioentries)
        return DB.to_pandas_frame(results, ["bioentry", "taxon", "ref_genome_bioentry"])


    def get_og_identity(self, og, ref_seqid):
        """
        For now need to have both an og and a ref_seqid.
        """

        query = (
            "SELECT id_1, id_2, identity "
            "FROM orthology_identity "
            "WHERE orthogroup = ? AND (id_1 = ? OR id_2 = ?);"
        )
        results = self.server.adaptor.execute_and_fetchall(query, (og, ref_seqid, ref_seqid))
        filtered_values = []
        for id_1, id_2, identity in results:
            if id_1==ref_seqid:
                filtered_values.append((id_2, identity))
            else:
                filtered_values.append((id_1, identity))
        df = DB.to_pandas_frame(filtered_values, ["seqid", "identity"])
        return df.set_index(["seqid"])


    def get_og_count(self, lookup_terms, search_on="taxid", plasmids=None, keep_taxid=False):
        """
        This function returns a pandas dataframe containing the orthogroup
        count for a set of taxon_ids. The user can differentiate between the chromosome
        and the plasmid.

        plasmids: if set to None, the function will not differentiate between
            plasmids and chromosomes. Otherwise, should contain a list of taxids whose plasmids
            will be taken into account in the search.
        lookup_term: the terms that will be search. Either a list of taxids, of orthogroup or of seqids.
        search_on: can be either orthogroup, seqid or taxid. Specifies what lookup_term is.
        keep_taxid: for queries using the seqid lookup term, also return the taxid in the results

        The index of the table depends on whether the diff_plasmid flag was set. 
        If the flag is set, MultiIndex(taxid, is_plasmid). If it is not set Index(taxid)
        """

        if not plasmids is None and search_on!="taxid":
            raise RuntimeError("Plasmid search only supported for taxid")

        entries = self.gen_placeholder_string(lookup_terms)

        add_plasmid = ""
        if not plasmids is None:
            add_plasmid = "CAST(is_plasmid.value AS int), "

        select   = f"SELECT entry.taxon_id, {add_plasmid} orthogroup, COUNT(*) "
        group_by = f"GROUP BY entry.taxon_id, {add_plasmid} orthogroup"
        where_clause = None
        add_plasmid_join = ""
        if not plasmids is None:
            plasmids_placeholder = self.gen_placeholder_string(plasmids)
            add_plasmid = "is_plasmid.value"
            add_plasmid_join = (
                "INNER JOIN bioentry_qualifier_value AS is_plasmid ON "
                "  is_plasmid.bioentry_id=entry.bioentry_id "
                "INNER JOIN term AS plasmid_term ON plasmid_term.term_id=is_plasmid.term_id "
                "  AND plasmid_term.name=\"plasmid\""
            )
            where_clause = (
                f" (entry.taxon_id IN ({entries}) AND is_plasmid.value=0) "
                f"OR (entry.taxon_id IN ({plasmids_placeholder}) AND is_plasmid.value=1)"
            )
        elif search_on=="taxid":
            where_clause = f"entry.taxon_id IN ({entries})"
        elif search_on=="orthogroup":
            where_clause = f"og.orthogroup IN ({entries})"
        elif search_on=="seqid":
            where_clause = f"feature.seqfeature_id IN ({entries})"
            select       = "SELECT feature.seqfeature_id, og.orthogroup "
            group_by     = ""
            if keep_taxid:
                select += ", entry.taxon_id "
        else:
            raise RuntimeError(f"Unsupported search {search_on}, must use orthogroup, bioentry or seqid")

        query = (
            f"{select} "
            "FROM bioentry AS entry "
            "INNER JOIN seqfeature AS feature ON entry.bioentry_id = feature.bioentry_id "
            "INNER JOIN og_hits AS og ON og.seqid = feature.seqfeature_id "
            f"{add_plasmid_join} "
            f"WHERE {where_clause} "
            f"{group_by};"
        )
        all_terms = lookup_terms
        if not plasmids is None:
            all_terms = lookup_terms + plasmids
        results = self.server.adaptor.execute_and_fetchall(query, all_terms)

        header = None
        if not plasmids is None:
            header = ["taxid", "plasmid", "orthogroup", "count"]
        elif search_on=="taxid" or search_on=="orthogroup":
            header = ["taxid", "orthogroup", "count"]
        elif search_on=="seqid":
            header = ["seqid", "orthogroup"]
            if keep_taxid:
                header.append("taxid")

        df = DB.to_pandas_frame(results, header)
        if len(df.index) == 0:
            return df
        if not plasmids is None:
            df = df.set_index(["taxid", "plasmid", "orthogroup"]).unstack(level=0, fill_value=0)
            return df.unstack(level=0, fill_value=0)
        elif search_on=="taxid" or search_on=="orthogroup":
            df = df.set_index(["taxid", "orthogroup"]).unstack(level=0, fill_value=0)
            df.columns = [col for col in df["count"].columns.values]
        elif search_on=="seqid":
            df = df.set_index(["seqid"])
        return df


    def get_translation(self, seqid):
        query = (
            "SELECT translation.value  "
            "FROM seqfeature_qualifier_value AS translation "
            "INNER JOIN term AS transl_term ON transl_term.term_id=translation.term_id "
            " AND transl_term.name=\"translation\" "
            f"WHERE translation.seqfeature_id = ?;"
        )
        results = self.server.adaptor.execute_and_fetchall(query, [seqid])
        if results==None or len(results)==0:
            raise RuntimeError("No translation")
        return results[0][0]


    def get_genes_from_og(self, orthogroups, taxon_ids=None, terms=["gene", "product"]):
        for term in terms:
            if term not in ["length", "gene", "product", "locus_tag"]:
                raise RuntimeError(f"Term not supported: {term}")

        og_entries = self.gen_placeholder_string(orthogroups)
        query_args = orthogroups

        db_terms = [t for t in terms if t!="length"]
        if "length" in terms:
            db_terms.append("translation")
        sel_terms = ",".join("\"" + f"{i}" + "\"" for i in db_terms)
        join, sel = "", ""
        if taxon_ids != None:
            taxon_id_query = self.gen_placeholder_string(taxon_ids)
            join = (
                "INNER JOIN seqfeature AS seq ON og.seqid=seq.seqfeature_id "
                "INNER JOIN bioentry AS entry ON seq.bioentry_id=entry.bioentry_id "
            )
            sel = f"AND entry.taxon_id IN ({taxon_id_query})"
            query_args = orthogroups+taxon_ids
        
        query = (
            f"SELECT feature.seqfeature_id, og.orthogroup, t.name, "
            "   CASE WHEN t.name==\"translation\" THEN "
            "   LENGTH(value) ELSE value END "
            "FROM og_hits AS og "
            "INNER JOIN seqfeature_qualifier_value AS feature ON og.seqid = feature.seqfeature_id "
            f"INNER JOIN term AS t ON t.term_id=feature.term_id AND t.name IN ({sel_terms})"
            f"{join}"
            f"WHERE og.orthogroup IN ({og_entries}) {sel};"
        )
        results = self.server.adaptor.execute_and_fetchall(query, query_args)
        df = DB.to_pandas_frame(results, columns=["seqid", "orthogroup", "term", "value"])
        df = df.set_index(["orthogroup", "seqid", "term"]).unstack("term", fill_value=None)
        if "value" in df.columns:
            df.columns = ["length" if col=="translation" else col for col in df["value"].columns.values]

        df = df.reset_index(level=0)
        for t in terms:
            if t not in df.columns:
                df[t] = None
        return df


    def get_cog_counts_per_category(self, taxon_ids):
        entries = self.gen_placeholder_string(taxon_ids)
        query = (
            "SELECT entry.taxon_id, cog.function "
            "FROM seqfeature AS feature "
            "INNER JOIN bioentry AS entry ON feature.bioentry_id=entry.bioentry_id "
            "INNER JOIN sequence_hash_dictionnary AS hsh on hsh.seqid = feature.seqfeature_id "
            "INNER JOIN cog_hits AS hit ON hit.hsh = hsh.hsh "
            "INNER JOIN cog_names AS cog ON cog.cog_id = hit.cog_id "
            f"WHERE entry.taxon_id IN ({entries});"
        )
        results = self.server.adaptor.execute_and_fetchall(query, taxon_ids)
        hsh_results = {}
        for line in results:
            entry_id, func = line
            if entry_id not in hsh_results:
                hsh_results[entry_id] = {func : 1}
            else:
                cnt = hsh_results[entry_id].get(func, 0)
                hsh_results[entry_id][func] = cnt+1
        return hsh_results


    # Get all cog hits for a given list of bioentries
    # The results are either indexed by the bioentry or by the seqid
    # NOTE: if indexing as bioentry or taxon_id, will return a dataframe of the form
    #           bioentry_1/taxon_1  bioentry_2/taxon_2  ...
    #   cog_1     cnt                   cnt
    #   cog_2     cnt                   cnt
    #   ...
    #
    # If the indexing is seqid, will return a dataframe with the format
    #   seqid1 cog1
    #   seqid2 cog2
    #   seqid3 cog3
    def get_cog_hits(self, ids, indexing="bioentry", search_on="bioentry", keep_taxid=False):
        entries = self.gen_placeholder_string(ids)

        if search_on=="bioentry":
            where_clause = f" entry.bioentry_id IN ({entries}) "
        elif search_on=="seqid":
            where_clause = f" hsh.seqid IN ({entries}) "
        elif search_on=="cog":
            where_clause = f" cogs.cog_id IN ({entries}) "
        elif search_on=="taxid":
            where_clause = f" entry.taxon_id IN ({entries}) "
        else:
            raise RuntimeError(f"Searching on {search_on} is not supported")

        if indexing=="seqid":
            index = "seqid.seqfeature_id"
            if keep_taxid:
                index += ", entry.taxon_id "
        elif indexing=="bioentry":
            index = "entry.bioentry_id"
        elif indexing=="taxid":
            index = "entry.taxon_id"
        else:
            raise RuntimeError(f"Indexing method not supported: {indexing}")

        query = (
            f"SELECT {index}, cogs.cog_id, COUNT(*) "
            "FROM bioentry AS entry "
            "INNER JOIN seqfeature AS seqid ON seqid.bioentry_id = entry.bioentry_id "
            "INNER JOIN sequence_hash_dictionnary AS hsh ON seqid.seqfeature_id = hsh.seqid "
            "INNER JOIN cog_hits AS cogs ON cogs.hsh = hsh.hsh "
            f"WHERE {where_clause} "
            f"GROUP BY {index}, cogs.cog_id;"
        )
        results = self.server.adaptor.execute_and_fetchall(query, ids)
        if indexing=="taxid" or indexing=="bioentry":
            column_names = [indexing, "cog", "count"]
            df = DB.to_pandas_frame(results, column_names)
            df = df.set_index([indexing, "cog"]).unstack(level=0, fill_value=0)
            df.columns = [col for col in df["count"].columns.values]
        elif indexing=="seqid":
            header = ["seqid", "cog"]
            if keep_taxid:
                header.append("taxid")
                results = ((seqid, cog, taxid) for seqid, taxid, cog, count in results)
            else:
                results = ((seqid, cog) for seqid, cog, count in results)

            df = DB.to_pandas_frame(results, header)
            df = df.set_index(["seqid"])
        return df


    def get_filenames_to_taxon_id(self):
        sql = (
            f"SELECT filename, taxon_id FROM filenames;"
        )
        hsh_filenames_to_entry = {}
        results = self.server.adaptor.execute_and_fetchall(sql)
        return {filename: entry_id for filename, entry_id in results}


    def load_chlamdb_config_tables(self, entries):
        sql = (
            "CREATE TABLE biodb_config"
            "(name varchar(200), type varchar(200), status BOOLEAN);"
        )
        self.server.adaptor.execute(sql)
        self.load_data_into_table("biodb_config", entries)


    def update_taxon_ids(self, update_lst):
        query = (
            "UPDATE bioentry SET taxon_id=? WHERE bioentry_id=?;"
        )
        for bioentry_id, taxon_id in update_lst:
            self.server.adaptor.execute(query, (taxon_id, bioentry_id))


    def update_plasmid_status(self, plasmid_bioentries):
        plasmid_term_id = self.get_term_id("plasmid", create_if_absent=True)
        sql = (
            "INSERT INTO bioentry_qualifier_value VALUES  "
            " (?, ?, ?, 0);"
        )
        data = [(bioentry_id, plasmid_term_id, is_plasmid)
                for bioentry_id, is_plasmid in plasmid_bioentries]
        self.server.adaptor.executemany(sql, data)

    
    def gen_placeholder_string(self, args):
        return ",".join(self.placeholder for _ in args)


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
        return DB.load_db(db_name, params)

