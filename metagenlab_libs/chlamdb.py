#!/usr/bin/env python

# This script is used to generate the skeleton of the database
# used by chlamdb. The entries will be filled by the annotation
# pipeline.
# 
# For now, this script creates a new mysql database and creates
# the tables of the biosql schema.
# It then pulls all kegg orthologs, their related modules and pathways
# and enters them in the database.
# Same with COGs.
#
# Bastian Marquis (bastian.marquis@protonmail.com)
# Date: 29.10.2020

import os
import db_utils
import argparse

import queue
import threading

from Bio.KEGG import REST
from Bio.KEGG import Gene

DEFAULT_PSSWD = ""
DEFAULT_TYPE = "sqlite"
DEFAULT_DB_NAME = "George"


def create_data_table(db):
    entry_list = [
        ("orthology", "mandatory", False),
        ("orthogroup_alignments", "mandatory", False),
        ("old_locus_table", "mandatory", False),
        ("reference_phylogeny", "mandatory", False),
        ("taxonomy_table", "mandatory", False),
        ("genome_statistics", "mandatory", False),
        ("BLAST_database", "optional", False),
        ("gene_phylogenies", "optional", False),
        ("interpro_data", "optional", False),
        ("interpro_comparative", "optional", False),
        ("interpro_comparative_accession", "optional", False),
        ("priam_data", "optional", False),
        ("priam_comparative", "optional", False),
        ("priam_comparative_accession", "optional", False),
        ("COG", "optional", False),
        ("KEGG", "optional", False),
        ("pfam_comparative", "optional", False),
        ("pfam_comparative_accession", "optional", False),       
        ("TCDB_data", "optional", False),
        ("psortb_data", "optional", False),
        ("T3SS_data", "optional", False),
        ("PDB_data", "optional", False),
        ("BLAST_refseq", "optional", False),
        ("BLAST_swissprot", "optional", False),
        ("BBH_phylogenies", "optional", False),
        ("GC_statistics", "optional", False),
        ("gene_clusters", "optional", False),
        ("phylogenetic_profile", "optional", False),
        ("synonymous_table", "optional", False),
        ("interpro_taxonomy", "optional", False), # interpro taxnonomy statistics
        ("pfam_taxonomy", "optional", False), #  taxnonomy statistics
        ("COG_taxonomy", "optional", False) # COG taxnonomy statistics
    ]
    db.load_chlamdb_config_tables(entry_list)
    db.commit()

def setup_biodb(kwargs):
    sqlpsw = kwargs.get("db_psswd", DEFAULT_PSSWD)
    db_type = kwargs.get("db_type", DEFAULT_TYPE)
    db_name = kwargs.get("db_name", DEFAULT_DB_NAME)
    schema_dir = kwargs.get("biosql_schema_dir", "biosql_schema")

    if db_type=="sqlite":
        import sqlite3
        conn = sqlite3.connect(db_name)
        cursor = conn.cursor()
        url_biosql_scheme = 'biosqldb-sqlite.sql'
        err_code = os.system(f"sqlite3 {db_name} < {schema_dir}/{url_biosql_scheme}")
        conn.execute("pragma journal_mode=wal")
    else:
        # NOTE: this part is untested!
        import MySQLdb
        conn = MySQLdb.connect(host="localhost", # your host, usually localhost
                                    user="root", # your username
                                    passwd=sqlpsw) # name of the data base
        cursor = conn.cursor()
        sql_db = f'CREATE DATABASE IF NOT EXISTS {db_name};'
        cursor.execute(sql_db,)
        conn.commit()
        cursor.execute(f"use {db_name};",)
        url_biosql_scheme = 'biosqldb-mysql.sql'
        err_code = os.system(f"mysql -uroot -p{sqlpsw} {db_name} < {schema_dir}/{url_biosql_scheme}")

    if err_code != 0:
        raise IOError("Problem loading sql schema:", err_code)

    # not really logical to me, but creating a database
    # from the biosql is necessary
    chlamdb_args = {"chlamdb.db_type": db_type,
            "chlamdb.db_name": db_name,
            "chlamdb.db_psswd": sqlpsw}
    db = db_utils.DB.load_db(chlamdb_args)
    db.create_biosql_database(chlamdb_args)
    db.commit()
    return db


# code imported from the Bio.KEGG module
class Record(object):
    def __init__(self):
        self.entry = ""
        self.name = []
        self.definition = ""
        self.modules = []
        self.pathways = []

    def simplified_entry(self):
        return int(self.entry[len("K"):])

    def __str__(self):
        acc = self.entry + " " + self.definition
        return acc


class Module(object):
    def __init__(self, entry, descr):
        self.entry = entry
        self.descr = descr
        self.classes = []
        self.definition = ""

    def simplified_entry(self):
        return int(self.entry[len("M"):])

    def __str__(self):
        return self.entry + ": " + self.descr

    def __hash__(self):
        return hash(self.entry)

    def __eq__(self, other):
        return self.entry == other.entry

    def parse(self, handle):
        for line in handle:
            if line[:3] == "///":
                return
            if line[:12] != "            ":
                keyword = line[:12]
            data = line[12:].strip()
            if keyword == "CLASS       ":
                tokens = data.split(";")
                self.classes = [token.strip() for token in tokens]
            if keyword == "DEFINITION  ":
                self.definition = data


class Pathway(object):
    def __init__(self, entry, descr):
        self.entry = entry
        self.descr  = descr

    def simplified_entry(self):
        return int(self.entry[len("ko"):])

    def __hash__(self):
        return hash(self.entry)

    def __eq__(self, other):
        return self.entry == other.entry

    def __str__(self):
        return self.entry + ": " + self.descr


def parse_gene(handle):
    record = Record()
    for line in handle:
        if line[:3] == "///":
            yield record
            record = Record()
            continue
        if line[:12] != "            ":
            keyword = line[:12]
        data = line[12:].strip()
        if keyword == "ENTRY       ":
            words = data.split()
            record.entry = words[0]
        elif keyword == "NAME        ":
            data = data.strip(";")
            record.name.append(data)
        elif keyword == "DEFINITION  ":
            record.definition = data
        elif keyword == "MODULE      ":
            module, *descr_tokens = data.split()
            descr = " ".join(descr_tokens)
            record.modules.append(Module(module, descr))
        elif keyword == "PATHWAY     ":
            pathway, *descr_tokens = data.split()
            descr = (" ".join(descr_tokens)).strip()
            record.pathways.append(Pathway(pathway, descr))


def load_KO_references(db, params):
    ko_string = REST.kegg_list("KO")
    ko_codes = []
    for line in ko_string:
        ko_codes.append(line.split()[0])

    genes = []
    for i, ko_code in enumerate(ko_codes):
        buf = REST.kegg_get(ko_code)
        for gene in parse_gene(buf):
            genes.append(gene)
        if i==15:
            # test version: do not download all keggs
            break

    modules_set = {module for gene in genes for module in gene.modules}
    pathway_set = {pathway for gene in genes for pathway in gene.pathways}
    for module in modules_set:
        mod_data = REST.kegg_get(module.entry)
        module.parse(mod_data)

    db.load_ko_module([(m.simplified_entry(), m.descr, m.definition) for m in modules_set])
    db.load_ko_pathway([(p.simplified_entry(), p.descr) for p in pathway_set])
    db.load_ko_def([(gene.simplified_entry(), gene.definition) for gene in genes])

    ko_to_path = []
    ko_to_module = []
    for gene in genes:
        simp = gene.simplified_entry()
        ko_to_path.extend((simp, path.simplified_entry()) for path in gene.pathways)
        ko_to_module.extend((simp, mod.simplified_entry()) for mod in gene.modules)
    db.load_ko_to_pathway(ko_to_path)
    db.load_ko_to_module(ko_to_module)
    db.commit()
    return genes


def setup_cog(db, cog_dir):
    cog2cdd_file = open(cog_dir+"/cog_corresp.tab", "r")
    cog2length_file = open(cog_dir+"/cog_length.tab", "r")
    fun_names_file = open(cog_dir+"/fun2003-2014.tab")
    cog_names_file = open(cog_dir+"/cognames2003-2014.tab")

    cdd_to_cog = []
    for line in cog2cdd_file:
        tokens = line.split("\t")
        cdd_to_cog.append( (int(tokens[1].strip()), int(tokens[0][len("COG"):])) )
    db.load_cdd_to_cog(cdd_to_cog)

    hsh_cog_to_length = {}
    for line in cog2length_file:
        tokens = line.split("\t")
        hsh_cog_to_length[int(tokens[0][3:])] = int(tokens[1])

    cog_ref_data = []
    # necessary to track, as some cogs listed in the CDD to COG mapping table
    # are not present in those descriptors
    hsh_cog_ids = {}
    for line_no, line in enumerate(cog_names_file):
        # pass header
        if line_no == 0:
            continue
        tokens = line.split("\t")
        cog_id = int(tokens[0][3:])
        hsh_cog_ids[cog_id] = True
        fun = tokens[1].strip()
        description = tokens[2].strip()
        cog_ref_data.append( (cog_id, fun, description) )
    db.load_cog_ref_data(cog_ref_data)

    cog_fun_data = []
    for line_no, line in enumerate(fun_names_file):
        if line_no==0:
            continue
        function, description = line.split("\t") 
        cog_fun_data.append( (function.strip(), description.strip()) )
    db.load_cog_fun_data(cog_fun_data)
    db.commit()


parser = argparse.ArgumentParser(description = "Creates a chlamdb database skeleton")

parser.add_argument("--db_name", nargs="?", default=DEFAULT_DB_NAME,
        help=f"name of the database (default name {DEFAULT_DB_NAME})")

parser.add_argument("--load_cog", action="store_true",
        help="load cog definitions (default no)")

parser.add_argument("--cog_dir", nargs="?", default="./",
        help="directory where the cog definitions files are (default current directory)")

parser.add_argument("--load_kegg", action="store_true",
        help="load kegg definitions (default yes)")

parser.add_argument("--db_type", nargs="?", default=DEFAULT_TYPE,
        help="database type (either sqlite or mysql)")

parser.add_argument("--db_psswd", nargs="+", default=DEFAULT_PSSWD,
        help="set db password (default none)")

args = vars(parser.parse_args())

print("Setting up the biosql schema")
db = setup_biodb(args)
create_data_table(db)

if args.get("load_cog", False):
    print("Loading COG tables")
    setup_cog(db, args.get("cog_dir", "."))

if args.get("load_kegg", False):
    print("Loading KEGG tables")
    genes = load_KO_references(db, args)
