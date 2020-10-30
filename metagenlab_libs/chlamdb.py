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

from Bio.KEGG import REST
from Bio.KEGG import Gene

DEFAULT_PSSWD = ""
DEFAULT_TYPE = "sqlite"
DEFAULT_DB_NAME = "George"


def create_data_table(db, kwargs):
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
    db = db_utils.DB.load_db(kwargs)
    db.create_biosql_database(kwargs)
    db.commit()
    return db


def simplify_pathway(raw_pathway):
    return int(raw_pathway[len("path:map"):])


def simplify_module(raw_module):
    return int(raw_module[len("md:M"):])


def parse_REST_result(string, c1_name, c2_name):
    entries = string.strip().split("\n")
    to_cols = pd.DataFrame([entry.split("\t") for entry in entries],
            columns=[c1_name, c2_name])
    return to_cols


def simplify_desc(string):
    tokens = string.split(";")
    return tokens[-1]


def parse_pathway(buff):
    pass


def parse_module(buff):
    pass


def load_KO_references(params):
    db = db_utils.DB.load_db(params)
    ko_string = REST.kegg_list("KO").read()
    ko_id_to_desc = parse_REST_result(ko_string, "ko", "desc")
    ko_id_to_desc["ko_simplified"] = ko_id_to_desc["ko"].apply(simplify_ko)
    ko_id_to_desc["desc"] = ko_id_to_desc["desc"].apply(simplify_desc)
    db.load_ko_def(ko_id_to_desc[["ko_simplified", "desc"]].values.tolist())

    path_string = REST.kegg_list("pathway").read()
    path_id_to_desc = parse_REST_result(path_string, "pathway", "desc")
    path_id_to_desc["pathway"] = path_id_to_desc["pathway"].apply(simplify_pathway)
    db.load_ko_pathway(path_id_to_desc.values.tolist())

    mod_string = REST.kegg_list("module").read()
    mod_id_to_desc = parse_REST_result(mod_string, "module", "desc")
    mod_id_to_desc["module"] = mod_id_to_desc["module"].apply(simplify_module)
    db.load_ko_module(mod_id_to_desc.values.tolist())

    ko_to_module_data = []
    ko_to_pathway_data = []
    for lst_ko in chunks(ko_id_to_desc["ko"].unique().tolist(), 200):
        ko_to_module = REST.kegg_link("module", lst_ko).read()
        if len(ko_to_module) > 1:
            ko_to_module_desc = parse_REST_result(ko_to_module, "ko", "module")
            ko_to_module_desc["ko"] = ko_to_module_desc["ko"].apply(simplify_ko)
            ko_to_module_desc["module"] = ko_to_module_desc["module"].apply(simplify_module)
            ko_to_module_data.append(ko_to_module_desc.drop_duplicates().values.tolist())

        ko_to_pathway = REST.kegg_link("pathway", lst_ko).read()
        if len(ko_to_pathway) > 1:
            ko_to_pathway_desc = parse_REST_result(ko_to_pathway, "ko", "pathway")
            ko_to_pathway_desc["ko"] = ko_to_pathway_desc["ko"].apply(simplify_ko)
            ko_to_pathway_desc = ko_to_pathway_desc[ko_to_pathway_desc.pathway.str.startswith("path:map")]
            ko_to_pathway_desc["pathway"] = ko_to_pathway_desc["pathway"].apply(simplify_pathway)
            ko_to_pathway_data.append(ko_to_pathway_desc.drop_duplicates().values.tolist())
    db.load_ko_to_pathway([elem for lst in ko_to_pathway_data for elem in lst])
    db.load_ko_to_module([elem for lst in ko_to_module_data for elem in lst])
    db.commit()


def setup_cog(params):
    db = db_utils.DB.load_db(params)

    cog2cdd_file = open(params["cog_dir"]+"/COG/cog_corresp.tab", "r")
    cog2length_file = open(params["cog_dir"]+"/COG/cog_length.tab", "r")
    fun_names_file = open(params["cog_dir"]+"/COG/fun2003-2014.tab")
    cog_names_file = open(params["cog_dir"]+"/COG/cognames2003-2014.tab")

    cdd_to_cog = []
    for line in cog2cdd_file:
        tokens = line.split("\t")
        cdd_to_cog.append( (int(tokens[1].strip()), int(tokens[0][3:])) )
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

parser.add_argument("--db_name", nargs="+", default=DEFAULT_DB_NAME,
        help="name of the database (default name George)")

parser.add_argument("--load_cog", nargs="?", default=True,
        help="load cog definitions (default yes)")

parser.add_argument("--cog_dir", nargs="?", default="./",
        help="directory where the cog definitions files are (default current directory)")

parser.add_argument("--load_kegg", nargs="?", default=True,
        help="load kegg definitions (default yes)")

parser.add_argument("--db_type", nargs="+", default=DEFAULT_TYPE,
        help="database type (either sqlite or mysql)")

parser.add_argument("--db_psswd", nargs="+", default=DEFAULT_PSSWD,
        help="set db password (default none)")

args = parser.parse_args()


db = setup_biodb(args)
create_data_table(db, args)
if args.load_cog:
    setup_cog(db, args)

if args.load_kegg:
    load_KO_references(db, args)
