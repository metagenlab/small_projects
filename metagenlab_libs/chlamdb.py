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

import sqlite3
import os
import MySQLdb

DEFAULT_PSSWD = ""
DEFAULT_TYPE = "sqlite"
DEFAULT_DB_NAME = "George"

def create_data_table(kwargs):
    db_type = kwargs.get("chlamdb.db_type", DEFAULT_TYPE)
    db_name = kwargs.get("chlamdb.db_name", DEFAULT_DB_NAME)
    if db_type=="sqlite":
        conn = sqlite3.connect(db_name)
        cursor = conn.cursor()
    else:
        sqlpsw = os.environ['SQLPSW']

        conn = MySQLdb.connect(host="localhost", # your host, usually localhost
                                    user="root", # your username
                                    passwd=sqlpsw, # your password
                                    db=db_name) # name of the data base
        cursor = conn.cursor()
    entry_list = [
        # Done
        ("orthology", "mandatory", False),

        # Done
        ("orthogroup_alignments", "mandatory", False),

        ("old_locus_table", "mandatory", False),

        # Done
        ("reference_phylogeny", "mandatory", False),

        # Done
        ("taxonomy_table", "mandatory", False),

        # Done
        ("genome_statistics", "mandatory", False),

        ############# Optional ###################
        ("BLAST_database", "optional", False),

        # Done
        ("gene_phylogenies", "optional", False),

        ("interpro_data", "optional", False),
        ("interpro_comparative", "optional", False),
        ("interpro_comparative_accession", "optional", False),
        ("priam_data", "optional", False),
        ("priam_comparative", "optional", False),
        ("priam_comparative_accession", "optional", False),

        # Done: will need to rewrite the queries
        ("COG", "optional", False),

        # Done, to be tested
        ("KEGG", "optional", False),

        ("pfam_comparative", "optional", False),
        ("pfam_comparative_accession", "optional", False),       
        ("TCDB_data", "optional", False),
        ("psortb_data", "optional", False),
        ("T3SS_data", "optional", False),
        ("PDB_data", "optional", False),
        ("BLAST_refseq", "optional", False),
        ("BLAST_swissprot", "optional", False),

        # Done
        ("BBH_phylogenies", "optional", False),
        ("GC_statistics", "optional", False),
        ("gene_clusters", "optional", False),
        ("phylogenetic_profile", "optional", False),
        ("synonymous_table", "optional", False),
        ("interpro_taxonomy", "optional", False), # interpro taxnonomy statistics
        ("pfam_taxonomy", "optional", False), #  taxnonomy statistics
        ("COG_taxonomy", "optional", False) # COG taxnonomy statistics
    ]
    
    sql = (
        "CREATE TABLE biodb_config"
        "(name varchar(200), type varchar(200), status BOOLEAN);"
    )
    
    cursor.execute(sql)
    conn.commit()
    
    sql = 'insert into biodb_config values ("%s", "%s", %s)'
    for row in entry_list:
        cursor.execute(sql % (row[0], row[1], row[2]),)
    conn.commit()

def setup_biodb(kwargs):
    sqlpsw = kwargs.get("chlamdb.sql_psswd", DEFAULT_PSSWD)
    db_type = kwargs.get("chlamdb.db_type", DEFAULT_TYPE)
    db_name = kwargs.get("chlamdb.db_name", DEFAULT_DB_NAME)
    schema_dir = kwargs.get("chlamdb.biosql_schema_dir", "biosql_scheme")

    if db_type=="sqlite":
        conn = sqlite3.connect(db_name)
        cursor = conn.cursor()
        url_biosql_scheme = 'biosqldb-sqlite.sql'
        err_code = os.system(f"sqlite3 {db_name} < {schema_dir}/{url_biosql_scheme}")
        conn.execute("pragma journal_mode=wal")
    else:
        # NOTE: this part is untested!
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
