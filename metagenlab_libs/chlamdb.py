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
###############################################################################
# NOTE: for this script to work, you either need to disable multithreading by 
# setting N_THREAD to 1 or you need to patch REST.py to use urllib3 (thread safe).
#
#
# TODO: as the REST API does not seem to be under active development, it would
# be a solution to just import the source files into this directory to solve this
# problem.
###############################################################################
#
# Bastian Marquis (bastian.marquis@protonmail.com)
# Date: 29.10.2020

import os
import db_utils
import argparse

import queue
import threading

# to be removed in favor of a local version
from Bio.KEGG import REST

DEFAULT_PSSWD = ""
DEFAULT_TYPE = "sqlite"
DEFAULT_DB_NAME = "George"

# from REST documentation, can get a max of 10 queries 
# in kegg_get
MAX_N_QUERIES = 10

# number of threads that will be used simultaneously
# to download the ko genes. 5 seems a good compromise
# between being blacklisted and speed
DEFAULT_N_THREADS = 5
DEFAULT_KO_DIR = "tmp"

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
    db = db_utils.DB.load_db(db_name, chlamdb_args)
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
    def __init__(self):
        self.entry = None
        self.descr = None
        self.classes = []
        self.definition = []

    def simplified_entry(self):
        return int(self.entry[len("M"):])

    def get_definition(self):
        if len(self.definition) == 1:
            return self.definition[0]
        return ",".join(f"({part})" for part in self.definition)

    def sub_category(self):
        return self.classes[-1]

    def category(self):
        return self.classes[-2]

    def is_signature(self):
        return self.classes[0] == "Signature modules"

    def __str__(self):
        return self.entry + ": " + self.descr

    def __hash__(self):
        return hash(self.entry)

    def __eq__(self, other):
        return self.entry == other.entry


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


def parse_module(handle):
    module = Module()
    for line in handle:
        if line[:3] == "///":
            yield module
            module = Module()
            continue
        if line[:12] != "            ":
            keyword = line[:12].strip()
        data = line[12:].strip()
        if keyword == "ENTRY":
            module.entry = data.split()[0]
        if keyword == "NAME":
            module.descr = data
        if keyword == "CLASS":
            tokens = data.split(";")
            module.classes = [token.strip() for token in tokens]
        if keyword == "DEFINITION":
            module.definition.append(data)


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
            record.modules.append(module)
        elif keyword == "PATHWAY     ":
            pathway, *descr_tokens = data.split()
            descr = (" ".join(descr_tokens)).strip()
            record.pathways.append(Pathway(pathway, descr))


def download_ko_genes(gene_queue, tmp_dir):
    done = False
    while not done:
        queries = []
        try:
            # non-blocking call, will raise an Empty exception
            # once there is no more elements to download
            for i in range(MAX_N_QUERIES):
                gene = gene_queue.get(False)
                queries.append(gene)
                size = gene_queue.qsize()
                if size % 100 == 0:
                    print(size)
        except Exception as e:
            if len(queries) == 0:
                break
            done = True

        try:
            buff = REST.kegg_get(queries)
            buff.reconfigure(encoding="Latin-1")
            text_version = buff.read()
            output_file = open(tmp_dir+"/"+queries[0], "w")
            output_file.write(text_version)
            output_file.close()
        except Exception as e:
            print(e)
            print(f"ERROR: Could not download one of the following kegg ")
            print(",".join(queries))


def download_ko_files(ko_dir, n_threads=DEFAULT_N_THREADS):
    ko_string = REST.kegg_list("KO")
    already_downloaded = []
    if not os.path.exists(ko_dir):
        os.mkdir(ko_dir)
    else:
        for ko in os.listdir(ko_dir):
            ko_file = open(ko_dir+"/"+ko, "r")
            for gene in parse_gene(ko_file):
                already_downloaded.append("ko:"+gene.entry)

    if len(already_downloaded)>0:
        print(f"Already {len(already_downloaded)} downloaded")

    gene_queue = queue.Queue()
    gene_list = []
    for line in ko_string:
        gene = line.split()[0]
        if gene in already_downloaded:
            continue
        gene_queue.put(gene)
        gene_list.append(gene)

    if len(gene_list) > 0:
        print(f"Preparing to download : {len(gene_list)}")
    else:
        print("Everything already downloaded")
        return

    threads = []
    for i in range(n_threads):
        thread = threading.Thread(target=download_ko_genes, args=(gene_queue, ko_dir))
        threads.append(thread)
        thread.start()

    # wait for completion
    for thread in threads:
        thread.join()

def chunk_modules(modules, size=MAX_N_QUERIES):
    curr = []
    for i, module in enumerate(modules):
        curr.append(module)
        if i % MAX_N_QUERIES == 0:
            yield curr
            curr = []
    if len(curr)>0:
        yield curr

def load_KO_references(db, params, ko_dir=DEFAULT_KO_DIR):
    genes = []

    for ko in os.listdir(ko_dir):
        ko_file = open(ko_dir+"/"+ko, "r")
        for gene in parse_gene(ko_file):
            genes.append(gene)

    # get only the module present in the genes
    modules_set = {module for gene in genes for module in gene.modules}
    hsh_modules = {}
    print("Downloading module data")
    for i, module_list in enumerate(chunk_modules(modules_set)):
        if i %10==0:
            print(i*MAX_N_QUERIES)
        mod_data = REST.kegg_get(module_list)
        for module in parse_module(mod_data):
            hsh_modules[module.entry] = module
    print("Done")

    pathway_set = set()
    category_set = set()
    for entry, module in hsh_modules.items():
        if len(module.classes) > 0:
            category_set.add(module.sub_category())
            category_set.add(module.category())
    for gene in genes:
        for pathway in gene.pathways:
            pathway_set.add(pathway)

    # really ugly, but spares some web requests
    module_classes = []
    hsh_class_to_id = {}
    for cat_id, cat in enumerate(category_set):
        module_classes.append( (cat_id, cat) )
        hsh_class_to_id[cat] = cat_id
    for entry, module in hsh_modules.items():
        cat_id = hsh_class_to_id[module.category()]
        subcat_id = hsh_class_to_id[module.sub_category()]
        module.cat_id = cat_id
        module.subcat_id = subcat_id

    db.load_ko_module_classes(module_classes)
    db.load_ko_module([(m.simplified_entry(), m.descr, m.get_definition(), m.is_signature(), m.cat_id, m.subcat_id)
        for m in hsh_modules.values()])
    db.load_ko_pathway([(p.simplified_entry(), p.descr) for p in pathway_set])
    db.load_ko_def([(gene.simplified_entry(), gene.definition) for gene in genes])

    ko_to_path = []
    ko_to_module = []
    for gene in genes:
        simp = gene.simplified_entry()
        ko_to_path.extend((simp, path.simplified_entry()) for path in gene.pathways)
        ko_to_module.extend((simp, hsh_modules[mod].simplified_entry())
                for mod in gene.modules)
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
        help="load kegg definitions (default no, must specify ko genes dir)")

parser.add_argument("--db_type", nargs="?", default=DEFAULT_TYPE,
        help="database type (either sqlite or mysql)")

parser.add_argument("--db_psswd", nargs="+", default=DEFAULT_PSSWD,
        help="set db password (default none)")

parser.add_argument("--download_ko_files", action="store_true",
        help="download ko file definition, necessary to do this before --load-kegg")

parser.add_argument("--ko_dir", nargs="?", default=DEFAULT_KO_DIR,
        help=f"ko directory, defaults to {DEFAULT_KO_DIR}")

parser.add_argument("--skip_biodb", action="store_true", default=False,
        help="skip setting up biodb, might be useful if using a pre-existing db")

args = vars(parser.parse_args())

if args.get("download_ko_files", False):
    print("Starting to download ko files to tmp")
    print("Will do this first as this tend to crash")
    print("Just insist until completion")
    ko_dir = args.get("ko_dir")
    download_ko_files(ko_dir)

db = None
if not args.get("skip_biodb"):
    print("Setting up the biosql schema")
    db = setup_biodb(args)
    create_data_table(db)

if db == None:
    db_name = args.get("db_name", DEFAULT_DB_NAME)
    db_type = args.get("db_type", DEFAULT_TYPE)
    db_psswd = args.get("db_psswd", DEFAULT_PSSWD)
    chlamdb_args = { "chlamdb.db_type" : db_type,
            "chlamdb.db_name" : db_name,
            "chlamdb.db_psswd" : db_psswd}
    db = db_utils.DB.load_db(db_name, chlamdb_args)

if args.get("load_cog", False):
    print("Loading COG tables")
    setup_cog(db, args.get("cog_dir", "."))

if args.get("load_kegg", False):
    print("Loading KEGG tables")
    ko_dir = args.get("ko_dir", DEFAULT_KO_DIR)
    genes = load_KO_references(db, args, ko_dir)
