#!/usr/bin/env python

import os, sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import KO_module
import db_utils


# to be manually set for tests
database_settings = {}


params = {"chlamdb.db_name": "George", "chlamdb.db_type": "sqlite", "chlamdb.db_psswd" : ""}
database = db_utils.DB.load_db("../George", params)

definitions = database.get_all_modules_definition()

for KO_id, definition in definitions:
    try:
        mod = KO_module.ModuleParser(definition)
        mod.parse()
    except Exception as e:
        print(str(e))
        print(f"Failed on KO{KO_id}" + definition)
