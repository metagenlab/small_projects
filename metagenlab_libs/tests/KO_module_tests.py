#!/usr/bin/env python
#
# Tests several things: 
# 1) that all module definitions in the database can be parsed (which seems
#    like a good sanity check before launching the database in production)
# 2) some unit tests for specific expression to ensure that the computation
#    of n_missing is actually correct     

import os, sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import KO_module
import db_utils


# to be manually set for tests
database_settings = {}
params = {"chlamdb.db_name": "George", "chlamdb.db_type": "sqlite", "chlamdb.db_psswd" : ""}
database = db_utils.DB.load_db("../George", params)

definitions = database.get_all_modules_definition()

passed_parsing = True
for KO_id, definition in definitions:
    try:
        mod = KO_module.ModuleParser(definition)
        mod.parse()
    except Exception as e:
        print(str(e))
        print(f"Failed on KO{KO_id}" + definition)
        passed_parsing = False

if not passed_parsing:
    print("Could not parse some modules def in the database")
    print("Stopped testing here. Please fix parsing.")

def test_expr(expression, expected_result, kos):
    parser = KO_module.ModuleParser(expression)
    tree = parser.parse()
    n_missing = tree.get_n_missing(kos)
    if n_missing == expected_result:
        print("Passed: ", expression)
    else:
        print("Failed: ", expression, " expected ", expected_result, " got ", n_missing)

# simple and
test_expr("K00001 K00002", 0, {1, 2})
test_expr("K00001 K00002", 1, {2})
test_expr("K00001 K00002", 1, {1})
test_expr("K00001 K00002", 2, {})

# simple or
test_expr("K00001,K00002", 0, {1, 2})
test_expr("K00001,K00002", 0, {2})
test_expr("K00001,K00002", 0, {2})
test_expr("K00001,K00002", 1, {})

# simple complex
test_expr("K00001+K00002", 0, {1, 2})
test_expr("K00001+K00002", 1, {2})
test_expr("K00001+K00002", 1, {1})
test_expr("K00001+K00002", 2, {})

# or / and + parentheses combination
test_expr("(K00001,K00002) K00003", 0, {1, 2, 3})
test_expr("(K00001,K00002) K00003", 0, {2, 3})
test_expr("(K00001,K00002) K00003", 0, {1, 3})
test_expr("(K00001,K00002) K00003", 1, {1, 2})
test_expr("(K00001,K00002) K00003", 1, {3})
test_expr("(K00001,K00002) K00003", 2, {})
