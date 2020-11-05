#!/usr/bin/env python

import os, sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import KO_module


test_string = "K18294 K18295+K18296-K08721"
foo = KO_module.Tokenizer(test_string)

for token in foo:
    print(token)
