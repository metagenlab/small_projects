#!/usr/bin/env python

import os, sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import KO_module


test_string = "(K01965+K01966,K11263+(K18472,K19312+K22568),K01964+K15036+K15037) K05606 (K01847,K01848+K01849)"

mod = KO_module.ModuleParser(test_string)
mod.parse()
