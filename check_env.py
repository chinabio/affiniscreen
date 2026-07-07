#!/usr/bin/env python
"""Convenience wrapper:  python check_env.py [--write-template]

Compatible with Python 2.7 / 3.5+ so old interpreters still get a useful
diagnostic instead of a stack trace.
"""

# -*- coding: utf-8 -*-
# Author:    AffiniScreen contributors
# Maintainer: AffiniScreen team
# Contact:   (see repository)
# Part of AffiniScreen.

from __future__ import print_function
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from amber_md.env import main
sys.exit(main(sys.argv[1:]))
