#!/usr/bin/env python3

import loompy
import random


loom_file = sys.argv[1]
n         = sys.argv[2]

with loompy.connect(loom_file) as ds:
	random.rows = random.sample(range(ds.shape[1]), n)
	random.rows.sort()
	matrix <- ds[:,random.rows]