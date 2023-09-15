#!/usr/bin/env python3

import sys
import pandas as pd

pd.read_parquet(sys.argv[1], engine="pyarrow").to_csv("/dev/stdout", sep="\t", index=False)
