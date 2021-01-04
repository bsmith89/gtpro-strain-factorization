#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sqlite3
import pandas as pd
import sys

if __name__ == "__main__":
    con = sqlite3.connect(sys.argv[1])
    with open(sys.argv[2]) as f:
        con.executescript(f.read())

    for input_string in sys.argv[3:]:
        table_name, input_path, header = input_string.split(":")
        header = bool(int(header))

        template = pd.read_sql(f"SELECT * FROM {table_name} LIMIT 0", con=con)
        columns = template.columns.to_list()

        print(f"{table_name} : has_header={header}")
        d = pd.read_csv(
            input_path,
            sep="\t",
            skiprows={True: 1, False: 0}[header],
            names=columns,
        )
        print(d.info())
        d.to_sql(table_name, con=con, if_exists="append", index=False)
    con.execute("PRAGMA foreign_keys = TRUE;")
    con.close()
