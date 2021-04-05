#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sqlite3
import pandas as pd
import sys
from dataclasses import dataclass, field
from typing import Any, Mapping
from lib.util import info


@dataclass
class DatabaseInput:
    table_name: str
    path: str
    has_header: bool
    constant_fields: Mapping[str, str] = field(default_factory=dict)

    def _constant_fields_to_args_list(self):
        return [
            f"{entry[0]}={entry[1]}" for entry in self.constant_fields.items()
        ]

    def to_arg(self) -> str:
        return ":".join(
            [self.table_name, self.path, str(int(self.has_header))]
            + self._constant_fields_to_args_list()
        )

    @classmethod
    def from_arg(cls, arg_string):
        table_name, path, has_header, *constant_args = arg_string.split(":")
        has_header = bool(int(has_header))
        constant_fields = {}
        for arg in constant_args:
            k, v = arg.split("=")
            constant_fields[k] = v
        return cls(table_name, path, has_header, constant_fields)


if __name__ == "__main__":
    db_path, script_path, *input_args = sys.argv[1:]
    con = sqlite3.connect(db_path)
    with open(script_path) as f:
        con.executescript(f.read())

    for arg in input_args:
        db_input = DatabaseInput.from_arg(arg)
        info(db_input)
        template = pd.read_sql(
            f"SELECT * FROM {db_input.table_name} LIMIT 0", con=con
        )
        columns = template.columns.to_list()

        d = pd.read_csv(
            db_input.path,
            sep="\t",
            skiprows={True: 1, False: 0}[db_input.has_header],
            names=columns,
        )
        d = d.assign(**db_input.constant_fields)
        info(d.info())
        d.to_sql(db_input.table_name, con=con, if_exists="append", index=False)
    con.execute("PRAGMA foreign_keys = TRUE;")
    con.close()
