# {{{1 Preamble

# {{{2 Imports

from lib.snake import (
    alias_recipe,
    noperiod_wc,
    integer_wc,
    single_param_wc,
    limit_numpy_procs_to_1,
    curl_recipe,
    limit_numpy_procs,
)
from lib.pandas import idxwhere
import pandas as pd
import math
from itertools import product
from textwrap import dedent as dd
import os.path as path

# {{{1 Configuration

# {{{2 General Configuration


configfile: "config.yaml"

local_config_path = "config_local.yaml"
if path.exists(local_config_path):
    configfile: local_config_path


if "container" in config:

    container: config["container"]


if "MAX_THREADS" in config:
    MAX_THREADS = config["MAX_THREADS"]
else:
    MAX_THREADS = 99
if "USE_CUDA" in config:
    USE_CUDA = config["USE_CUDA"]
else:
    USE_CUDA = 0

# {{{2 Data Configuration

if "_mgen" in config and path.exists(config["_mgen"]):
    _mgen = pd.read_table(config["_mgen"], index_col="mgen_id")
    config["mgen"] = {}
    for mgen_id, row in _mgen.iterrows():
        config["mgen"][mgen_id] = {}
        config["mgen"][mgen_id]["r1"] = row["filename_r1"]
        config["mgen"][mgen_id]["r2"] = row["filename_r2"]
else:
    warn(dd("""
            Could not load config from config["_mgen"].
            Check that path is defined and file exists.
            """))
    config["mgen"] = {}

if "_mgen_x_mgen_group" in config and path.exists(config["_mgen_x_mgen_group"]):
    _mgen_x_mgen_group = pd.read_table(config["_mgen_x_mgen_group"])
    config["mgen_group"] = {}
    for mgen_group, d in _mgen_x_mgen_group.groupby("mgen_group"):
        config["mgen_group"][mgen_group] = d.mgen_id.tolist()
else:
    warn(dd("""
            Could not load config from config["_mgen_x_mgen_group"].
            Check that path is defined and file exists.
            """))
    config["mgen_group"] = {}



# {{{2 Sub-pipelines


include: "snake/template.smk"
include: "snake/util.smk"
include: "snake/general.smk"
if path.exists("snake/local.smk"):
    include: "snake/local.smk"


wildcard_constraints:
    r="r|r1|r2|r3",
    group=noperiod_wc,
    mgen=noperiod_wc,
    species=noperiod_wc,
    strain=noperiod_wc,
    compound_id=noperiod_wc,
    hmm_cutoff="XX",
    model=noperiod_wc,
    param=single_param_wc,
    params=noperiod_wc,


# {{{1 Default actions


rule all:
    input:
        ["sdata/database.db"],


# {{{1 Checkpoint rules


rule gather_all_mgen_read_pairs_from_mgen_group:
    output:
        touch("data/{group}.a.{stem}.ALL_MGEN_PAIRS.flag"),
    input:
        r1=lambda w: [
            f"data/{mgen}.r1.{{stem}}" for mgen in config["mgen_group"][w.group]
        ],
        r2=lambda w: [
            f"data/{mgen}.r2.{{stem}}" for mgen in config["mgen_group"][w.group]
        ],
    shell:
        "touch {output}"


localrules:
    gather_all_mgen_read_pairs_from_mgen_group,


rule gather_all_mgen_from_mgen_group:
    output:
        touch("data/{group}.a.{stem}.ALL_MGEN.flag"),
    input:
        lambda w: [f"data/{mgen}.r.{{stem}}" for mgen in config["mgen_group"][w.group]],
    shell:
        "touch {output}"


localrules:
    gather_all_mgen_from_mgen_group,


# {{{1 Database


database_inputs = [
    # Metadata
    DatabaseInput("subject", "smeta/subject.tsv", True),
    DatabaseInput("sample", "meta/sample.tsv", True),
    DatabaseInput("mgen", "meta/mgen.tsv", True),
    DatabaseInput("mgen_x_mgen_group", "meta/mgen_x_mgen_group.tsv", True),
]


rule build_db:
    output:
        "sdata/database.db",
    input:
        script="scripts/build_db.py",
        schema="schema.sql",
        inputs=[entry.path for entry in database_inputs],
    params:
        args=[entry.to_arg() for entry in database_inputs],
    shell:
        dd(
            r"""
        {input.script} {output} {input.schema} {params.args}
        """
        )
