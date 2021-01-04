# {{{1 Preamble

# {{{2 Imports

from lib.snake import (
    alias_recipe,
    noperiod_wc,
    single_param_wc,
    limit_numpy_procs_to_1,
    curl_recipe,
    limit_numpy_procs,
)
from lib.util import idxwhere
import pandas as pd
import math
from itertools import product

# {{{1 Configuration

# {{{2 General Configuration


container: "docker://bsmith89/compbio@sha256:5755ae1efa04bbb4f7994fbcab55d3d386a35e44d4546c8a1c01c7a62d7270cb"


configfile: "config.yaml"
configfile: "config_local.yaml"


if "MAX_THREADS" in config:
    MAX_THREADS = config["MAX_THREADS"]
else:
    MAX_THREADS = 99
if "USE_CUDA" in config:
    USE_CUDA = config["USE_CUDA"]
else:
    USE_CUDA = 0

# {{{2 Data Configuration

_mgen_library = pd.read_table(config["_mgen_library"], index_col="mgen_library_id")
config["mgen_library"] = {}
for mgen_library_id, row in _mgen_library.iterrows():
    config["mgen_library"][mgen_library_id] = {}
    config["mgen_library"][mgen_library_id]["r1"] = row["filename_r1"]
    config["mgen_library"][mgen_library_id]["r2"] = row["filename_r2"]

_mgen_library_group = pd.read_table(config["_mgen_library_group"])
config["mgen_library_group"] = {}
for mgen_library_group, d in _mgen_library_group.groupby("mgen_library_group"):
    config["mgen_library_group"][mgen_library_group] = d.mgen_library_id.tolist()


# {{{2 Sub-pipelines


include: "snake/local.smk"
include: "snake/mgen_taxonomic.smk"
include: "snake/mgen_functional.smk"
include: "snake/mgen_preprocess.smk"
include: "snake/rrs_taxonomic.smk"
include: "snake/util.smk"
include: "snake/general.smk"


wildcard_constraints:
    r="r1|r2|r3",
    group=noperiod_wc,
    library=noperiod_wc,
    species=noperiod_wc,
    strain=noperiod_wc,
    compound_id=noperiod_wc,
    hmm_cutoff="XX",
    model=noperiod_wc,
    param=single_param_wc,
    params=noperiod_wc,


# {{{2 Default actions


rule all:
    output:
        ["sdata/database.db"],


# {{{1 Database


rule build_db:
    output:
        "sdata/database.db",
    input:
        script="scripts/build_db.py",
        schema="schema.sql",
        subject="smeta/subject.tsv",
        visit="meta/visit.tsv",
        sample="meta/sample.tsv",
        rrs_library="meta/rrs_library.tsv",
        mgen_library="meta/mgen_library.tsv",
        mgen_library_group="meta/mgen_library_group.tsv",
        rrs_otu="data/otu_details.tsv",
        rrs_otu_count="data/otu_count.stack.tsv",
        genome="ref/uhgg_genome_to_taxonomy.tsv",
        # species='ref/iggtools_species_to_uhgg_rep_genome.tsv',
        # mgen_species_coverage='data/core.iggtools.species.tsv',
        # mgen_strain_coverage='data/core.uhgg_phyeco-sf20-strains.tsv',
        metab_chemical="meta/metab_chemical.tsv",
        metab_sample="meta/metab_sample.tsv",
        bile_acid="meta/bile_acid.tsv",
        metab_peak_size="data/core.metab.peak_size.tsv",
        mgen_ko="data/core.uhgp90-blastqx.ko-tally.merged.tsv",
    shell:
        r"""
        {input.script} {output} {input.schema} \
                _subject:{input.subject}:1 \
                visit:{input.visit}:1 \
                _sample:{input.sample}:1 \
                rrs_library:{input.rrs_library}:1 \
                rrs_otu:{input.rrs_otu}:1 \
                rrs_otu_count:{input.rrs_otu_count}:1 \
                genome:{input.genome}:0 \
                mgen_library:{input.mgen_library}:1 \
                mgen_library_group:{input.mgen_library_group}:1 \
                metab_chemical:{input.metab_chemical}:1 \
                metab_sample:{input.metab_sample}:1 \
                bile_acid:{input.bile_acid}:1 \
                metab_peak_size:{input.metab_peak_size}:1 \
                mgen_ko:{input.mgen_ko}:0 \

        """
