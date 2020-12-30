# {{{1 Preamble

# {{{2 Imports

from lib.snake import alias_recipe
from lib.util import idxwhere
import pandas as pd

# {{{2 Constants

limit_numpy_procs = \
        """
        export MKL_NUM_THREADS={threads}
        export OPENBLAS_NUM_THREADS={threads}
        export NUMEXPR_NUM_THREADS={threads}
        export OMP_NUM_THREADS={threads}
        export VECLIB_MAXIMUM_THREADS={threads}

        """
limit_numpy_procs_to_1 = limit_numpy_procs.format(threads=1)

independent_theano_compiledir = \
        """
        # Circumvent theano compiledir locking.
        compiledir=$(mktemp -d)
        # tar --strip-components 1 -xzf raw/theano.tgz -C $compiledir
        export THEANO_FLAGS="base_compiledir=$compiledir"
        echo $THEANO_FLAGS

        """

# Utility wildcard constrains
noperiod_wc = '[^.]+'
integer_wc = '[0-9]+'
float_noperiod_wc = '[0-9]+(e[0-9]+)?'
single_param_wc = '[^.-]+'
params_wc = noperiod_wc

wildcard_constraints:
    r='[12]',
    library=noperiod_wc,
    library_group=noperiod_wc,

# {{{1 Configuration

# {{{2 General Configuration

container: 'docker://bsmith89/compbio@sha256:5755ae1efa04bbb4f7994fbcab55d3d386a35e44d4546c8a1c01c7a62d7270cb'

configfile: 'config.yaml'
configfile: 'config_local.yaml'

if 'MAX_THREADS' in config:
    MAX_THREADS = config['MAX_THREADS']
else:
    MAX_THREADS = 99
if 'USE_CUDA' in config:
    USE_CUDA = config['USE_CUDA']
else:
    USE_CUDA = 0

# {{{2 Project Configuration

# FIXME: Refactor this
config['library_group'] = {}
_library_group = pd.read_table(
    config['_library_group'], index_col=['library_id'], squeeze=True
)
for library_group in _library_group.unique():
    config['library_group'][library_group] = idxwhere(_library_group == library_group)

# {{{2 Configure default actions

rule all:
    output: []

# {{{2 Sub-pipelines

include: 'snake/local.snake'

# {{{1 Utility rules

rule initialize_project:
    shell:
        '''
        git config --local filter.dropoutput_ipynb.clean scripts/ipynb_output_filter.py
        git config --local filter.dropoutput_ipynb.smudge cat
        git config --local diff.daff-csv.command "daff.py diff --git"
        git config --local merge.daff-csv.name "daff.py tabular merge"
        git config --local merge.daff-csv.driver "daff.py merge --output %A %O %A %B"
        echo 'Activate your python environment and install packages in `requirements.txt`.'
        '''

rule start_jupyter:
    threads: MAX_THREADS
    params:
        port=config['jupyter_port']
    shell: limit_numpy_procs + 'jupyter notebook --config=nb/jupyter_notebook_config.py --notebook-dir=nb/ --port={params.port}'

rule start_ipython:
    threads: MAX_THREADS
    shell: limit_numpy_procs + 'ipython'

rule start_shell:
    shell: 'bash'

rule visualize_rulegraph:
    output: 'data/rulegraph.dot'
    input: 'Snakefile'
    shell:
        """
        snakemake --rulegraph all > {output}
        """

rule generate_report:
    output: 'fig/report.html'
    input: 'Snakefile'
    shell:
        """
        snakemake --forceall --report {output} all
        """

rule dot_to_pdf:
    output: 'fig/{stem}.pdf'
    input: 'data/{stem}.dot'
    shell:
        """
        dot -Tpdf < {input} > {output}
        """

# {{{1 Generalizable rules

# {{{1 Database

# {{{1 Download and organize reference data

# {{{1 Organize raw data

# # TODO: Convert path appending to python pathlib.
# rule link_raw_gtpro:
#     output: directory('raw/gtpro_output')
#     input: config['raw_mgen_source']
#     shell: 'ln -sT {input} {output}'

# This rule is used to configure the project based on the contents of
# the GT-PRO output dir.
rule extract_read_metadata_from_dir:
    output: 'data/meta.reads.{r}.tsv'
    input: 'raw/gtpro_output/'
    params: r=lambda w: w.r
    shell:
        r"""
find raw/gtpro_output/ -name '*_{params.r}.parse.tsv.bz2' \
        | sed 's:raw/gtpro_output/\(.*\)_{params.r}.parse.tsv.bz2:\1\t\1_{params.r}:' \
    > {output}
        """

rule join_read_metadata_into_library_meta:
    output: 'data/meta.library.tsv'
    input: r1='data/meta.reads.1.tsv', r2='data/meta.reads.2.tsv'
    shell:
        r"""
cat <<EOF | sqlite3 > {output}
.separator \t
CREATE TABLE r1 (
library_id, r1_path
);
.import {input.r1} r1

CREATE TABLE r2 (
library_id, r2_path
);
.import {input.r2} r2

CREATE VIEW library AS
SELECT DISTINCT library_id
FROM (
    SELECT library_id FROM r1
    UNION
    SELECT library_id FROM r2
)
;

SELECT library_id, r1_path, r2_path
FROM library
LEFT JOIN r1 USING (library_id)
LEFT JOIN r2 USING (library_id)
;

EOF
        """

rule extract_libraries_into_library_group_table:
    output: 'data/meta.library_group.tsv'
    input: 'meta/library.tsv'
    shell:
        """
awk -v OFS='\t' '$2!="" && $3!="" {{print $1,"core"}}' < {input} > {output}
        """


# {{{1 Process data

rule gtpro_process_raw_ucfmt_reads:
    output: 'data/{stem}.gtpro.gz'
    input: 'data/{stem}.fq.gz'
    params:
        db_l=32,
        db_m=36,
        db_stem=config['gtpro_dbstem']
    threads: 4
    resources:
        mem_mb=60000,
        pmem=60000 // 4,
    shell:
        """
        zcat {input} \
                | gt_pro -t {threads} -l {params.db_l} -m {params.db_m} -d {params.db_stem} \
                | gzip -c \
            > {output}
        """

# Helper rule that pre-formats paths from library_id to r1 and r2 paths.
rule helper_build_library_path_table:
    output: 'data/{library_group}/library_to_path.tsv',
    run:
        with open(output[0], 'w') as f:
            for library in config['library_group'][wildcards.library_group]:
                print(
                        library,
                        f'raw/gtpro_output/{library}_1.parse.tsv.bz2',
                        f'raw/gtpro_output/{library}_2.parse.tsv.bz2',
                        sep='\t',
                        file=f
                    )

# FIXME: Is there a better way to work with such a long list of filenames?
rule count_species_lines_from_one_read:
    output: 'data/{library_group}/gtpro.species_site_tally.r{r}.tsv.bz2',
    input:
        script='scripts/tally_gtpro_species_lines.sh',
        gtpro=lambda w: [f'raw/gtpro_output/{library}_{w.r}.parse.tsv.bz2'
                         for library in config['library_group'][w.library_group]],
        argsfile='data/{library_group}/library_to_path.tsv',
    params:
        path_column=lambda w: int(w.r) + 1
    threads: MAX_THREADS
    shell:
        r"""
        parallel --colsep='\t' --bar -j {threads} \
                bash {input.script} :::: <(cut -f1,{params.path_column} {input.argsfile}) \
            | bzip2 -c \
            > {output}

        """

# Selects a single species from every file and concatenates.
rule concatenate_all_libraries_one_read_count_data_from_one_species:
    output: 'data/{library_group}/{species}/gtpro.read_r{r}.tsv.bz2'
    input:
        script='scripts/select_gtpro_species_lines.sh',
        gtpro=lambda w: [f'raw/gtpro_output/{library}_{w.r}.parse.tsv.bz2'
                         for library in config['library_group'][w.library_group]],
        argsfile='data/{library_group}/library_to_path.tsv',
    params:
        species=lambda w: w.species,
        path_column=lambda w: int(w.r) + 1
    threads: MAX_THREADS
    shell:
        """
        parallel --colsep='\t' --bar -j {threads} \
                {input.script} {params.species} :::: <(cut -f1,{params.path_column} {input.argsfile}) \
            | bzip2 -c \
            > {output}

        """

rule merge_both_reads_species_count_data:
    output: 'data/{library_group}/{species}/gtpro.nc'
    input:
        script='scripts/merge_gtpro_to_netcdf.py',
        r1='data/{library_group}/{species}/gtpro.read_r1.tsv.bz2',
        r2='data/{library_group}/{species}/gtpro.read_r2.tsv.bz2',
    shell:
        """
        {input.script} {input.r1} {input.r2} {output}
        """

# rule coverage_filter_allele_counts:
#     output: '{stem}.merge_reads.filt.tsv.bz2'
#     input:
#         script='scripts/coverage_filter_allele_counts.py',
#         tsv='{stem}.merge_reads.tsv',
#     shell:
#         """
#         {input.script} {input.tsv} > {output}
#
        # """

rule factorize_strains:
    output: '{stem}gtpro.sfacts.nc'
    input:
        script='scripts/strain_facts.py',
        pileup='{stem}gtpro.nc',
    params:
        device={0: 'cpu', 1: 'cuda'}[config['USE_CUDA']]
    threads: 4
    resources:
        gpu={0: 0, 1: 1}[config['USE_CUDA']],
        vmem={0: 0, 1: 10000}[config['USE_CUDA']],
    shell:
        """
        scripts/strain_facts.py \
                --cvrg-thresh 0.05 \
                --nstrain 3000 \
                --npos 2000 \
                --device {params.device} \
                --learning-rate 5e-2 \
                --stop-at 1 \
                --outpath {output} \
                {input.pileup}

        """
