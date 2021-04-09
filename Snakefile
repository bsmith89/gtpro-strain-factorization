# {{{1 Preamble

# {{{2 Imports

from lib.snake import alias_recipe
from lib.pandas import idxwhere
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

container: 'docker://bsmith89/compbio@sha256:82e4cd0171e7b8ff9851108aeb625abfce6a0dcb1d3364f6af58c339ec560c6c'

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
    params:
        port=config['jupyter_port']
    shell: 'jupyter lab --port={params.port}'

rule start_jupyter_notebook:
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

rule qsub_gpu_test:
    output: 'qsub_gpu_test.{n}.env'
    threads: 4
    resources:
        gpu=1,
        gpu_mem=10000,
    shell:
        """
        env > {output}
        """


# {{{1 Generalizable rules

# {{{1 Database

# {{{1 Download and organize reference data

# {{{1 Organize raw data

# # This rule is used to configure the project based on the contents of
# # the GT-PRO output dir.
# rule extract_read_metadata_from_dir:
#     output: 'data/meta.reads.{r}.tsv'
#     input: 'raw/gtpro_output/'
#     params: r=lambda w: w.r
#     shell:
#         r"""
# find raw/gtpro_output/ -name '*_{params.r}.parse.tsv.bz2' \
#         | sed 's:raw/gtpro_output/\(.*\)_{params.r}.parse.tsv.bz2:\1\t\1_{params.r}:' \
#     > {output}
#         """



# rule extract_libraries_into_library_group_table:
#     output: 'data/meta.library_group.tsv'
#     input: 'meta/library.tsv'
#     shell:
#         """
# awk -v OFS='\t' '$2!="" && $3!="" {{print $1,"core"}}' < {input} > {output}
#         """


# {{{1 Process data

rule gtpro_processed_reads:
    output: 'data/gtpro/{stem}.m.proc.r{r}.gtpro.gz'
    input: 'data/reads/{stem}.m.proc.r{r}.fq.gz'
    params:
        db_l=32,
        db_m=36,
        db_stem=config['gtpro_db_src'] + config['gtpro_db_stem']
    threads: 4
    resources:
        mem_mb=60000,
        pmem=60000 // 4,
        walltime_hr=4,
    shell:
        """
        gt_pro -t {threads} -l {params.db_l} -m {params.db_m} -d {params.db_stem} -f -o {output} {input}
        mv {output}.tsv.gz {output}
        """

# # NOTE: Comment-out to speed up execution.
# rule gtpro_finish_processing_reads:
#     output: 'data/gtpro/{library}_{r}.parse.tsv.bz2'
#     input:
#         script='scripts/gtp_parse.py',
#         gtpro='data/gtpro/{library}.m.proc.r{r}.gtpro.gz',
#     params:
#         db=config['gtpro_db_src'] + 'variants_main.covered.hq.snp_dict.tsv',
#     shell:
#         """
#         {input.script} --dict {params.db} --in <(zcat {input.gtpro}) --v2 \
#                 | bzip2 -c \
#             > {output}
#         """

# # NOTE: This rule should be run for everything *except* UCFMT libs.
# # NOTE: It's much faster to run this process manually e.g.
# #       `find raw/gtpro_output/ -name '*.parse.tsv.bz2' -exec ln -rfs {} data/gtpro/ \;`
# # NOTE: Comment out this rule to speed up execution.
# rule link_precomputed_libs:
#     output: 'data/gtpro/{library}_{r}.parse.tsv.bz2'
#     input: 'raw/gtpro_output/{library}_{r}.parse.tsv.bz2'
#     shell: alias_recipe
# localrules: link_precomputed_libs

# Helper rule that pre-formats paths from library_id to r1 and r2 paths.
rule helper_build_library_path_table:
    output: temp('data/{library_group}.library_to_gtpro_path.tsv'),
    run:
        with open(output[0], 'w') as f:
            for library in config['library_group'][wildcards.library_group]:
                print(
                        library,
                        f'data/gtpro/{library}_1.parse.tsv.bz2',
                        f'data/gtpro/{library}_2.parse.tsv.bz2',
                        sep='\t',
                        file=f
                    )

rule count_species_lines_from_one_read:
    output: 'data/{library_group}.gtpro-site_tally.r{r}.tsv.bz2',
    input:
        script='scripts/tally_gtpro_species_lines.sh',
        gtpro=lambda w: [f'data/gtpro/{library}_{w.r}.parse.tsv.bz2'
                         for library in config['library_group'][w.library_group]],
        argsfile='data/{library_group}.library_to_gtpro_path.tsv',
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

# NOTE: Comment out this rule to speed up DAG evaluation.
# Selects a single species from every file and concatenates.
rule concatenate_all_libraries_one_read_count_data_from_one_species:
    output: 'data/{library_group}.sp-{species}.gtpro-pileup.read_r{r}.tsv.bz2'
    input:
        script='scripts/select_gtpro_species_lines.sh',
        gtpro=lambda w: [f'data/gtpro/{library}_{w.r}.parse.tsv.bz2'
                         for library in config['library_group'][w.library_group]],
        argsfile='data/{library_group}.library_to_gtpro_path.tsv',
    params:
        path_column=lambda w: int(w.r) + 1,
        species=lambda w: w.species,
    threads: MAX_THREADS
    resources:
        walltime_hr=6,
    shell:
        """
        parallel --colsep='\t' --bar -j {threads} \
                {input.script} {params.species} :::: <(cut -f1,{params.path_column} {input.argsfile}) \
            | bzip2 -c \
            > {output}

        """

rule merge_both_reads_species_count_data:
    output: 'data/{library_group}.sp-{species}.gtpro-pileup.nc'
    input:
        script='scripts/merge_gtpro_to_netcdf.py',
        r1='data/{library_group}.sp-{species}.gtpro-pileup.read_r1.tsv.bz2',
        r2='data/{library_group}.sp-{species}.gtpro-pileup.read_r2.tsv.bz2',
    threads: 4
    resources:
        walltime_hr=4,
        mem_mb=100000,
        pmem=lambda w, threads: 100000 // threads,
    shell:
        """
        {input.script} {input.r1} {input.r2} {output}
        """

rule append_ucfmt_pileup:
    output: 'data/{group}_and_ucfmt.sp-{species}.gtpro-pileup.nc'
    input:
        script='scripts/concatenate_gtpro_pileups.py',
        pileups=[
            'data/{group}.sp-{species}.gtpro-pileup.nc',
            'data/ucfmt.sp-{species}.gtpro-pileup.nc',
        ],
    shell:
        """
        {input.script} {output} {input.pileups}
        """
ruleorder: append_ucfmt_pileup > merge_both_reads_species_count_data

rule filter_pileup_positions:
    output: '{stem}.gtpro-pileup.filt.nc'
    input:
        script='scripts/filter_gtpro_pileup_positions.py',
        pileup='{stem}.gtpro-pileup.nc',
    params:
        minor_allele_thresh=0.02,
        position_thresh=0.35,
        npos_subsample=5000,
        dist_thresh=0.10,
        clust_size_thresh=3,
        clust_pos_frac_thresh=0.51,
        frac_clust_thresh=0.95,
    shell:
        r"""
        {input.script} {input.pileup} \
                {params.minor_allele_thresh} \
                {params.position_thresh} \
                {params.npos_subsample} \
                {params.dist_thresh} \
                {params.clust_size_thresh} \
                {params.clust_pos_frac_thresh} \
                {params.frac_clust_thresh} \
                {output}
        """

def estimate_strain_facts_gpu_mem_mb(n, g, s, constant=64):
    return int(((n * s) + 2 * (n * g) + (g * s)) * constant / 1e6)

rule factorize_strains:
    output: 'data/{group}.sp-{species}.gtpro-pileup.filt.sfacts-s{nstrain}-g{npos}-gamma{gamma_hyper}-rho{rho_hyper}-pi{pi_hyper}-eps{epsilon}-alph{alpha}.nc'
    input:
        script='scripts/strain_facts.py',
        pileup='data/{group}.sp-{species}.gtpro-pileup.filt.nc',
    params:
        nstrain=lambda w: int(w.nstrain),
        npos=lambda w: int(w.npos),
        gamma_hyper=lambda w: 10**(-int(w.gamma_hyper)),
        rho_hyper=lambda w: 10**(-int(w.rho_hyper)),
        pi_hyper=lambda w: 10**(-int(w.pi_hyper)),
        epsilon=lambda w: float(w.epsilon) * 1e-5,
        alpha=lambda w: int(w.alpha),
        cvrg_thresh=0.05,
        device={0: 'cpu', 1: 'cuda'}[USE_CUDA],
        learning_rate=0.050,
        stop_at=10,
        max_iter=int(1e4),
        collapse=0.005,
    threads: 8
    resources:
        gpu={0: 0, 1: 1}[USE_CUDA],
        gpu_mem_mb=lambda w: estimate_strain_facts_gpu_mem_mb(10000, int(w.npos), int(w.nstrain)),
        # n=10000 is a guess of the maximum number of libs.
        walltime_hr=3,
        mem_mb=30000,
        pmem=lambda w, threads: 30000 // threads,
    shell:
        limit_numpy_procs + r"""
        scripts/strain_facts.py \
                --cvrg-thresh 0.05 \
                --nstrain {params.nstrain} \
                --npos {params.npos} \
                --gamma-hyper {params.gamma_hyper} \
                --rho-hyper {params.rho_hyper} \
                --pi-hyper {params.pi_hyper} \
                --epsilon {params.epsilon} \
                --alpha {params.alpha} \
                --collapse {params.collapse} \
                --device {params.device} \
                --learning-rate {params.learning_rate} \
                --stop-at {params.stop_at} \
                --max-iter {params.max_iter} \
                --outpath {output} \
                {input.pileup}

        """
