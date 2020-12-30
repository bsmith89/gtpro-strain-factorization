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

# config['library_group'] = {}
# _library_group = pd.read_table(
#     config['_library_group'], index_col=['library_id'], squeeze=True
# )
# for library_group in _library_group.unique():
#     config['library_group'][library_group] = idxwhere(_library_group == library_group)

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

# {{{1 Process data
