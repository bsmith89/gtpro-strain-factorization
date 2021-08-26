rule start_jupyter:
    threads: MAX_THREADS
    params:
        port=config["jupyter_port"],
    shell:
        "jupyter lab --port={params.port}"


rule start_ipython:
    threads: MAX_THREADS
    shell:
        limit_numpy_procs + "ipython"


rule start_shell:
    shell:
        "bash"


rule visualize_rulegraph:
    output:
        "data/rulegraph.dot",
    input:
        "Snakefile",
    shell:
        dd(
            """
        snakemake --rulegraph all > {output}
        """
        )


rule generate_report:
    output:
        "fig/report.html",
    input:
        "Snakefile",
    shell:
        dd(
            """
        snakemake --forceall --report {output} all
        """
        )


rule dot_to_pdf:
    output:
        "fig/{stem}.pdf",
    input:
        "data/{stem}.dot",
    shell:
        dd(
            """
        dot -Tpdf < {input} > {output}
        """
        )


rule processed_notebook_to_html:
    output:
        "build/{stem}.ipynb.html",
    input:
        "build/{stem}.ipynb",
    shell:
        dd(
            """
        jupyter nbconvert -t html {input} {output}
        """
        )


rule query_db:
    output:
        "data/{db}.select_{query}.tsv",
    input:
        db="data/{db}.db",
        query="scripts/query/{query}.sql",
    shell:
        dd(
            """
        sqlite3 -header -separator '\t' {input.db} < {input.query} > {output}
        """
        )


rule config_debug:
    output:
        "config_debug.{config_key}",
    params:
        meta=lambda w: config[w.config_key],
    shell:
        """
        echo {wildcards.config_key}
        echo {params.meta}
        false  # Recipe never succeeds.
        """
