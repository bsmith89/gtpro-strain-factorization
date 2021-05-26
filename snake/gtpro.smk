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
