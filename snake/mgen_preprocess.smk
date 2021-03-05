# {{{1 Download and organize reference data


rule download_illumina_adapters:
    output:
        "raw/ref/illumina_adapters.fn",
    params:
        url="https://raw.githubusercontent.com/vsbuffalo/scythe/master/illumina_adapters.fa",
    resources:
        network_connections=1,
    shell:
        curl_recipe


rule link_illumina_adapters:
    output:
        "ref/illumina_adapters.fn",
    input:
        "raw/ref/illumina_adapters.fn",
    shell:
        alias_recipe


rule download_GRCh38_index:
    output:
        "raw/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz",
    params:
        url="ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz",
    resources:
        network_connections=1,
    shell:
        curl_recipe


rule unpack_GRCh38_index:
    output:
        "raw/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.1.bt2",
        "raw/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.2.bt2",
        "raw/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.3.bt2",
        "raw/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.4.bt2",
        "raw/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.rev.1.bt2",
        "raw/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.rev.2.bt2",
    input:
        "raw/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz",
    shell:
        "tar -C raw/ref/ -xzf {input}"


rule alias_GRCh38_index_file:
    output:
        "ref/GRCh38.{suffix}",
    input:
        "raw/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.{suffix}",
    shell:
        alias_recipe


ruleorder: alias_GRCh38_index_file > bowtie_index_build


# {{{1 Organize raw data


# This rule is intended to enable debugging on platforms without the raw data present.
# It should always fail if actually run.
rule dummy_raw_read_source:
    output:
        "raw/mgen/{raw_mgen_filename}",
    shell:
        dd(
            """
        echo Raw data file {output} does not exist on this computer.  This is a dummy rule.
        false
        """
        )


rule alias_raw_read_r1:
    output:
        "sdata/{mgen}.r1.fq.gz",
    input:
        lambda wildcards: "raw/mgen/{}".format(config["mgen"][wildcards.mgen]["r1"]),
    shell:
        alias_recipe


localrules:
    alias_raw_read_r1,


rule alias_raw_read_r2:
    output:
        "sdata/{mgen}.r2.fq.gz",
    input:
        lambda wildcards: "raw/mgen/{}".format(config["mgen"][wildcards.mgen]["r2"]),
    shell:
        alias_recipe


localrules:
    alias_raw_read_r2,


# {{{1 Process data
# {{{2 Metagenomic reads


rule qc_raw_reads:
    output:
        directory("data/{group}.fastqc.d"),
    input:
        r1=lambda w: [
            f"sdata/{mgen}.r1.fq.gz" for mgen in config["mgen_group"][w.group]
        ],
        r2=lambda w: [
            f"sdata/{mgen}.r2.fq.gz" for mgen in config["mgen_group"][w.group]
        ],
    threads: MAX_THREADS
    shell:
        dd(
            """
        fastqc -t {threads} -o {output} {input}
        """
        )


rule qc_processed_reads:
    output:
        directory("data/{group}.proc.fastqc.d"),
    input:
        r1=lambda w: [
            f"data/{mgen}.r1.proc.fq.gz" for mgen in config["mgen_group"][w.group]
        ],
        r2=lambda w: [
            f"data/{mgen}.r2.proc.fq.gz" for mgen in config["mgen_group"][w.group]
        ],
    threads: MAX_THREADS
    shell:
        dd(
            """
        fastqc -t {threads} -o {output} {input}
        """
        )


rule deduplicate_reads:
    output:
        r1=temp("sdata/{mgen}.r1.dedup.fq.gz"),
        r2=temp("sdata/{mgen}.r2.dedup.fq.gz"),
    input:
        script="scripts/fastuniq_wrapper.sh",
        r1="sdata/{mgen}.r1.fq.gz",
        r2="sdata/{mgen}.r2.fq.gz",
    resources:
        mem_mb=80000,
    shell:
        "{input.script} {input.r1} {input.r2} {output.r1} {output.r2}"


rule trim_adapters:
    output:
        fq=temp("sdata/{stem}.deadapt.fq.gz"),
    input:
        adapters="ref/illumina_adapters.fn",
        fq="sdata/{stem}.fq.gz",
    log:
        "log/{stem}.scythe.log",
    threads: 2
    shell:
        dd(
            """
        scythe -a {input.adapters} {input.fq} >{output.fq} 2>{log}
        ! grep -Fxq 'Blank FASTA header or sequence in adapters file.' {log}
        """
        )


rule quality_trim_reads:
    output:
        r1=temp("sdata/{mgen}.r1.{stem}.qtrim.fq.gz"),
        r2=temp("sdata/{mgen}.r2.{stem}.qtrim.fq.gz"),
        r3=temp("sdata/{mgen}.r3.{stem}.qtrim.fq.gz"),
    input:
        r1="sdata/{mgen}.r1.{stem}.fq.gz",
        r2="sdata/{mgen}.r2.{stem}.fq.gz",
    params:
        qual_type="sanger",
        qual_thresh=20,
    shell:
        dd(
            """
        sickle pe -t {params.qual_type} -q {params.qual_thresh} --gzip-output \
            --pe-file1 {input.r1} --pe-file2 {input.r2} \
            --output-pe1 {output.r1} --output-pe2 {output.r2} \
            --output-single {output.r3}
        """
        )


# NOTE: Input from sdata/ output to data/
rule filter_out_host:
    output:
        r1="data/{mgen}.r1.{stem}.hfilt.fq.gz",
        r2="data/{mgen}.r2.{stem}.hfilt.fq.gz",
    input:
        script="scripts/filter_out_mapping_reads.sh",
        r1="sdata/{mgen}.r1.{stem}.fq.gz",
        r2="sdata/{mgen}.r2.{stem}.fq.gz",
        index=[
            "ref/GRCh38.1.bt2",
            "ref/GRCh38.2.bt2",
            "ref/GRCh38.3.bt2",
            "ref/GRCh38.4.bt2",
            "ref/GRCh38.rev.1.bt2",
            "ref/GRCh38.rev.2.bt2",
        ],
    params:
        index="ref/GRCh38",
    shell:
        dd(
            """
        {input.script} {params.index} {input.r1} {input.r2} {output.r1} {output.r2}
        """
        )


# NOTE: Hub-rule: Comment out this rule to reduce DAG-building time
# once all libraries have been processed
rule alias_cleaned_reads:
    output:
        "data/{stem}.proc.fq.gz",
    input:
        "data/{stem}.dedup.deadapt.qtrim.hfilt.fq.gz",
    shell:
        alias_recipe


localrules:
    alias_cleaned_reads,
