rule initialize_project:
    input:
        local_config_files=[ancient(p) for p in ["config_local.yaml", "env_local", "snake/local.smk"]],
    shell:
        dd(
            """
        git config --local filter.dropoutput_ipynb.clean scripts/ipynb_output_filter.py
        git config --local filter.dropoutput_ipynb.smudge cat
        git config --local diff.daff-csv.command "daff.py diff --git"
        git config --local merge.daff-csv.name "daff.py tabular merge"
        git config --local merge.daff-csv.driver "daff.py merge --output %A %O %A %B"
        echo 'Add local configuration to {input.local_config_files}'
        echo 'Please activate your environment and then run `pip install -r requirements.txt` or analogous.'
        """
        )

rule build_empty_local_config_yaml:
    output:"config_local.yaml"
    shell: "echo 'DUMMY_: 1' > {output}"

rule build_empty_local_snakefile:
    output:"snake/local.smk"
    shell: "touch {output}"

rule build_empty_local_env:
    output:"env_local"
    shell: "touch {output}"
