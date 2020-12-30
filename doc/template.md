# Description:

A Snakemake project template.

# Converting template into a project:

1. Clone this repository.
2. `git remote rename origin template`
3. `git checkout -b master`
4. `ln -fs doc/notes.md README.md` (optional)
6. Create conda environment and install Snakemake (at a minimum).
7. Edit `env_local` to activate this environment.
7. `source env`
5. `smake -j1 initialize_project`
10. Install any other required packages/software, e.g. `conda/default.yaml`

# TODO:

-   [ ] Make the default conda environment truly generic.
-   [ ] Write a base project initialization recipe in `Snakefile`.
-   [ ] Alternatively: switch to [cookiecutter](https://cookiecutter.readthedocs.io/en/latest/).
