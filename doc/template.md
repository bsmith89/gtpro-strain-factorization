# Description:

A Snakemake project template.

# Converting template into a project:

1. Clone this repository.
2. `git remote rename origin template`
3. `git checkout -b master`
4. `ln -s doc/NOTES.md README.md` (optional)
5. `smake -j1 initialize_project`
6. Create conda environment and install necessary packages.

# TODO:

-   [ ] Make the default conda environment truly generic.
-   [ ] Write a base project initialization recipe in `Snakefile`.
-   [ ] Alternatively: switch to [cookiecutter](https://cookiecutter.readthedocs.io/en/latest/).
