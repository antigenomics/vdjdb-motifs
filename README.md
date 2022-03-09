## VDJdb motif inference scripts

One should clone [vdjdb-db repo](https://github.com/antigenomics/vdjdb-db) and [this repo](https://github.com/antigenomics/vdjdb-motifs) to the same directory, say ``~/vcs/`` and then navigate to ``~/vcs/vdjdb-db`` and run ``bash release.sh`` which will build the VDJdb database, navigate to motifs repository and run motif inference.

This folder lacks ``pools/`` directory with control TCR sequences that is ~500Mb compressed, it will be automatically downloaded from this [Zenodo record](https://zenodo.org/record/6339774) when running ``compute_vdjdb_motifs.Rmd`` or ``release.sh``.