# Snakemake workflow: Identification of potential cancer neoepitopes from DNA-seq data
This workflow detects genomic variants with [Strelka](https://github.com/Illumina/strelka) and and tries to incorporate germline and somatic variants into a sample-specific peptidome using [Microphaser](https://github.com/koesterlab/microphaser). Somatic neopeptides are filtered against germline peptides in terms of similarity and MHC affinity.
Affinity prediction is performed using ([netMHCpan](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?netMHCpan), [netMHC2pan](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?netMHCIIpan)).
The workflow allows easy definition of tumor-normal pairs, where multiple tumor samples (e.g. metastasis) can be grouped with the same normal sample.

**Software requirements**
Unfortunately, the use of netMHCpan and netMHCIIpan is restricted to academic use only. To use those tools in the pipeline, please download them separately following instructions at https://services.healthtech.dtu.dk/software.php.
