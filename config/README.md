# General settings
To configure this workflow, modify ``config/config.yaml`` according to your needs, following the explanations provided in the file.

# Sample sheet
Add samples to `config/samples.tsv`. For each sample, the columns `sample`, `type` and `platform` have to be defined. The columns `matched_normal` and `purity` are only needed for tumor samples.
* The `sample` column contains the name of a specific biological sample.
* The `type` (normal, tumor) will be used to differentiate between multiple conditions of the same biological entity (e.g. a sinlge patient).
* `matched_normal` should contain the sample name of the corresponding normal sample to a tumor sample.
* `purity` should contain an estimated purity value for the tumor sample in the range of `]0.0, 1.0]`.

# Unit sheet
For each sample, add one or more sequencing units (runs, lanes or replicates) to the unit sheet `config/units.tsv`. By activating or deactivating `mergeReads` in the `config/config.yaml`, you can decide wether to merge replicates or run them individually.
* Each unit has a column `unit` which can be e.g. a running number, or an actual run, lane or replicate id.
* Each unit has a column `sample` which refers to the biological sample of the unit.
* For each unit, define either one (column `fq1`) or two (columns `fq1`, `fq2`) FASTQ files (these can point to anywhere in your system). Alternatively, you can define an SRA (sequence read archive) accession (starting with e.g. ERR or SRR) by using a column `sra`. In the latter case, the pipeline will automatically download the corresponding paired end reads from SRA. If both local files and SRA accession are available, the local files will be preferred.
* Define adapters in the adapters column, by putting cutadapt arguments in quotation marks (e.g. "-a ACGCGATCG -A GCTAGCGTACT").

Missing values can be specified by empty columns or by writing `NA`. Lines can be commented out with `#`.

# netMHC configuration
For the use of netMHCpan and netMHCIIpan, you need to install the software manually by downloading them separately following instructions at https://services.healthtech.dtu.dk/software.php.
Once you have downloaded and extracted the provided folders, follow the instructions to customise the `netMHCpan` or `netMHCIIpan` file. To call the correct installation of netMHC, please give the path to the netMHC folder in corresponding field of the `config.yaml`, e.g. `path/to/netMHCpan-4.1`.
