Neoantigen Prediction from Whole Exome Sequencing (and RNA-Seq)

Neoantigen candidate prediction from WES data on six melanoma samples with one healthy normal sample.
RNA information was available for Samples Mel-103a and Mel-103b. For those samples, neoantigens can be filtered by TPM (Transcript per Million) value.


Reads were mapped onto {{ snakemake.config["reference"]["build"] }} with `BWA mem`_, and both optical and PCR duplicates were removed with Picard_.
The Strelka_ variant caller was used to call somatic and germline variants for all melanoma samples. Germline variants were also called for the normal sample.
Somatic variants were annotated using SnpEff_ to predict and report variant effects.

The variant calls were phased on every melanoma sample using Microphaser_. Resulting neopeptides with a lenght of nine aminoacids were filtered for self-similarity against all normal peptides.

HLA-types of all samples were predicted using OptiType_.

The MHC binding affinity for all neoantigen candidates and their unmutated counterparts was predicted using netMHCpan_. The developers define cut-off values for binding peptides as percentage-rank < 0.5% for strong binders and percentage-rank < 2% for weak binders.

For available RNA-Seq data, quantification was performed using kallisto_. The TPM values for every transcript were added to the output table.

.. _BWA mem: http://bio-bwa.sourceforge.net/
.. _SnpEff: http://snpeff.sourceforge.net
.. _MultiQC: http://multiqc.info/
.. _Strelka: https://github.com/Illumina/strelka
.. _Picard: https://broadinstitute.github.io/picard/
.. _FastQC: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
.. _netMHCpan: http://www.cbs.dtu.dk/services/NetMHCpan/index.php
.. _OptiType: https://github.com/FRED-2/OptiType
.. _Microphaser: https://github.com/koesterlab/microphaser
.. _kallisto: https://pachterlab.github.io/kallisto/

