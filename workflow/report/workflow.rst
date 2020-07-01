Neoantigen Prediction from Whole Exome/Genome Sequencing (and RNA-Seq)

This workflow predicts tumor-specific neoantigens towards HLA class I and II from NGS data. It incorporates SNVs, MNVs and short InDels.
Adapters were removed with Cutadapt_. Reads were mapped with `BWA MEM`_, PCR and optical duplicates were removed with Picard_.
Variants are called two-fold. Traditional variant calling (somatic and germline) was perfromed using Strelka_. Furthermore, candidate variant discovery was performed with Freebayes_. Statisticall assessment of those candidate variants was conducted with Varlociraptor_.
Somatic variants were annotated using SnpEff_ to predict and report variant effects.
Neopeptide generation was performed using the small-scale phasing approach from Microphaser_.
HLA typing was conducted using Optitype_ and HLA_LA_.
MHC binding affinity scores were predicted using netMHCpan_ and netMHCIIpan_.

Transcript expression levels were analysed using kallisto_ on available RNASeq data.


.. _Varlociraptor: https://varlociraptor.github.io
.. _Cutadapt: https://cutadapt.readthedocs.io
.. _Picard: https://broadinstitute.github.io/picard
.. _Freebayes: https://github.com/ekg/freebayes
.. _BWA mem: http://bio-bwa.sourceforge.net/
.. _SnpEff: http://snpeff.sourceforge.net
.. _MultiQC: http://multiqc.info/
.. _Strelka: https://github.com/Illumina/strelka
.. _netMHCpan: http://www.cbs.dtu.dk/services/NetMHCpan/index.php
.. _netMHCIIpan: http://www.cbs.dtu.dk/services/NetMHCIIpan/index.php
.. _OptiType: https://github.com/FRED-2/OptiType
.. _Microphaser: https://github.com/koesterlab/microphaser
.. _kallisto: https://pachterlab.github.io/kallisto/
.. _HLA_LA: https://github.com/DiltheyLab/HLA-LA
