Neoantigens and corresponding normal peptides as phased and determined by
microphaser, with elution ligand / binding affinity predictions by netMHC
(netMHCpan and netMHCIIpan) to the HLA alleles determined by HLA-LA.

===================
Column descriptions
===================

* **id**: Peptide ID assigned by microphaser.
* **pep_seq**: Sequence of the full peptide that was given to netMHC(II)pan. Amino acids that
    are different between the normal and the tumor sample are highlighted in lower case.
* **pos_in_id_seq**: Position in pep_seq of the peptide that was used for prediction.
    It seems like netMHCpan positions start at 0 and netMHCIIpan positions at 1.
* **alias**: Indicator of the type of sample (normal vs. some tumor sample).
* **num_binders**: Total number of peptide-HLA allele pairs from pep_seq that are considered binders,
    either weak or strong. Cutoffs are applied to el_rank and are:
    * netMHCpan 4.1: <0.5% (strong binder), <2.0% (weak binder)
    * netMHCIIpan 4.1: <1.0% (strong binder), <5.0% (weak binder)
* **freq**: Allelic frequency of the peptide as predicted by microphaser. For a
    credible allele frequency interval, see column freq_credible_interval.
* **depth**: Read depth at the peptide position.
* **num_var_sites**: Number of variant sites on the peptides haplotype in that sample (alias).
* **num_var_in_pep**: Number of variant sites within the peptide sequence.
* **top_el_rank_allele**: HLA allele with the best eluted ligand prediction score percentile rank.
* **top_el_rank_bind_core**: Binding core of the peptide for the HLA allele with the best
    elution ligand score percentile rank.
* **top_el_rank_el_rank**: Percentile rank of the elution ligand score. The rank of the predicted binding
    score when compared to a set of random natural peptides. This measure is not
    affected by inherent bias of certain molecules towards higher or lower mean
    predicted affinities. It is the recommended value for determining likely
    binders / neoantigens of interest. Cutoffs recommended by netMHC(II)pan authors
    are:
    * netMHCpan 4.1: <0.5% (strong binder), <2.0% (weak binder)
    * netMHCIIpan 4.1: <1.0% (strong binder), <5.0% (weak binder)
* **top_el_rank_el_score**: The raw eluted ligand prediction score
* **ave_el_score**: Average across the eluted ligand prediction scores of all alleles for this
    peptide in the particular sample (alias).
* **top_ba_rank_allele**: HLA allele with the best binding affinity prediction score percentile rank.
* **top_ba_rank_bind_core**: Binding core of the peptide for the HLA allele with the best
    binding affinity score percentile rank.
* **top_ba_rank_ba_rank**: Percentile rank of the predicted binding affinity compared to a set of 100.000
    random natural peptides. This measure is not affected by inherent bias of certain
    molecules towards higher or lower mean predicted affinities.
* **top_ba_rank_ba_score**: Predicted binding affinity in log-scale.
* **aa_changes**: List of aa changes. For normal samples, only germline changes are listed. For
    tumor samples, only somatic tumor changes are listed, even though the germline
    changes also affect the tumor peptide.
* **genomic_pos**: Genomic position of the nucleotide change.
* **nt_seq**: Nucleotide sequence underlying the peptide. Nucleotide changes underlying the
    amino acid changes are highlighted as lower case letters.    
* **gene_name**: Common gene name / gene symbol of the peptide's gene of origin.
* **gene_id**: Ensembl gene id of the peptide's gene of origin.
* **transcript**: List of Ensembl transcript ids in which the peptide occurs.
* **chrom**: Chromosome on which the peptide's gene of origin is located.    
* **offset**: Chromosomal position of the peptide's gene of origin.    
* **frame**: Open reading frame that the peptide originates from. 0 indicates the regular
    reading frame, non-zero values indicate frame shifts.
* **strand**: Strand of the gene / transcript.
* **freq_credible_interval**: Credible interval for freq, the allelic frequency of the peptide as predicted
    by microphaser.
