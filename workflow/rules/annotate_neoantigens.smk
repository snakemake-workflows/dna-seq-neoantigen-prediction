rule prepare_neo_fox_config_and_resources:
    output:
        config="resources/neo_fox/neo_fox_config.txt",
        references=directory("resources/neo_fox/references/"),
    conda:
        "../envs/neo_fox_deps.yaml"
    params:
        hla_alleles=config["params"]["neo_fox"]["hla_alleles"],
    shell:
        """
        # environment variables necessary for neofox-configure
        # NOTE: we have to provide all binaries with hard-coded
        # paths, because NeoFox checks that the file at the path
        # exists (and not simply that a binary exists):
        # https://github.com/TRON-Bioinformatics/neofox/blob/629443b637fc41b1ab81f4f770e7a8a1c976d3f2/neofox/references/references.py#L90
        CONDA_BIN=$CONDA_PREFIX/bin

        ## pre-installed via conda
        export NEOFOX_MAKEBLASTDB=$CONDA_BIN/makeblastdb
        echo 'NEOFOX_MAKEBLASTDB=$CONDA_BIN/makeblastdb' > {output.config}
        export NEOFOX_RSCRIPT=$CONDA_BIN/Rscript
        echo 'NEOFOX_RSCRIPT=$CONDA_BIN/Rscript' >> {output.config}

        ## pre-installed into conda environment via post-deploy script
        export NEOFOX_NETMHCPAN=$CONDA_BIN/netMHCpan
        echo 'NEOFOX_NETMHCPAN=$CONDA_BIN/netMHCpan' >> {output.config}
        export NEOFOX_NETMHC2PAN=$CONDA_BIN/netMHCIIpan
        echo 'NEOFOX_NETMHC2PAN=$CONDA_BIN/netMHCIIpan' >> {output.config}

        ## specification of hla_allele link via config.yaml
        export NEOFOX_HLA_DATABASE={params.hla_alleles}

        neofox-configure --reference-folder {output.references}
        echo 'NEOFOX_REFERENCE_FOLDER={output.references}' >> {output.config}

        # further environment variables needed for the config file

        ## pre-installed via conda
        echo 'NEOFOX_BLASTP=$CONDA_BIN/blastp' >> {output.config}
        
        ## pre-installed into conda environment via post-deploy script
        echo 'NEOFOX_MIXMHCPRED=$CONDA_BIN/MixMHCpred' >> {output.config}
        echo 'NEOFOX_MIXMHC2PRED=$CONDA_BIN/MixMHC2pred_unix' >> {output.config}
        echo 'NEOFOX_PRIME=$CONDA_BIN/PRIME' >> {output.config}
        """


rule adjust_microphaser_output_for_neo_fox:
    input:
        candidates="results/microphaser/info/filtered/{group}.{tumor_alias}.merged_tumor_normal.pep_len_{peptide_length}.tsv",
    output:
        candidates="results/neo_fox/candidates/{group}.{tumor_alias}.merged_tumor_normal.pep_len_{peptide_length}.tsv",
    threads: 1
    conda:
        "../envs/polars.yaml"
    script:
        "../scripts/adjust_microphaser_output_for_neo_fox.py"


rule neo_fox:
    input:
        config="resources/neo_fox/neo_fox_config.txt",
        references=directory("resources/neo_fox/references/"),
        candidates=expand(
            "results/neo_fox/candidates/{{group}}.{{tumor_alias}}.merged_tumor_normal.pep_len_{peptide_length}.tsv",
            peptide_length=config["params"]["neo_fox"]["peptide_len"],
        ),
        patient_annotation="results/neo_fox/patient_data/{group}.{tumor_alias}.hla_alleles.tumor_type.tsv",
    output:
        tsv="results/neo_fox/annotated/{group}.{tumor_alias}.annotated_neoantigens.tsv",
        json="results/neo_fox/annotated/{group}.{tumor_alias}.annotated_neoantigens.json",
        meta_json="results/neo_fox/annotated/{group}.{tumor_alias}.meta_annotations.json",
    threads: 8
    conda:
        "../envs/neo_fox_deps.yaml"
    params:
        folder=lambda wc, output: path.dirname(output.annotated),
        prefix=lambda wc, output: path.plitext(path.basename(output.annotated))[0],
        organism="human" if config["ref"]["species"]=="homo_sapiens" else "mouse" if config["ref"]["species"]=="mus_musculus" else "unsupported",
    shell:
        "(neofox "
        "  --num_cpus {threads} "
        "  --config {input.config} "
        "  --candidate-file {input.candidates} "
        "  --patient-data {input.patient_annotation} "
        "  --with-table "
        "  --with-json "
        "  --organism {params.organism} "
        "  --output-folder {params.folder} "
        "  --output-prefix {params.prefix} ; "
        " mv {params_folder}/{params.prefix}_neoantigen_candidates_annotated.tsv {output.tsv}; "
        " mv {params_folder}/{params.prefix}_neoantigen_candidates_annotated.json {output.json}; "
        " mv {params_folder}/{params.prefix}_neoantigen_features.json {output.meta_json}; "
        ") 2> {log} "
