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

        ## pre-installed via conda
        export NEOFOX_MAKEBLASTDB=makeblastdb
        echo 'NEOFOX_MAKEBLASTDB=makeblastdb' > {output.config}
        export NEOFOX_RSCRIPT=Rscript
        echo 'NEOFOX_RSCRIPT=Rscript' >> {output.config}

        ## pre-installed into conda environment via post-deploy script
        export NEOFOX_NETMHCPAN=netMHCpan
        echo 'NEOFOX_NETMHCPAN=netMHCpan' >> {output.config}
        export NEOFOX_NETMHC2PAN=netMHCIIpan
        echo 'NEOFOX_NETMHC2PAN=netMHCIIpan' >> {output.config}

        ## specification of hla_allele link via config.yaml
        export NEOFOX_HLA_DATABASE={params.hla_alleles}

        neofox-configure --reference-folder {output.references}
        echo 'NEOFOX_REFERENCE_FOLDER={output.references}' >> {output.config}

        # further environment variables needed for the config file

        ## pre-installed via conda
        echo 'NEOFOX_BLASTP=blastp' >> {output.config}
        
        ## pre-installed into conda environment via post-deploy script
        echo 'NEOFOX_MIXMHCPRED=MixMHCpred' >> {output.config}
        echo 'NEOFOX_MIXMHC2PRED=MixMHC2pred_unix' >> {output.config}
        echo 'NEOFOX_PRIME=PRIME' >> {output.config}
        """


rule neo_fox:
    input:
        config="resources/neo_fox/neo_fox_config.txt",
        references=directory("resources/neo_fox/references/"),
        candidates="results/microphaser/info/filtered/{group}.{tumor_alias}.merged_tumor_normal.neo_fox.tsv",
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
