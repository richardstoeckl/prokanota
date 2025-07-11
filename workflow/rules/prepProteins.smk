"""

Copyright Richard Stöckl 2025.
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE or copy at 
https://www.boost.org/LICENSE_1_0.txt)

"""

rule predict_features:
    input:
        fasta=lambda wildcards: path[wildcards.id],
    output:
        faa=os.path.join(RESULTPATH, "{id}", "features", "{id}.faa"), 
        gff=os.path.join(RESULTPATH, "{id}", "features", "{id}.gff"),
        gbk=os.path.join(RESULTPATH, "{id}", "features", "{id}.gbk"),
        fna=os.path.join(RESULTPATH, "{id}", "features", "{id}.fna"),
        tsv=os.path.join(RESULTPATH, "{id}", "features", "{id}.tsv"),
        rna=os.path.join(RESULTPATH, "{id}", "features", "{id}_rna.tsv"),
        genome=os.path.join(RESULTPATH, "{id}", "features", "{id}.fasta"),
    log:
        os.path.join(LOGPATH, "{id}", "logs", "{id}_predict_features.log"),
    message:
        "Predicting Features for {wildcards.id}",
    threads: 4
    params:
        sample_id="{id}",
        meta=False,
        closed=False,
        translation_table=11,
        write_faa=True,
        write_gff=True,
        write_gbk=True,
        write_fna=True,
        write_tsv=True,
        scriptpath= os.path.join(workflow.basedir, "scripts","features.py")
    conda:
        os.path.join(workflow.basedir, "envs","features.yaml"),
    shell:
        """
        python {params.scriptpath} \
        {params.sample_id} {input.fasta} \
        --faa_path {output.faa} \
        --gff_path {output.gff} \
        --gbk_path {output.gbk} \
        --fna_path {output.fna} \
        --tsv_path {output.tsv} \
        --genome_path {output.genome} \
        --run_rrna --run_trna \
        --rna_tsv_path {output.rna} \
        --translation_table {params.translation_table} \
        --threads {threads} > {log} 2>&1
        """