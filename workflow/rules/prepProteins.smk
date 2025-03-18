"""

Copyright Richard Stöckl 2025.
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE or copy at 
https://www.boost.org/LICENSE_1_0.txt)

"""

rule predict_cds:
    input:
        fasta=lambda wildcards: path[wildcards.id],
    output:
        faa=os.path.join(RESULTPATH, "{id}", "cds", "{id}.faa"), 
        gff=os.path.join(RESULTPATH, "{id}", "cds", "{id}.gff"),
        gbk=os.path.join(RESULTPATH, "{id}", "cds", "{id}.gbk"),
        fna=os.path.join(RESULTPATH, "{id}", "cds", "{id}.fna"),
        tsv=os.path.join(RESULTPATH, "{id}", "cds", "{id}.tsv"),
        rrna=os.path.join(RESULTPATH, "{id}", "cds", "{id}_rrna.tsv"),
    log:
        os.path.join(LOGPATH, "{id}", "logs", "{id}_predict_cds.log"),
    message:
        "Predicting CDS for {wildcards.id}",
    threads: 1
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
        scriptpath= os.path.join(workflow.basedir, "scripts","cds.py")
    conda:
        os.path.join(workflow.basedir, "envs","cds.yaml"),
    shell:
        """
        python {params.scriptpath} {params.sample_id} {input.fasta} {output.faa} {output.gff} {output.gbk} {output.fna} {output.tsv} \
        --translation_table {params.translation_table} \
        --write_faa --write_gff --write_gbk --write_fna --write_tsv \
        --run_rrna --rrna_tsv {output.rrna} > {log} 2>&1
        """