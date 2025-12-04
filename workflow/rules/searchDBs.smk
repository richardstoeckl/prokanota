"""

Copyright Richard Stöckl 2025.
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE or copy at 
https://www.boost.org/LICENSE_1_0.txt)

"""

rule searchCDDs:
    input:
        allFaa=rules.predict_features.output.faa, # rules.concatinateFaa.output.allFaa,
        db=rules.prepareCDDdb.output.db,
        version=rules.prepareCDDdb.output.version,
    output:
        tab=os.path.join(INTERIMPATH, "{id}", "cdd", "results_all.tsv"),
        dbversion=os.path.join(RESULTPATH, "{id}", "versions", "CDD_DB.txt"),
        toolversion=os.path.join(RESULTPATH, "{id}", "versions", "CDD_rpsblast.txt"),
    message:
        "Searching CDDs for {wildcards.id}",
    threads:
        max(1, int(workflow.cores * 0.5))  # Use safe threads calculation
    log: 
        os.path.join(LOGPATH, "{id}", "logs", "searchCDD.log"),
    benchmark: 
        os.path.join(LOGPATH, "{id}", "benchmarks", "searchCDD.tsv"),
    params:
        dbPath=os.path.join(rules.prepareCDDdb.output.db, "Cdd"),
    conda:
        os.path.join(workflow.basedir, "envs","environment.yaml"),
    shell:
        """
        cp {input.version} {output.dbversion}
        rpsblast -version > {output.toolversion}
        rpsblast -query {input.allFaa} -db {params.dbPath} -out {output.tab} -seg no -comp_based_stats 1 -outfmt 6 -evalue 0.001 -num_threads {threads} > {log} 2>&1
        """
 
rule searchCOGs:
    input:
        allFaa=rules.predict_features.output.faa, # rules.concatinateFaa.output.allFaa,
        db=rules.prepareCOGdb.output.db,
        version=rules.prepareCOGdb.output.version,
    output:
        tab=os.path.join(INTERIMPATH, "{id}", "cogs", "results_all.tsv"),
        dbversion=os.path.join(RESULTPATH, "{id}", "versions", "COG_DB.txt"),
        toolversion=os.path.join(RESULTPATH, "{id}", "versions", "COG_rpsblast.txt"),
    message:
        "Searching COGs for {wildcards.id}",
    threads:
        max(1, int(workflow.cores * 0.5))  # Use safe threads calculation
    log: 
        os.path.join(LOGPATH, "{id}", "logs", "searchCOGs.log")
    benchmark: 
        os.path.join(LOGPATH, "{id}", "benchmarks", "searchCOGs.tsv")
    params:
        dbPath=os.path.join(rules.prepareCOGdb.output.db, "Cog"),
    conda:
        os.path.join(workflow.basedir, "envs","environment.yaml"),
    shell:
        """
        cp {input.version} {output.dbversion}
        rpsblast -version > {output.toolversion}
        rpsblast -query {input.allFaa} -db {params.dbPath} -out {output.tab} -seg no -comp_based_stats 1 -outfmt 6 -evalue 0.001 -num_threads {threads} > {log} 2>&1
        """

rule searchArCOGs:
    input:
        allFaa=rules.predict_features.output.faa, #rules.concatinateFaa.output.allFaa,
        db=rules.prepareARCOGdb.output.h3m,
        version=rules.downloadARCOGdb.output.version,
    output:
        tblout=os.path.join(INTERIMPATH, "{id}", "arcogs", "sequence_results.txt"),
        # allresults=os.path.join(RESULTPATH, "{id}", "arcogs", "results_all.txt"),
        # domtblout=os.path.join(RESULTPATH, "{id}", "arcogs", "domain_results.txt"),
        dbversion=os.path.join(RESULTPATH, "{id}", "versions", "ArCOG_DB.txt"),
        toolversion=os.path.join(RESULTPATH, "{id}", "versions", "ArCOG_hmmsearch.txt"),
    message:
        "Searching ArCOGs for {wildcards.id}",
    threads:
        max(1, int(workflow.cores * 0.5))  # Use safe threads calculation
    log: 
        os.path.join(LOGPATH, "{id}", "logs", "searchArCOGs.log")
    benchmark: 
        os.path.join(LOGPATH, "{id}", "benchmarks", "searchArCOGs.tsv")
    params:
        script=os.path.join(workflow.basedir, "scripts", "pyhmmer_search.py"),
    conda:
        os.path.join(workflow.basedir, "envs","pyhmmer.yaml"),
    shell:
        """
        python {params.script} \
        --db {input.db} \
        --faa {input.allFaa} \
        --tblout {output.tblout} \
        --toolversion {output.toolversion} \
        --threads {threads} > {log} 2>&1 && cp {input.version} {output.dbversion}
        """

rule searchPGAP:
    input:
        allFaa=rules.predict_features.output.faa, #rules.concatinateFaa.output.allFaa,
        db=rules.downloadPGAPdb.output.db,
        version=rules.downloadPGAPdb.output.version,
    output:
        tblout=os.path.join(INTERIMPATH, "{id}", "pgap", "sequence_results.txt"),
        # allresults=os.path.join(RESULTPATH, "{id}", "pgap", "results_all.txt"),
        # domtblout=os.path.join(RESULTPATH, "{id}", "pgap", "domain_results.txt"),
        dbversion=os.path.join(RESULTPATH, "{id}", "versions", "PGAP_DB.txt"),
        toolversion=os.path.join(RESULTPATH, "{id}", "versions", "PGAP_hmmsearch.txt"),
    message:
        "Searching PGAP for {wildcards.id}",
    threads:
        max(1, int(workflow.cores * 0.5))  # Use safe threads calculation
    log: 
        os.path.join(LOGPATH, "{id}", "logs", "searchPGAP.log")
    benchmark: 
        os.path.join(LOGPATH, "{id}", "benchmarks", "searchPGAP.tsv")
    params:
        script=os.path.join(workflow.basedir, "scripts", "pyhmmer_search.py"),
    conda:
        os.path.join(workflow.basedir, "envs","pyhmmer.yaml"),
    shell:
        """
        python {params.script} \
        --db {input.db} \
        --faa {input.allFaa} \
        --tblout {output.tblout} \
        --toolversion {output.toolversion} \
        --threads {threads} > {log} 2>&1 && cp {input.version} {output.dbversion}
        """
