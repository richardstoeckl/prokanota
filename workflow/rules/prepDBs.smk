"""

Copyright Richard Stöckl 2025.
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE or copy at 
https://www.boost.org/LICENSE_1_0.txt)

"""

rule downloadCDDdb:
    output:
        db=protected(os.path.join(DBPATH, "CDD","Cdd_LE.tar.gz")),
        mapping=protected(os.path.join(DBPATH, "CDD", "cddid_all.tbl.gz")),
        version=os.path.join(DBPATH, "CDD", "cdd.info"),
        done=os.path.join(DBPATH, "CDD", "download.done"),
    log:
        os.path.join(LOGPATH, "common", "downloadCDDdb.log"),
    message:
        "Downloading CDD database to {output.db}",
    conda:
        os.path.join(workflow.basedir, "envs","environment.yaml"),
    shell:
        """
        curl -L -o {output.db} https://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/Cdd_LE.tar.gz
        curl -L -o {output.mapping} https://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddid_all.tbl.gz
        curl -L -o {output.version} https://ftp.ncbi.nih.gov/pub/mmdb/cdd/cdd.info
        touch {output.done}
        """

rule prepareCDDdb:
    input:
        db=rules.downloadCDDdb.output.db,
        mapping=rules.downloadCDDdb.output.mapping,
        version=rules.downloadCDDdb.output.version,
        done=rules.downloadCDDdb.output.done,
    output:
        db=directory(os.path.join(DBPATH, "CDD","db")),
        mapping=os.path.join(DBPATH, "CDD", "cddid_all.tbl"),
        version=os.path.join(DBPATH, "CDD", "cdd_info.txt"),
    message:
        "Preparing CDD database",
    shell:
        """
        mkdir -p {output.db}
        tar -xzf {input.db} -C {output.db}
        gunzip -c {input.mapping} > {output.mapping}
        cp {input.version} {output.version}
        """

rule downloadCOGdb:
    output:
        db=protected(os.path.join(DBPATH, "COG","Cog_LE.tar.gz")),
        mapping=protected(os.path.join(DBPATH, "COG", "cddid_all.tbl.gz")),
        version=os.path.join(DBPATH, "COG", "cdd.info"),
        done=os.path.join(DBPATH, "COG", "download.done"),
    log:
        os.path.join(LOGPATH, "common", "downloadCOGdb.log"),
    message:
        "Downloading COG database to {output.db}",
    conda:
        os.path.join(workflow.basedir, "envs","environment.yaml"),
    shell:
        """
        curl -L -o {output.db} https://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/Cog_LE.tar.gz
        curl -L -o {output.mapping} https://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddid_all.tbl.gz
        curl -L -o {output.version} https://ftp.ncbi.nih.gov/pub/mmdb/cdd/cdd.info
        touch {output.done}
        """

rule prepareCOGdb:
    input:
        db=rules.downloadCOGdb.output.db,
        mapping=rules.downloadCOGdb.output.mapping,
        version=rules.downloadCOGdb.output.version,
        done=rules.downloadCOGdb.output.done,
    output:
        db=directory(os.path.join(DBPATH, "COG","db")),
        mapping=os.path.join(DBPATH, "COG", "cddid_all.tbl"),
        version=os.path.join(DBPATH, "COG", "cdd_info.txt"),
    message:
        "Preparing COG database",
    shell:
        """
        mkdir -p {output.db}
        tar -xzf {input.db} -C {output.db}
        gunzip -c {input.mapping} > {output.mapping}
        cp {input.version} {output.version}
        """

rule downloadPGAPdb:
    output:
        db=protected(os.path.join(DBPATH, "PGAP","hmm_PGAP.LIB")),
        mapping=protected(os.path.join(DBPATH, "PGAP", "hmm_PGAP.tsv")),
        version=os.path.join(DBPATH, "PGAP", "RELEASE_NOTES.txt"),
    log:
        os.path.join(LOGPATH, "common", "downloadPGAPdb.log"),
    message:
        "Downloading PGAP database to {output.db}",
    conda:
        os.path.join(workflow.basedir, "envs","environment.yaml"),
    shell:
        """
        curl -L -o {output.db} https://ftp.ncbi.nlm.nih.gov/hmm/17.0/hmm_PGAP.LIB
        curl -L -o {output.mapping} https://ftp.ncbi.nlm.nih.gov/hmm/17.0/hmm_PGAP.tsv
        curl -L -o {output.version} https://ftp.ncbi.nlm.nih.gov/hmm/17.0/RELEASE_NOTES.txt
        """

rule preparePGAPdb:
    input:
        db=rules.downloadPGAPdb.output.db,
        mapping=rules.downloadPGAPdb.output.mapping,
        version=rules.downloadPGAPdb.output.version,
    output:
        mapping=os.path.join(DBPATH, "PGAP", "hmm_PGAP_clean.tsv"),
    message:
        "Preparing PGAP database",
    conda:
        os.path.join(workflow.basedir, "envs","environment.yaml"),
    shell:
        """
        sed 's/ /_/g' {input.mapping} | gawk 'BEGIN {{ FS = OFS = "\t" }} {{ for(i=1; i<=NF; i++) if($i == "-") $i="*"}}; 1' | sed '1d' > {output.mapping}
        """

rule downloadARCOGdb:
    output:
        db=os.path.join(DBPATH, "arCOG","zip.aliar14.tgz"),
        mapping=os.path.join(DBPATH, "arCOG", "arCOGdef.tab"),
        version=os.path.join(DBPATH, "arCOG", "arCOG.info"),
        done=os.path.join(DBPATH, "arCOG", "download.done"),
    message:
        "Downloading arCOG database to {output.db}",
    log:
        os.path.join(LOGPATH, "common", "downloadARCOGdb.log"),
    conda:
        os.path.join(workflow.basedir, "envs","environment.yaml"),
    shell:
        """
        curl -L -o {output.db} https://ftp.ncbi.nih.gov/pub/wolf/COGs/arCOG/zip.aliar14.tgz 2>{log}
        curl -L -o {output.mapping} https://ftp.ncbi.nlm.nih.gov/pub/wolf/COGs/arCOG/tmp.ar18/arCOGdef.tab 2>{log}
        echo 'arCOG DB downloaded via ftp.ncbi.nih.gov/pub/wolf/COGs/arCOG/zip.aliar14.tgz. Mapping file downloaded via ftp.ncbi.nlm.nih.gov/pub/wolf/COGs/arCOG/tmp.ar18/arCOGdef.tab' > {output.version} 2> {log}
        touch {output.done}
        """

rule prepareARCOGdb:
    input:
        db=rules.downloadARCOGdb.output.db,
        mapping=rules.downloadARCOGdb.output.mapping,
        done=rules.downloadARCOGdb.output.done,
    output:
        # dbdir=directory(os.path.join(DBPATH, "arCOG", "raw")),
        # srfiles=os.path.join(DBPATH, "arCOG", "raw","{i}.sr"),
        # hmm=protected(os.path.join(DBPATH, "arCOG", "arCOGs_combined.hmm")),
        h3f=protected(os.path.join(DBPATH, "arCOG", "arCOGs_combined.h3f")),
        h3i=protected(os.path.join(DBPATH, "arCOG", "arCOGs_combined.h3i")),
        h3m=protected(os.path.join(DBPATH, "arCOG", "arCOGs_combined.h3m")),
        h3p=protected(os.path.join(DBPATH, "arCOG", "arCOGs_combined.h3p")),
        # mapping=os.path.join(DBPATH, "arCOG", "arCOGdef_clean.tab"),
    conda:
        os.path.join(workflow.basedir, "envs","pyhmmer.yaml"),
    message:
        "Preparing arCOG database",
    threads:
        max(1, int(workflow.cores * 1))  # Use safe threads calculation
    log:
        os.path.join(LOGPATH, "common", "prepareARCOGdb.log"),
    params:
        script=os.path.join(workflow.basedir, "scripts", "prepareArCOGdb.py"),
        outdir=os.path.join(DBPATH, "arCOG"),
    shadow: "shallow"
    shell:
        """
        python {params.script} --tar {input.db} --output-dir {params.outdir} --workers {threads} > {log} 2>&1
        """

rule prepareARCOGdbMapping:
    input:
        mapping=rules.downloadARCOGdb.output.mapping,
    output:
        mapping=os.path.join(DBPATH, "arCOG", "arCOGdef_clean.tab"),
    conda:
        os.path.join(workflow.basedir, "envs","environment.yaml"),
    message:
        "Preparing arCOG mapping file",
    shell:
        """
        csvtk tab2csv {input.mapping} | csvtk csv2tab | gawk 'BEGIN {{ FS = OFS = "\t" }} {{ for(i=1; i<=NF; i++) if($i == "-") $i="*"}}; 1' > {output.mapping}
        """
