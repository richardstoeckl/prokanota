"""

Copyright Richard Stöckl 2025.
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE or copy at 
https://www.boost.org/LICENSE_1_0.txt)

"""

# Names of the Main Result Dataframe:
# A - Starting
# B - CDD
# C - NCBI COGs
# D - arCOGs
# E - PGAP

if config["searches"]["cdd"]:
    rule parseCDDResults:
        input:
            tblout=rules.searchCDDs.output.tab,
            prevDF=rules.predict_cds.output.tsv, # rules.createStartingAnnotDF.output.mainProkka,
            mapping=rules.prepareCDDdb.output.mapping,
        output:
            finalDF=os.path.join(INTERIMPATH, "{id}", "annotation", "mainDF_AB.tsv"),
            ecutoff=temp(os.path.join(INTERIMPATH, "{id}", "cogs", "sequence_results_red_e_cutoff.tsv")),
            allResults=temp(os.path.join(INTERIMPATH, "{id}", "cogs", "All_CDD.tsv")),
            tmp1=temp(os.path.join(INTERIMPATH, "{id}", "cogs", "temp1.tsv")),
            finalDFtmp=temp(os.path.join(INTERIMPATH, "{id}", "cogs", "cogs_final_tmp.tsv")),
            finalOut=temp(os.path.join(INTERIMPATH, "{id}", "cogs", "cogs_final.tsv")),
        message:
            "Parsing CDD results for {wildcards.id}",
        log:
            os.path.join(LOGPATH, "{id}", "logs", "parseCOGResults.log")
        benchmark: 
            os.path.join(LOGPATH, "{id}", "benchmarks", "parseCOGResults.tsv")
        threads: 1
        priority: 10
        params:
            # mappingFile=config["mappingPaths"]["cog_mapping"],
            mainDFdir=lambda w, output:os.path.split(os.path.splitext(output[0])[0])[0],
        conda:
            os.path.join(workflow.basedir, "envs","environment.yaml"),
        shell:
            #format the full table and only select sequences above a certain e-value
            #get best hit based on bit score, and then e-value
            #merge with mapped names
            #add in an extra column that lists whether hits have a high confidence score
            #add in header
            #combine with previous dataframe
            #control lines
            """
            sed 's/ \\+ /\t/g' {input.tblout} | sed '/^#/d'| sed 's/ /\t/g' | sed 's/CDD://g' | gawk -F'\t' -v OFS='\t' '{{print $1, $2, $12, $11}}' | gawk -F'\t' -v OFS='\t' '($4 + 0) <= 1e-3'  > {output.ecutoff}
            sort -t$'\t' -k3,3gr -k4,4g {output.ecutoff} | sort -t$'\t' --stable -u -k1,1  | sort -t$'\t' -k3,3gr -k4,4g > {output.allResults}        
            LC_ALL=C join -a1 -1 2 -2 1 -e'*' -t $'\t' -o1.1,2.2,2.3,2.4,1.4 <(LC_ALL=C sort -k2 {output.allResults}) <(LC_ALL=C sort -k1 {input.mapping}) | LC_ALL=C  sort > {output.tmp1}        
            echo -e "gene_id\tNCBI_CDD\tNCBI_CDD_GeneID\tNCBI_CDD_Description\tNCBI_CDD_evalue" | cat - {output.tmp1} > {output.finalOut}
            gawk 'BEGIN{{FS="\t";OFS="\t"}}FNR==NR{{a[$1]=$2"\t"$3"\t"$4"\t"$5;next}}{{print $0,a[$5]}}' {output.finalOut} {input.prevDF} > {output.finalDFtmp}
            gawk 'BEGIN {{ FS = OFS = "\t" }} {{ for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = "*\t*\t*\t*" }}; 1' {output.finalDFtmp} > {output.finalDF}
            wc -l {params.mainDFdir}/* > {log}
            """

else:
    rule skipCDDResults:
        input:
            prevDF=rules.predict_cds.output.tsv, #rules.createStartingAnnotDF.output.mainProkka,
        output:
            finalDF=os.path.join(INTERIMPATH, "{id}", "annotation", "mainDF_AB.tsv"),
        message:
            "Skipping CDD annotation for {wildcards.id} as specified in the config file",
        conda:
            os.path.join(workflow.basedir, "envs","environment.yaml"),
        shell:
            """
            gawk 'BEGIN{{FS=OFS="\t"}} NR==1{{print $0, "NCBI_CDD", "NCBI_CDD_GeneID", "NCBI_CDD_Description", "NCBI_CDD_evalue"}} NR>1{{print $0, "*", "*", "*","*"}}' {input.prevDF} > {output.finalDF}
            """

if config["searches"]["cog"]:
    rule parseCOGResults:
        input:
            tblout=rules.searchCOGs.output.tab, # rules.searchCOGs.output.tblout,
            prevDF=rules.parseCDDResults.output.finalDF if config["searches"]["cdd"] else rules.skipCDDResults.output.finalDF,
            mapping=rules.prepareCOGdb.output.mapping,
        output:
            finalDF=os.path.join(INTERIMPATH, "{id}", "annotation", "mainDF_ABC.tsv"),
            ecutoff=temp(os.path.join(INTERIMPATH, "{id}", "cogs", "sequence_results_red_e_cutoff.tsv")),
            allResults=temp(os.path.join(INTERIMPATH, "{id}", "cogs", "All_COG.tsv")),
            tmp1=temp(os.path.join(INTERIMPATH, "{id}", "cogs", "temp1.tsv")),
            finalDFtmp=temp(os.path.join(INTERIMPATH, "{id}", "cogs", "cogs_final_tmp.tsv")),
            finalOut=temp(os.path.join(INTERIMPATH, "{id}", "cogs", "cogs_final.tsv")),
        message:
            "Parsing COG results for {wildcards.id}",
        log:
            os.path.join(LOGPATH, "{id}", "logs", "parseCOGResults.log")
        benchmark: 
            os.path.join(LOGPATH, "{id}", "benchmarks", "parseCOGResults.tsv")
        threads: 1
        priority: 10
        params:
            # mappingFile=config["mappingPaths"]["cog_mapping"],
            mainDFdir=lambda w, output:os.path.split(os.path.splitext(output[0])[0])[0],
        conda:
            os.path.join(workflow.basedir, "envs","environment.yaml"),
        shell:
            #format the full table and only select sequences above a certain e-value
            #get best hit based on bit score, and then e-value
            #merge with mapped names
            #add in an extra column that lists whether hits have a high confidence score
            #add in header
            #combine with previous dataframe
            #control lines
            """
            sed 's/ \\+ /\t/g' {input.tblout} | sed '/^#/d'| sed 's/ /\t/g' | sed 's/CDD://g' | gawk -F'\t' -v OFS='\t' '{{print $1, $2, $12, $11}}' | gawk -F'\t' -v OFS='\t' '($4 + 0) <= 1e-3'  > {output.ecutoff}
            sort -t$'\t' -k3,3gr -k4,4g {output.ecutoff} | sort -t$'\t' --stable -u -k1,1  | sort -t$'\t' -k3,3gr -k4,4g > {output.allResults}        
            LC_ALL=C join -a1 -1 2 -2 1 -e'*' -t $'\t' -o1.1,2.2,2.3,2.4,1.4 <(LC_ALL=C sort -k2 {output.allResults}) <(LC_ALL=C sort -k1 {input.mapping}) | LC_ALL=C  sort > {output.tmp1}        
            echo -e "gene_id\tNCBI_COG\tNCBI_COG_GeneID\tNCBI_COG_Description\tNCBI_COG_evalue" | cat - {output.tmp1} > {output.finalOut}
            gawk 'BEGIN{{FS="\t";OFS="\t"}}FNR==NR{{a[$1]=$2"\t"$3"\t"$4"\t"$5;next}}{{print $0,a[$5]}}' {output.finalOut} {input.prevDF} > {output.finalDFtmp}
            gawk 'BEGIN {{ FS = OFS = "\t" }} {{ for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = "*\t*\t*\t*" }}; 1' {output.finalDFtmp} > {output.finalDF}
            wc -l {params.mainDFdir}/* > {log}
            """

else:
    rule skipCOGResults:
        input:
            prevDF=rules.parseCDDResults.output.finalDF if config["searches"]["cdd"] else rules.skipCDDResults.output.finalDF,
        output:
            finalDF=os.path.join(RESULTPATH, "{id}", "annotation", "mainDF_ABC.tsv"),
        message:
            "Skipping COG annotation for {wildcards.id} as specified in the config file",
        conda:
            os.path.join(workflow.basedir, "envs","environment.yaml"),
        shell:
            """
            gawk 'BEGIN{{FS=OFS="\t"}} NR==1{{print $0, "NCBI_COG", "NCBI_COG_GeneID", "NCBI_COG_Description", "NCBI_COG_evalue"}} NR>1{{print $0, "*", "*", "*","*"}}' {input.prevDF} > {output.finalDF}
            """

if config["searches"]["arcog"]:
    rule parseArCOGResults:
        input:
            tblout=rules.searchArCOGs.output.tblout,
            prevDF=rules.parseCOGResults.output.finalDF if config["searches"]["cog"] else rules.skipCOGResults.output.finalDF,
            mapping=rules.prepareARCOGdbMapping.output.mapping, # "/mnt/DATA/common/AnnotateGenomesDBs/arCOG/arCOGdef.tab", # rules.prepareARCOGdb.output.mapping,
        output:
            finalDF=os.path.join(INTERIMPATH, "{id}", "annotation", "mainDF_ABCD.tsv"),
            ecutoff=temp(os.path.join(INTERIMPATH, "{id}", "arcogs", "sequence_results_red_e_cutoff.tsv")),
            allResults=temp(os.path.join(INTERIMPATH, "{id}", "arcogs", "All_arCOG.tsv")),
            tmp1=temp(os.path.join(INTERIMPATH, "{id}", "arcogs", "temp1.tsv")),
            tmp2=temp(os.path.join(INTERIMPATH, "{id}", "arcogs", "temp2.tsv")),
            finalDFtmp=temp(os.path.join(INTERIMPATH, "{id}", "arcogs", "arcogs_final_tmp.tsv")),
            finalOut=temp(os.path.join(INTERIMPATH, "{id}", "arcogs", "arcogs_final.tsv")),
        log:
            os.path.join(LOGPATH, "{id}", "logs", "parseArCOGResults.log")
        message:
            "Parsing arCOG results for {wildcards.id}",
        benchmark: 
            os.path.join(LOGPATH, "{id}", "benchmarks", "parseArCOGResults.tsv")
        threads: 1
        priority: 9
        params:
            mainDFdir=lambda w, output:os.path.split(os.path.splitext(output[0])[0])[0],
        conda:
            os.path.join(workflow.basedir, "envs","environment.yaml"),
        shell:
            #format the full table and only select sequences above a certain e-value
            #get best hit based on bit score, and then e-value
            #separate header with arcog and ID
            #merge with mapped names
            #add in an extra column that lists whether hits have a high confidence score
            #add in header
            #combine with previous dataframe
            #add "-" as empty character
            #control lines
            """
            sed 's/ \\+ /\t/g' {input.tblout} | sed '/^#/d'| sed 's/ /\t/g'| gawk -F'\t' -v OFS='\t' '{{print $1, $3, $5, $4}}' | gawk -F'\t' -v OFS='\t' '($4 + 0) <= 1e-3'  > {output.ecutoff}
            sort -t$'\t' -k3,3gr -k4,4g {output.ecutoff} | sort -t$'\t' --stable -u -k1,1  | sort -t$'\t' -k3,3gr -k4,4g > {output.allResults}        
            gawk -F'\t' -v OFS='\t' '{{split($2,a,"."); print $1, a[1], $3,$4}}' {output.allResults} | LC_ALL=C sort > {output.tmp1}
            LC_ALL=C join -a1 -1 2 -2 1 -e'*' -t $'\t' -o1.1,0,2.3,2.4,2.2,1.4 <(LC_ALL=C sort -k2 {output.tmp1}) <(LC_ALL=C sort -k1 {input.mapping}) | LC_ALL=C  sort > {output.tmp2}        
            echo -e "gene_id\tarcogs\tarcogs_geneID\tarcogs_Description\tarcogs_Pathway\tarcogs_evalue" | cat - {output.tmp2} > {output.finalOut}
            gawk 'BEGIN{{FS="\t";OFS="\t"}}FNR==NR{{a[$1]=$2"\t"$3"\t"$4"\t"$5"\t"$6;next}}{{print $0,a[$5]}}' {output.finalOut} {input.prevDF} > {output.finalDFtmp}
            gawk 'BEGIN {{ FS = OFS = "\t" }} {{ for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = "*\t*\t*\t*\t*" }}; 1' {output.finalDFtmp} > {output.finalDF}
            wc -l {params.mainDFdir}/* > {log}
            """

else:
    rule skipARCOGResults:
        input:
            prevDF=rules.parseCOGResults.output.finalDF if config["searches"]["cog"] else rules.skipCOGResults.output.finalDF,
        output:
            finalDF=os.path.join(INTERIMPATH, "{id}", "annotation", "mainDF_ABCD.tsv"),
        message:
            "Skipping arCOG annotation for {wildcards.id} as specified in the config file",
        conda:
            os.path.join(workflow.basedir, "envs","environment.yaml"),
        shell:
            """
            gawk 'BEGIN{{FS=OFS="\t"}} NR==1{{print $0, "arcogs", "arcogs_Description", "arcogs_Pathway", "arcogs_evalue"}} NR>1{{print $0, "*", "*", "*", "*"}}' {input.prevDF} > {output.finalDF}
            """

if config["searches"]["pgap"]:
    rule parsePGAPResults:
        input:
            tblout=rules.searchPGAP.output.tblout,
            prevDF=rules.parseArCOGResults.output.finalDF if config["searches"]["arcog"] else rules.skipARCOGResults.output.finalDF,
            mapping=rules.preparePGAPdb.output.mapping,
        output:
            finalDF=os.path.join(RESULTPATH, "{id}", "annotation", "{id}_finalAnnotation.tsv"),
            ecutoff=temp(os.path.join(INTERIMPATH, "{id}", "pgap", "sequence_results_red_e_cutoff.tsv")),
            allResults=temp(os.path.join(INTERIMPATH, "{id}", "pgap", "All_pgap.tsv")),
            tmp1=temp(os.path.join(INTERIMPATH, "{id}", "pgap", "temp1.tsv")),
            finalDFtmp=temp(os.path.join(INTERIMPATH, "{id}", "pgap", "pgap_final_tmp.tsv")),
            finalOut=temp(os.path.join(INTERIMPATH, "{id}", "pgap", "pgap_final.tsv")),
        log:
            os.path.join(LOGPATH, "{id}", "logs", "parsePGAPResults.log")
        message:
            "Parsing PGAP results for {wildcards.id}",
        benchmark: 
            os.path.join(LOGPATH, "{id}", "benchmarks", "parsePGAPResults.tsv")
        threads: 1
        priority: 6
        params:
            # mappingFile=config["mappingPaths"]["pgap_mapping"],
            mainDFdir=lambda w, output:os.path.split(os.path.splitext(output[0])[0])[0],
        conda:
            os.path.join(workflow.basedir, "envs","environment.yaml"),
        shell:
            #format the full table and only select sequences above a certain e-value
            #get best hit based on bit score, and then e-value
            #merge with pgap names
            #add in header
            #combine with previous dataframe
            #control lines
            """
            sed 's/ \\+ /\t/g' {input.tblout} | sed '/^#/d'| sed 's/ /\t/g'| gawk -F'\t' -v OFS='\t' '{{print $1, $3, $5, $4}}' | gawk -F'\t' -v OFS='\t' '($4 + 0) <= 1e-3'  > {output.ecutoff}
            sort -t$'\t' -k3,3gr -k4,4g {output.ecutoff} | sort -t$'\t' --stable -u -k1,1  | sort -t$'\t' -k3,3gr -k4,4g > {output.allResults}
            LC_ALL=C join -a1 -1 2 -2 1 -e'*' -t $'\t' -o1.1,0,2.11,2.13,1.4  <(LC_ALL=C sort -k2 {output.allResults}) <(LC_ALL=C sort -k1 {input.mapping}) | LC_ALL=C sort > {output.tmp1}
            echo -e "gene_id\tPGAP\tPGAP_description\tPGAP_EC\tPGAP_Evalue" | cat - {output.tmp1} > {output.finalOut}
            gawk 'BEGIN{{FS="\t";OFS="\t"}}FNR==NR{{a[$1]=$2"\t"$3"\t"$4"\t"$5;next}}{{print $0,a[$5]}}' {output.finalOut} {input.prevDF} > {output.finalDFtmp}
            gawk 'BEGIN {{ FS = OFS = "\t" }} {{ for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = "*\t*\t*\t*" }}; 1' {output.finalDFtmp} > {output.finalDF}
            wc -l {params.mainDFdir}/* > {log}
            """
else:
    rule skipPGAPResults:
        input:
            prevDF=rules.parseArCOGResults.output.finalDF if config["searches"]["arcog"] else rules.skipARCOGResults.output.finalDF,
        output:
            finalDF=os.path.join(RESULTPATH, "{id}", "annotation", "{id}_finalAnnotation.tsv"),
        message:
            "Skipping PGAP annotation for {wildcards.id} as specified in the config file",
        conda:
            os.path.join(workflow.basedir, "envs","environment.yaml"),
        shell:
            """
            gawk 'BEGIN{{FS=OFS="\t"}} NR==1{{print $0, "PGAP", "PGAP_description", "PGAP_EC", "PGAP_Evalue"}} NR>1{{print $0, "*", "*", "*", "*"}}' {input.prevDF} > {output.finalDF}
            """

rule collectMasterTable:
    input:
        finalDF= expand(os.path.join(RESULTPATH, "{id}", "annotation", "{id}_finalAnnotation.tsv"), id=NAMES),
    output:
        finalAnnotation=os.path.join(RESULTPATH, "common", "annotation", "finalAnnotation.tsv")
    message:
        "Collecting all annotations into a single table"
    conda:
        os.path.join(workflow.basedir, "envs","environment.yaml"),
    shell:
        """
        gawk 'FNR==1 && NR!=1 {{next}} {{print}}' {input.finalDF} > {output.finalAnnotation}
        """
