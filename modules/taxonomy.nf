
process SAVONT_ASV {
    container 'docker.io/aangeloo/nxf-savont:latest'
    tag "${reads.name}"
    publishDir "${params.outdir}/01-taxonomy/asvs", mode: 'copy', pattern: "*.asvs.fasta"

    input:
        path reads

    output:
        path "${reads.simpleName}.asvs.fasta", optional: true
        path "${reads.simpleName}", emit: ch_asvs, type: 'directory'
        path "${reads.simpleName}.filtered_reads.txt", emit: counts
        path "versions.txt", emit: versions

    script:
    """
    savont asv $reads -o ${reads.simpleName} -t $task.cpus --min-read-length ${params.minreadlength} --max-read-length ${params.maxreadlength}
    
    # Copy final_asvs.fasta to task root with sample name for publishDir
    [ -f "${reads.simpleName}/final_asvs.fasta" ] && cp "${reads.simpleName}/final_asvs.fasta" "${reads.simpleName}.asvs.fasta" || true
    
    # Extract filtered read counts from log
    log_file=\$(ls ${reads.simpleName}/savont_*.log | head -n 1)
    echo "filtered_reads\\tasvs" > ${reads.simpleName}.filtered_reads.txt
    if [ -f "\$log_file" ]; then
        counts=\$(grep "Number of valid reads" "\$log_file" | awk -F' - ' '{print \$2}' | awk -F'.' '{print \$1}' | sed 's/ //g')
        asvs=\$(grep "Final consensus count after EM refinement:" "\$log_file" | awk -F': ' '{print \$2}' | awk -F'.' '{print \$1}' | sed 's/ //g')
        echo "\$counts\\t\$asvs" >> ${reads.simpleName}.filtered_reads.txt
    else
        echo "0\\t0" >> ${reads.simpleName}.filtered_reads.txt
    fi

    cat <<EOF > versions.txt
    ${task.process}: savont \$(savont --version 2>&1 | sed 's/savont //')
    EOF
    """
}

process SAVONT_CLASSIFY {
    container 'docker.io/aangeloo/nxf-savont:latest'
    tag "${asvs.simpleName}"
    publishDir "${params.outdir}/01-taxonomy", mode: 'copy', pattern: "*.tsv"

    input:
        path asvs

    output:
        path "${asvs.simpleName}_rel-abundance.tsv", emit: ch_abundance
        path "versions.txt", emit: versions

    script:
    def db = params.db == "emu" ? "--emu-db /databases/emu_default" : "--silva-db /databases/silva_db"
    """
    if [ -s "${asvs}/final_asvs.fasta" ]; then
        savont classify -i $asvs -o . $db -t $task.cpus --species-threshold ${params.species_threshold} --genus-threshold ${params.genus_threshold}
        mv species_abundance.tsv ${asvs.simpleName}_rel-abundance.tsv
    else
        # Create empty abundance file with headers if no ASVs found
        echo -e "abundance\\ttax_id\\tspecies\\tgenus\\tfamily\\torder\\tclass\\tphylum\\tclade\\tsuperkingdom" > ${asvs.simpleName}_rel-abundance.tsv
    fi

    cat <<EOF > versions.txt
    ${task.process}: savont \$(savont --version 2>&1 | sed 's/savont //')
    EOF
    """
}

process SAVONT_COMBINE {
    container 'docker.io/aangeloo/nxf-tgs:latest'
    publishDir "${params.outdir}/01-taxonomy", mode: 'copy', pattern: "*.tsv"

    input:
        path "abundances/*"

    output:
        path "savont-combined-species.tsv", emit: ch_combined
        path "versions.txt", emit: versions

    script:
    """
    combine_savont.py savont-combined-species.tsv abundances/*.tsv

    cat <<EOF > versions.txt
    ${task.process}: python \$(python3 --version | awk '{print \$2}')
    EOF
    """
}

process TAXONKIT_LINEAGE {
    container 'docker.io/aangeloo/nxf-savont:latest'
    publishDir "${params.outdir}/01-taxonomy/lineage", mode: 'copy', pattern: "*.tsv"

    input:
        path abundance_file

    output:
        path "${abundance_file.simpleName}_lineage.tsv", emit: ch_lineage
        path "versions.txt", emit: versions

    script:
    """
    cut -f1,2 $abundance_file | taxonkit lineage -R --data-dir /opt/taxonkit > ${abundance_file.simpleName}_lineage.tsv
    
    cat <<EOF > versions.txt
    ${task.process}: taxonkit \$(taxonkit -h | grep Version 2>&1 | sed 's/Version: //')
    EOF
    """
}
