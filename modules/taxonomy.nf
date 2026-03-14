
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
    """
}

process SAVONT_CLASSIFY {
    container 'docker.io/aangeloo/nxf-savont:latest'
    tag "${asvs.simpleName}"
    publishDir "${params.outdir}/01-taxonomy", mode: 'copy'

    input:
        path asvs

    output:
        path "${asvs.simpleName}_rel-abundance.tsv", emit: ch_abundance

    script:
    def db = params.db == "emu" ? "emu_default" : "silva_db"
    """
    if [ -s "${asvs}/final_asvs.fasta" ]; then
        savont classify -i $asvs -o . --emu-db /databases/$db -t $task.cpus
        mv species_abundance.tsv ${asvs.simpleName}_rel-abundance.tsv
    else
        # Create empty abundance file with headers if no ASVs found
        echo -e "abundance\\ttax_id\\tspecies\\tgenus\\tfamily\\torder\\tclass\\tphylum\\tclade\\tsuperkingdom" > ${asvs.simpleName}_rel-abundance.tsv
    fi
    """
}

process SAVONT_COMBINE {
    container 'docker.io/aangeloo/nxf-tgs:latest'
    publishDir "${params.outdir}/01-taxonomy", mode: 'copy'

    input:
        path "abundances/*"

    output:
        path "savont-combined-species.tsv", emit: ch_combined

    script:
    """
    combine_savont.py savont-combined-species.tsv abundances/*.tsv
    """
}

/* process EMU_ABUNDANCE {
    container 'docker.io/aangeloo/emu:latest'
    //errorStrategy 'ignore'
    //maxForks 1 // avoid memory issues due to cpus x forks
    tag "${reads.name}"
    publishDir "${params.outdir}/01-taxonomy", mode: 'copy', pattern: '**rel-abundance.tsv'

    input:
        path reads

    output:
        path "*rel-abundance.tsv", emit: ch_abundance
        path "*.mapping_stats.tsv", emit: counts

    script:
    """
    # Run EMU but allow it to fail
    set +e
    emu abundance --db \$EMU_DATABASE_DIR/${params.db} --keep-counts --N ${params.N} --K ${params.K} --output-dir . --threads $task.cpus $reads
    EMU_EXIT=\$?
    set -e

    if [ "\$EMU_EXIT" -ne 0 ]; then
        # Check if the error was due to 0 mapped reads by capturing the log or just by checking if the output file is missing.
        if [ ! -f "${reads.simpleName}_rel-abundance.tsv" ]; then
            # Create a dummy abundance file representing 0 mapped reads
            echo -e "tax_id\tabundance\testimated counts\ttax_name\tsuperkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies" > ${reads.simpleName}_rel-abundance.tsv
            
            # Since 0 mapped reads, set mapped to 0 and unmapped to total sequences if possible, or just mock it all as 0 for simplicity.
            mapped=0
            unmapped=0
            mapped_filtered=0
            mapped_unclassified=0
            
            # Since there is no easy way to get total reads without counting the fastq here, we'll set unmapped to 0 too or count it.
            # Counting reads in fastq (assuming 4 lines per read)
            if [[ "$reads" == *.gz ]]; then
                total_reads=\$(zcat "1.fastq.gz" | wc -l | awk '{print \$1/4}')
            else
                total_reads=\$(wc -l < "$reads" | awk '{print \$1/4}')
            fi
            unmapped=\$total_reads
            
            # Write out the three required lines if EMU expects them or to match the downstream script
            echo -e "unmapped\t1.0\t\$unmapped\tunmapped\t\t\t\t\t\t\t" >> ${reads.simpleName}_rel-abundance.tsv
            echo -e "mapped_filtered\t0.0\t0\tmapped_filtered\t\t\t\t\t\t\t" >> ${reads.simpleName}_rel-abundance.tsv
            echo -e "mapped_unclassified\t0.0\t0\tmapped_unclassified\t\t\t\t\t\t\t" >> ${reads.simpleName}_rel-abundance.tsv
        else
            # Some other error happened, fail properly
            exit \$EMU_EXIT
        fi
    else
        REL_ABUNDANCE=\$(ls *rel-abundance.tsv)
        
        # Get total mapped (estimated counts column is the last column)
        mapped=\$(head -n -3 \$REL_ABUNDANCE | tail -n +2 | awk '{sum+=\$NF} END {print sum}')
        unmapped=\$(grep "^unmapped" \$REL_ABUNDANCE | awk '{print \$NF}')
        mapped_filtered=\$(grep "^mapped_filtered" \$REL_ABUNDANCE | awk '{print \$NF}')
        mapped_unclassified=\$(grep "^mapped_unclassified" \$REL_ABUNDANCE | awk '{print \$NF}')
    fi

    echo -e "mapped\\tunmapped\\tmapped_filtered\\tmapped_unclassified" > ${reads.simpleName}.mapping_stats.tsv
    echo -e "\$mapped\\t\$unmapped\\t\$mapped_filtered\\t\$mapped_unclassified" >> ${reads.simpleName}.mapping_stats.tsv
    """
}

process EMU_COMBINE {
    container 'docker.io/aangeloo/emu:latest'
    publishDir "${params.outdir}/01-taxonomy", mode: 'copy', pattern: 'emu-combined-species.tsv'

    input:
        path abundance_dir

    output:
        path "emu-combined-species.tsv", emit: ch_combined

    script:
    """
    emu combine-outputs . species
    """
} */