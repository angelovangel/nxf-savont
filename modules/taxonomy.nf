
process EMU_ABUNDANCE {
    container 'docker.io/aangeloo/emu:latest'
    //errorStrategy 'ignore'
    tag "${reads.name}"
    publishDir "${params.outdir}/01-taxonomy", mode: 'copy', pattern: '**rel-abundance.tsv'

    input:
        path reads

    output:
        path "*rel-abundance.tsv", emit: ch_abundance
        path "*.mapping_stats.tsv", emit: counts

    script:
    """
    emu abundance --db \$EMU_DATABASE_DIR/${params.db} --keep-counts --output-dir . --threads $task.cpus $reads
    
    # Extract stats from rel-abundance.tsv
    # mapped: sum of estimated counts (excluding the last 3 rows)
    # unmapped, mapped_filtered, mapped_unclassified: from the last 3 rows
    
    REL_ABUNDANCE=\$(ls *rel-abundance.tsv)
    
    # Get total mapped (estimated counts column is the last column)
    mapped=\$(head -n -3 \$REL_ABUNDANCE | tail -n +2 | awk '{sum+=\$NF} END {print sum}')
    unmapped=\$(grep "^unmapped" \$REL_ABUNDANCE | awk '{print \$NF}')
    mapped_filtered=\$(grep "^mapped_filtered" \$REL_ABUNDANCE | awk '{print \$NF}')
    mapped_unclassified=\$(grep "^mapped_unclassified" \$REL_ABUNDANCE | awk '{print \$NF}')
    
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
}