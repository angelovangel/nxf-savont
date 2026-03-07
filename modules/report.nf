process MAKE_REPORT {
    container 'docker.io/aangeloo/nxf-tgs:latest'
    publishDir "${params.outdir}", mode: 'copy'

    input:
        path "raw/*"
        path "filtered/*"
        path "mapping/*"
        path "taxonomy/*"

    output:
        path "summary_counts.tsv"
        path "report.html"

    script:
    """
    echo -e "sample\\traw_reads\\tfiltered_reads\\tmapped\\tunmapped\\tmapped_filtered\\tmapped_unclassified" > summary_counts.tsv
    
    # Iterate through raw read count files
    for raw_file in raw/*.reads.txt; do
        sample=\$(basename \$raw_file .reads.txt)
        
        raw_count=\$(cat \$raw_file)
        
        # Check if filtered count exists
        if [ -f "filtered/\${sample}.filtered_reads.txt" ]; then
            filtered_count=\$(cat "filtered/\${sample}.filtered_reads.txt")
        else
            filtered_count=\$raw_count
        fi
        
        # Check if mapping stats exist
        # Emu might have different naming if filtered (sample.filtered.mapping_stats.tsv)
        mapping_file=\$(ls mapping/\${sample}*.mapping_stats.tsv 2>/dev/null | head -n 1)
        
        if [ -n "\$mapping_file" ]; then
            mapping_data=\$(tail -n 1 \$mapping_file)
        else
            mapping_data="0\\t0\\t0\\t0"
        fi
        
        echo -e "\$sample\\t\$raw_count\\t\$filtered_count\\t\$mapping_data" >> summary_counts.tsv
    done

    # Generate HTML report
    make_html_report.py summary_counts.tsv taxonomy/*.tsv
    """
}
