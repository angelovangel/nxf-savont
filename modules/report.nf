process MAKE_REPORT {
    container 'amancevice/pandas:slim'
    publishDir "${params.outdir}", mode: 'copy'

    input:
        path "raw_stats/*"
        path "norm_stats/*"
        path "filtered/*"
        path "taxonomy/*"
        path "lineage/*"
        path "combined_species.tsv"

    output:
        path "nxf-savont-report.html"

    script:
    def wf_cmd = workflow.commandLine.replace('"', '\\"')
    """
    echo -e "sample\\traw_reads\\traw_n50\\traw_Q20\\tnormalised_reads\\tfiltered_reads\\tasvs" > summary_counts.tsv
    
    # Iterate through raw stats files
    for raw_stats in raw_stats/*.readstats.tsv; do
        sample=\$(basename \$raw_stats .readstats.tsv)
        
        # Extract count and N50 from raw stats
        raw_data=\$(awk 'NR==2' \$raw_stats)
        raw_count=\$(echo "\$raw_data" | awk '{print \$2}')
        raw_n50=\$(echo "\$raw_data" | awk '{print \$11}')
        raw_q20=\$(echo "\$raw_data" | awk '{print \$12}')
        
        # Check if norm/subsampled count exists
        if ls norm_stats/\${sample}.*readstats.tsv 1> /dev/null 2>&1; then
            norm_file=\$(ls norm_stats/\${sample}.*readstats.tsv | head -n 1)
            norm_count=\$(awk 'NR==2' "\$norm_file" | awk '{print \$2}')
        else
            norm_count=\$raw_count
        fi

        # Check if filtered count exists (from SAVONT_ASV)
        if [ -f "filtered/\${sample}.filtered_reads.txt" ]; then
            filtered_count=\$(awk 'NR==2' "filtered/\${sample}.filtered_reads.txt" | awk '{print \$1}')
            asvs=\$(awk 'NR==2' "filtered/\${sample}.filtered_reads.txt" | awk '{print \$2}')
        else
            filtered_count=\$raw_count
            asvs=0
        fi
        
        echo -e "\$sample\\t\$raw_count\\t\$raw_n50\\t\$raw_q20\\t\$norm_count\\t\$filtered_count\\t\$asvs" >> summary_counts.tsv
    done

    # Write workflow info CSV
    cat > wfinfo.csv << EOF
    Pipeline,${workflow.manifest.name ?: 'nxf-savont'}
    Workflow version,${workflow.manifest.version ?: 'dev'}
    Workflow revision,${workflow.scriptId.take(10)}
    Run name,${workflow.runName}
    Nextflow version,${nextflow.version}
    Command line,"${wf_cmd}"
    EOF

    # Gather readstats files for rarefaction (use norm_stats if exist, else raw_stats)
    READSTATS_FILES=""
    for raw_stats in raw_stats/*.readstats.tsv; do
        sample=\$(basename \$raw_stats .readstats.tsv)
        if ls norm_stats/\${sample}.*readstats.tsv 1> /dev/null 2>&1; then
            norm_file=\$(ls norm_stats/\${sample}.*readstats.tsv | head -n 1)
            READSTATS_FILES="\$READSTATS_FILES \$norm_file"
        else
            READSTATS_FILES="\$READSTATS_FILES \$raw_stats"
        fi
    done

    # Generate rarefaction JSON
    rarefaction.py --lineage combined_species.tsv --readstats \$READSTATS_FILES --rarefaction rarefaction.json

    # Generate HTML report
    make_html_report.py --summary summary_counts.tsv --combined combined_species.tsv --abundances taxonomy/*.tsv --lineages lineage/*.tsv --wfinfo wfinfo.csv --rarefaction rarefaction.json
    """
}

process VERSIONS {
    publishDir "${params.outdir}/logs", mode: 'copy'
    container 'docker.io/aangeloo/nxf-tgs:latest'
    
    input:
    path 'versions*.txt'
    path summary_file
    
    output:
    path "nxf-savont-execution-summary.txt"

    script:
    """
    echo -e "\nSoftware Versions" > software_versions.txt
    echo -e "-----------------------------------------" >> software_versions.txt
    
    for f in \$(ls versions*.txt | sort -V); do
        awk -F ': ' '{ if (NF>1) printf "%-25s: %s\\n", \$1, \$2; else print \$0 }' \$f >> software_versions.txt
    done
    
    cat $summary_file software_versions.txt > nxf-savont-execution-summary.txt
    """
}
