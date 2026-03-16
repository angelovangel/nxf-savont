process MAKE_REPORT {
    container 'docker.io/aangeloo/nxf-tgs:latest'
    publishDir "${params.outdir}", mode: 'copy'

    input:
        path "raw_stats/*"
        path "filtered/*"
        path "taxonomy/*"
        path "combined_species.tsv"

    output:
        path "nxf-savont-report.html"

    script:
    def wf_cmd = workflow.commandLine.replace('"', '\\"')
    """
    echo -e "sample\\traw_reads\\traw_n50\\traw_Q20\\tfiltered_reads\\tasvs" > summary_counts.tsv
    
    # Iterate through raw stats files
    for raw_stats in raw_stats/*.readstats.tsv; do
        sample=\$(basename \$raw_stats .readstats.tsv)
        
        # Extract count and N50 from raw stats
        raw_data=\$(awk 'NR==2' \$raw_stats)
        raw_count=\$(echo "\$raw_data" | awk '{print \$2}')
        raw_n50=\$(echo "\$raw_data" | awk '{print \$11}')
        raw_q20=\$(echo "\$raw_data" | awk '{print \$12}')
        
        # Check if filtered count exists (from SAVONT_ASV)
        if [ -f "filtered/\${sample}.filtered_reads.txt" ]; then
            filtered_count=\$(awk 'NR==2' "filtered/\${sample}.filtered_reads.txt" | awk '{print \$1}')
            asvs=\$(awk 'NR==2' "filtered/\${sample}.filtered_reads.txt" | awk '{print \$2}')
        else
            filtered_count=\$raw_count
            asvs=0
        fi
        
        echo -e "\$sample\\t\$raw_count\\t\$raw_n50\\t\$raw_q20\\t\$filtered_count\\t\$asvs" >> summary_counts.tsv
    done

    # Write workflow info CSV
    cat > wfinfo.csv << EOF
    Pipeline,${workflow.manifest.name ?: 'nxf-savont'}
    Version,${workflow.manifest.version ?: 'dev'}
    Nextflow version,${nextflow.version}
    Workflow script ID,${workflow.scriptId.take(10)}
    Run name,${workflow.runName}
    Command line,"${wf_cmd}"
    EOF

    # Generate HTML report
    make_html_report.py --summary summary_counts.tsv --combined combined_species.tsv --abundances taxonomy/*.tsv --wfinfo wfinfo.csv
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
