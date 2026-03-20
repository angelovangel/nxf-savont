#!/usr/bin/env nextflow

process CONVERT_EXCEL {
    container 'docker.io/aangeloo/nxf-tgs:latest'

    input:
        path(excelfile)

    output:
        path("*.csv")
        path "versions.txt", emit: versions

    script:
    """
    convert_excel.R $excelfile

    cat <<EOF > versions.txt
    ${task.process}: R \$(R --version | head -n 1 | awk '{print \$3}')
    EOF
    """
}

// takes in csv, checks for duplicate barcodes etc.
process VALIDATE_SAMPLESHEET {
    container 'docker.io/aangeloo/nxf-tgs:latest'
    //publishDir "$params.outdir", mode: 'copy', pattern: '00-validated-samplesheet.csv'

    input: 
    path(csv)

    output:
    path("samplesheet-validated.csv")
    path "versions.txt", emit: versions


    script:
    """
    validate_samplesheet.R $csv

    cat <<EOF > versions.txt
    ${task.process}: R \$(R --version | head -n 1 | awk '{print \$3}')
    EOF
    """
}

process MERGE_READS {
    container 'docker.io/aangeloo/nxf-tgs:latest'
    errorStrategy {
        if (task.exitStatus == 42) {
            println "MERGE_READS: [WARNING] ${samplename} (${barcode}) defined in the samplesheet but does not exist in the data. Ignoring."
            return 'ignore'
        }
        return 'ignore'
    }
    tag "${samplename} - ${barcode}"

    //publishDir "$params.outdir/00-basecall/processed", mode: 'copy', pattern: '*{fastq.gz,fastq,bam}'

    input:
        tuple val(samplename), val(barcode), path(bam_pass)
    
    output: 
        path('*{fastq.gz,fastq,bam}')
        path "versions.txt", emit: versions
    
    script:
    """
    if [ -d "${bam_pass}/${barcode}" ]; then
        samtools cat ${bam_pass}/${barcode}/*.bam > ${samplename}.bam
    else
        exit 42
    fi

    cat <<EOF > versions.txt
    ${task.process}: samtools \$(samtools --version | head -n 1 | awk '{print \$2}')
    EOF
    """
}

process CONVERT_READS {
    container 'docker.io/aangeloo/nxf-tgs:latest'
    publishDir "${params.outdir}/00-basecall", mode: 'copy', pattern: '*{fastq.gz,fastq}'
    tag "${reads.simpleName}, ${reads.extension} file"
    
    input:
        path reads

    output:
        path "*fastq.gz", includeInputs: true
        path "versions.txt", emit: versions
    
    script:
    """
    if [[ ${reads.extension} == "bam" ]]; then
        samtools fastq --threads ${task.cpus} ${reads} > ${reads.simpleName}.fastq.gz
    elif [[ ${reads.extension} == "fastq" || ${reads.extension} == "fq" ]]; then
        pigz -c ${reads} > ${reads.simpleName}.fastq.gz
    fi

    cat <<EOF > versions.txt
    ${task.process}: samtools \$(samtools --version | head -n 1 | awk '{print \$2}')
    ${task.process}: pigz \$(pigz --version 2>&1 | awk '{print \$2}')
    EOF
    """
}

process READ_STATS {
    container 'docker.io/aangeloo/nxf-tgs:latest'
    publishDir "${params.outdir}/00-basecall/readqc", mode: 'copy', pattern: '*readstats.tsv'
    tag "${reads.simpleName}"

    input:
        path reads

    output:
        path "*readstats.tsv", emit: stats
        path "versions.txt", emit: versions

    script:
    
    """
    faster -t ${reads} > ${reads.simpleName}.readstats.tsv

    cat <<EOF > versions.txt
    ${task.process}: faster \$(faster --version 2>&1 | awk '{print \$2}')
    EOF
    """
}

process SUBSAMPLE_READS {
    container 'docker.io/aangeloo/nxf-tgs:latest'
    tag "${reads.simpleName}"
    
    input:
        path reads
    
    output:
        path "*.fastq.gz", emit: reads
        path "*.subsampled.readstats.tsv", emit: stats
        path "versions.txt", emit: versions
    
    script:
    """
    faster -p ${params.subsample} ${reads} > ${reads.simpleName}.subsampled.fastq
    pigz -c ${reads.simpleName}.subsampled.fastq > ${reads.simpleName}.subsampled.fastq.gz
    
    # Get stats
    faster -t ${reads.simpleName}.subsampled.fastq.gz > ${reads.simpleName}.subsampled.readstats.tsv

    cat <<EOF > versions.txt
    ${task.process}: faster \$(faster --version 2>&1 | awk '{print \$2}')
    EOF
    """
}

process GET_MIN_READS {
    container 'docker.io/aangeloo/nxf-tgs:latest'
    
    input:
        path stats_files
        
    output:
        env MIN_READS, emit: min_reads
        
    script:
    """
    MIN_READS=\$(awk 'FNR>1 {print \$2}' *.readstats.tsv | sort -n | head -n1)
    """
}

process NORMALIZE_READS {
    container 'docker.io/aangeloo/nxf-tgs:latest'
    publishDir "${params.outdir}/00-basecall/readqc", mode: 'copy', pattern: '*norm.readstats.tsv'
    tag "${reads.simpleName}"
    
    input:
        path reads
        val min_reads
    
    output:
        path "*.norm.fastq.gz", emit: reads
        path "*.norm.readstats.tsv", emit: stats
        path "versions.txt", emit: versions
    
    script:
    """
    # Get current read count
    CURRENT_READS=\$(faster -t ${reads} | awk 'FNR>1 {print \$2}')
    
    # Calculate proportion (min / current) reliably using POSIX locale
    PROP=\$(LC_ALL=C awk -v min=${min_reads} -v curr=\$CURRENT_READS 'BEGIN { printf "%.6f", min / curr }')
    
    seqkit sample -p \${PROP} ${reads} > ${reads.simpleName}.norm.fastq
    pigz -c ${reads.simpleName}.norm.fastq > ${reads.simpleName}.norm.fastq.gz
    
    # Get stats for normalized reads
    faster -t ${reads.simpleName}.norm.fastq.gz > ${reads.simpleName}.norm.readstats.tsv
    
    cat <<EOF > versions.txt
    ${task.process}: seqkit \$(seqkit -h | grep "Version: " | awk '{print \$2}')
    EOF
    """
}

process READ_HIST {
    container 'docker.io/aangeloo/nxf-tgs:latest'
    tag "${reads.simpleName}"

    input:
        path reads

    output:
        path "*.hist"
        path "versions.txt", emit: versions

    script:
    """
    if [[ ${reads.extension} == "bam" ]]; then
        samtools fastq ${reads} | faster2 --len -  | bincount.awk -v type=len  | sort -n > ${reads.simpleName}.len.hist
        samtools fastq ${reads} | faster2 --gc -   | bincount.awk -v type=gc   | sort -n > ${reads.simpleName}.gc.hist
        samtools fastq ${reads} | fasterplot -q - | sort -n > ${reads.simpleName}.qual.hist
    else
        faster2 --len ${reads}  | bincount.awk -v type=len  | sort -n > ${reads.simpleName}.len.hist
        faster2 --gc ${reads}   | bincount.awk -v type=gc   | sort -n > ${reads.simpleName}.gc.hist
        fasterplot -q ${reads} | sort -n > ${reads.simpleName}.qual.hist
    fi

    cat <<EOF > versions.txt
    ${task.process}: samtools \$(samtools --version | head -n 1 | awk '{print \$2}')
    ${task.process}: faster2 \$(faster2 --version 2>&1 | awk '{print \$2}')
    ${task.process}: fasterplot \$(fasterplot --version 2>&1 | awk '{print \$2}')
    EOF
    """
}
