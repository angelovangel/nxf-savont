#!/usr/bin/env nextflow

process CONVERT_EXCEL {
    container 'docker.io/aangeloo/nxf-tgs:latest'

    input:
        path(excelfile)

    output:
        path("*.csv")

    script:
    """
    convert_excel.R $excelfile
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


    script:
    """
    validate_samplesheet.R $csv
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
    
    script:
    """
    if [ -d "${bam_pass}/${barcode}" ]; then
        samtools cat ${bam_pass}/${barcode}/*.bam > ${samplename}.bam
    else
        exit 42
    fi
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
    
    script:
    """
    if [[ ${reads.extension} == "bam" ]]; then
        samtools fastq --threads ${task.cpus} ${reads} > ${reads.simpleName}.fastq.gz
    elif [[ ${reads.extension} == "fastq" || ${reads.extension} == "fq" ]]; then
        pigz -c ${reads} > ${reads.simpleName}.fastq.gz
    fi
    """
}

process READ_STATS {
    container 'docker.io/aangeloo/nxf-tgs:latest'
    publishDir "${params.outdir}/00-basecall", mode: 'copy', pattern: '*readstats.tsv'
    tag "${reads.simpleName}"

    input:
        path reads

    output:
        path "*readstats.tsv"
        path "*.reads.txt", emit: counts

    script:
    
    """
    faster -t ${reads} > ${reads.simpleName}.readstats.tsv
    awk 'NR==2 {print \$2}' ${reads.simpleName}.readstats.tsv > ${reads.simpleName}.reads.txt
    """
}

process FILTER_READS {
    container 'docker.io/aangeloo/nxf-tgs:latest'
    tag "${reads.simpleName}"
    //publishDir "${params.outdir}/00-basecall/filtered", mode: 'copy', pattern: '*readstats.tsv'
    errorStrategy {
        if (task.exitStatus == 42) {
            println "FILTER_READS: [WARNING] ${reads} has no reads after filtering. Ignoring."
            return 'ignore'
        }
        return 'ignore'
    }
    
    input:
        path reads

    output:
        path "*fastq.gz", emit: ch_filtered_reads
        path "*readstats.tsv"
        path "*.filtered_reads.txt", emit: counts
    
    script:
    """
    faster --filterl 1000 $reads | faster --filterl -2000 - > ${reads.simpleName}.filtered.fastq
    if [ \$(wc -l < ${reads.simpleName}.filtered.fastq) -eq 0 ]; then
        exit 42
    fi
    pigz -c ${reads.simpleName}.filtered.fastq > ${reads.simpleName}.filtered.fastq.gz
    faster -t ${reads.simpleName}.filtered.fastq.gz > ${reads.simpleName}.filtered.readstats.tsv
    awk 'NR==2 {print \$2}' ${reads.simpleName}.filtered.readstats.tsv > ${reads.simpleName}.filtered_reads.txt
    """
}

process READ_HIST {
    container 'docker.io/aangeloo/nxf-tgs:latest'
    tag "${reads.simpleName}"

    input:
        path reads

    output:
        path "*.hist"

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
    """
}
