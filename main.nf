include {DORADO_BASECALL; DORADO_BASECALL_BARCODING} from "./modules/basecall.nf"
include {MERGE_READS; READ_STATS; READ_HIST; CONVERT_EXCEL; VALIDATE_SAMPLESHEET; CONVERT_READS; SUBSAMPLE_READS} from './modules/reads.nf'
include {SAVONT_ASV; SAVONT_CLASSIFY; SAVONT_COMBINE; TAXONKIT_LINEAGE} from './modules/taxonomy.nf'

include {MAKE_REPORT; VERSIONS} from './modules/report.nf'

if (params.kit && !params.samplesheet) {
    error "If --kit is specified, --samplesheet must also be provided."
}

if (!params.kit && params.samplesheet) {
    error "If --samplesheet is provided, --kit must also be specified."
}

// Log execution environment
def summary = """
    NXF-SAVONT Execution Summary
    =========================================
    Start                    : ${workflow.start.format('yyyy-MM-dd HH:mm:ss')}
    Profile                  : ${workflow.profile}
    Container Engine         : ${workflow.containerEngine ?: 'local'}
    Workflow version         : ${workflow.manifest.version}
    Workflow script ID       : ${workflow.scriptId.take(10)}
    NXF version              : ${workflow.nextflow.version}
    -----------------------------------------
    Parameters
    -----------------------------------------
    """.stripIndent()

params.each { name, value -> summary += "${name.padRight(25)}: ${value}\n" }
summary += "-----------------------------------------"

def out_dir = file(params.outdir)
if( !out_dir.exists() ) out_dir.mkdirs()
def ch_summary_file = Channel.of(summary)
    .collectFile(name: 'execution-summary.txt', newLine: false)


workflow basecall {
    ch_pod5 = Channel.fromPath(params.pod5, checkIfExists: true)
    ch_samplesheet = params.samplesheet ? Channel.fromPath(params.samplesheet, checkIfExists: true) : null
    ch_versions = Channel.empty()
    
    if (params.kit) {    
        if (params.samplesheet.endsWith('.xlsx')) {
            CONVERT_EXCEL(ch_samplesheet)
            ch_samplesheet_validated = VALIDATE_SAMPLESHEET(CONVERT_EXCEL.out[0])
            ch_versions = ch_versions.mix(CONVERT_EXCEL.out.versions.first(), VALIDATE_SAMPLESHEET.out.versions.first())
        
        } else if (params.samplesheet.endsWith('.csv')) {
            ch_samplesheet_validated = VALIDATE_SAMPLESHEET(ch_samplesheet)
            ch_versions = ch_versions.mix(VALIDATE_SAMPLESHEET.out.versions.first())
        
        } else {
            error "Unsupported samplesheet format. Use .xlsx or .csv"
        }
        
        DORADO_BASECALL_BARCODING(ch_pod5)  
        ch_versions = ch_versions.mix(DORADO_BASECALL_BARCODING.out.versions.first())

        ch_samplesheet_validated[0]
        .splitCsv(header:true)
        .filter{ it -> it.barcode =~ /^barcode*/ }
        .map { row -> tuple( row.sample, row.barcode ) }
        .combine( DORADO_BASECALL_BARCODING.out.ch_bam_pass )
        | MERGE_READS 
        ch_versions = ch_versions.mix(MERGE_READS.out.versions.first())

    } else {
        DORADO_BASECALL(ch_pod5)
        ch_versions = ch_versions.mix(DORADO_BASECALL.out.versions.first())
    }

    emit:
    ch_basecall = params.kit ? MERGE_READS.out[0] : DORADO_BASECALL.out[0]
    ch_versions = ch_versions
}

workflow {
    ch_versions = Channel.empty()
    if (params.reads) {        
        if ( file(params.reads).isDirectory() ) {
            pattern = "*.{bam,fasta,fastq,fastq.gz,fq,fq.gz}"
            ch_reads = Channel.fromPath(params.reads + "/" + pattern, type: 'file', checkIfExists: true)           
        } else {
            ch_reads = Channel.fromPath(params.reads, checkIfExists: true)        
        }
    } else {
        def bc_out = basecall()
        ch_reads = bc_out.ch_basecall
        ch_versions = ch_versions.mix(bc_out.ch_versions)
    }
    
    
    CONVERT_READS(ch_reads) 
    READ_STATS(CONVERT_READS.out[0])
    ch_reads_conv = CONVERT_READS.out[0]
    ch_versions = ch_versions.mix(CONVERT_READS.out.versions.first(), READ_STATS.out.versions.first())
    
    if (params.subsample) {
        SUBSAMPLE_READS(ch_reads_conv)
        ch_reads_conv = SUBSAMPLE_READS.out.reads
        ch_versions = ch_versions.mix(SUBSAMPLE_READS.out.versions.first())
    }
    
    SAVONT_ASV(ch_reads_conv)
    SAVONT_CLASSIFY(SAVONT_ASV.out.ch_asvs)
    SAVONT_COMBINE(SAVONT_CLASSIFY.out.ch_abundance.collect())
    TAXONKIT_LINEAGE(SAVONT_CLASSIFY.out.ch_abundance)
    ch_versions = ch_versions.mix(SAVONT_ASV.out.versions.first(), SAVONT_CLASSIFY.out.versions.first(), SAVONT_COMBINE.out.versions)
    
    MAKE_REPORT(
        READ_STATS.out.stats.collect(),
        SAVONT_ASV.out.counts.collect(),
        SAVONT_CLASSIFY.out.ch_abundance.collect().ifEmpty([]),
        SAVONT_COMBINE.out.ch_combined
    )

    VERSIONS(ch_versions.collect(), ch_summary_file)
}
