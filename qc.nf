#!/usr/bin/env nextflow
/*
========================================================================================
                                PHW QC Pipeline
========================================================================================
----------------------------------------------------------------------------------------
*/



def helpMessage() {
    log.info"""

    Usage:
    nextflow run qc.nf --fq ['*_R{1,2}.fastq.gz'] --outdir [output_directory]

    Mandatory arguments:
      --fq                          Path to paired-end input reads, with fileglob (see: https://www.nextflow.io/docs/latest/channel.html#fromfilepairs)
      --outdir                      The output directory where the results will be saved
      
    Other options:
      --db_uri                      Centrifuge database URI pointing to .tar.gz containing Centrifuge DB files (*.{1,2,3}.cf). Can be "http://", "https://", "ftp://" or "file://". 
                                    Default: "${params.db_uri}"
    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}

// RUN NAME (ASSIGNED OR PROVIDED)
if (params.name) {
    custom_runName = params.name 
    } else {
    custom_runName = workflow.runName
}

// INPUT FASTQ
if (params.fq) {
    fastqGlob = params.fq
} else {
    println("Please supply a path to some fastq files\n")
    helpMessage()
    exit 0
}

// OUTPUT DIRECTORY
if (params.outdir) {
    outputDir = [ "${params.outdir}".replace(/\/$/, ""), custom_runName ].join('/')
} else {
    println("Please supply an output directory\n")
    helpMessage()
    exit 0
}


// INPUT CENTRIFUGEDB URI
centrifugeDbUri = params.db_uri



// INPUT CHANNELS
Channel.fromFilePairs( "${fastqGlob}" , flat: true)
       .set{ ch_inputReads }


Channel.from( "${ centrifugeDbUri }" )
       .set{ ch_centrifugeDbUri }



process TRIMREADS_TRIMGALORE {
    tag { sample_id }

    publishDir "${outputDir}/trimmed_reads", pattern: '*_val_{1,2}.fq.gz', mode: 'copy'
    publishDir "${outputDir}/fastqc", pattern: '*_fastqc.{zip,html}', mode: 'copy'
    publishDir "${outputDir}/trim_galore" , pattern: '*_trimming_report.txt', mode: 'copy'

    cpus 2

    input: 
    set sample_id, file(forward), file(reverse) from ch_inputReads
 
    output:
    set sample_id, file("*_val_1.fq.gz"), file("*_val_2.fq.gz") optional true into ch_trimReads_readLength, ch_trimReads_insertSize, ch_trimReads_sampleComposition
    set file("*trimming_report.txt"), file("*_fastqc.{zip,html}") optional true into ch_trimReads_qcSummary

    script:
      """
      if [[ \$(zcat ${forward} | head -n4 | wc -l) -eq 0 ]]; then
        exit 0
      else
        trim_galore --fastqc --paired ${forward} ${reverse}
      fi
      """
}

process INSERTSIZE_BBMERGE {
    tag { sample_id }

    cpus 2 

    input:
    set sample_id, file(forward), file(reverse) from ch_trimReads_insertSize

    output:
    file("${sample_id}.ihist") into ch_insertSize_qcSummary
    
    script:
    """
    bbmerge.sh strict=t reads=1000000 ihist=${sample_id}.ihist in=${forward} in2=${reverse}
    """
}

process PREPAREDB_CENTRIFUGE {
    
    cpus 1

    input:
    val db_archive_uri from ch_centrifugeDbUri

    output:
    file("*.cf") into ch_prepareDb_sampleComposition

    script:
    """
    curl -fsSL "${db_archive_uri}" | tar -xz
    find . -name "*.cf" -exec mv {} . \\;
    """
}

process SAMPLECOMPOSITION_CENTRIFUGE {
    tag { sample_id }

    publishDir "${outputDir}/centrifuge", pattern: '${sample_id}.tab', mode: 'copy'

    cpus 8

    input:
    set sample_id, file(forward), file(reverse), file(database) from ch_trimReads_sampleComposition.combine(ch_prepareDb_sampleComposition.toList())

    output:
    file("${sample_id}.tab") into ch_sampleComposition_compositionSummary

    script:
    db_prefix = database[0].toString().replace(".1.cf", "")
    
    """
    centrifuge --mm -q -p ${task.cpus} -x ${db_prefix}  -1 ${forward} -2 ${reverse} -S /dev/null --report-file ${sample_id}.tab
    """
}

process QCSUMMARY_MULTIQC {

    tag { custom_runName }

    publishDir "${outputDir}/", pattern: "${custom_runName}_multiqc*", mode: 'copy'

    input:
    file("*") from ch_insertSize_qcSummary.collect()
    file("*") from ch_trimReads_qcSummary.flatten().collect()

    output:
    file "${custom_runName}_multiqc.html"
    file "${custom_runName}_multiqc_data"

    script:
     """
     multiqc -m bbmap -m cutadapt -m fastqc -i "${custom_runName}" -n ${custom_runName}_multiqc.html .
     """
}

process SAMPLECOMPSUMMARY_KRONA {
    tag { proctag }

    publishDir "${outputDir}/", mode: 'copy'

    input:
    file("centrifuge_reports/*") from ch_sampleComposition_compositionSummary.collect()

    output:
    file "centrifuge_report.html"

    script:
    """
    ktImportTaxonomy -o ${runName}_krona.html -m 5 -s 7 centrifuge_reports/*
    """
}
