#!/usr/bin/env nextflow
/*
========================================================================================
                                PHW QC Pipeline
========================================================================================
----------------------------------------------------------------------------------------
*/



def helpMessage() {
    log.info"""

    \nUsage:
    nextflow run qc.nf --fq ['*_R{1,2}.fastq.gz'] --outdir [output_directory]

    Mandatory arguments:
      --fq                          Path to paired-end input reads, with fileglob (see: https://www.nextflow.io/docs/latest/channel.html#fromfilepairs)

    Other options:
      --outdir                      The output directory where the results will be saved
                                    Default: "${params.outdir}"
      
      --name                        Assign a name to this run. Nextflow will make one up if not supplied (this run is: "${workflow.runName}")

      --db_uri                      Centrifuge database URI pointing to .tar.gz containing Centrifuge DB files (*.{1,2,3}.cf). 
                                    Can be "http://", "https://" or "ftp://". 
                                    Default: "${params.db_uri}"

      --hg_uri                      URI pointing to human genome reference fasta. Optionally gzipped. 
                                    Can be "http://", "https://" or "ftp://". 
                                    Default: "${params.hg_uri}"
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
    exit 1
}

// OUTPUT DIRECTORY
if (params.outdir) {
    outputDir = [ "${params.outdir}".replace(/\/$/, ""), custom_runName ].join('/')
} else {
    println("Please supply an output directory\n")
    helpMessage()
    exit 1
}


// INPUT CENTRIFUGEDB URI
centrifugeDbUri = params.db_uri

// INPUT HUMAN GENOME URI
humanGenomeFastaUri = params.hg_uri


// INPUT CHANNELS
// Fastq files
Channel.fromFilePairs( "${fastqGlob}" , flat: true)
       .set{ ch_inputReads }

// Centrifuge DB (set local or remote).
if ( centrifugeDbUri.startsWith("file://") ) {
    remoteCentrifugeDb = false
    Channel.fromPath( "${centrifugeDbUri}".replace("file:\\/\\/", "") )
           .set{ ch_centrifugeDbUri }
} else if ( centrifugeDbUri.startsWith("http://") || centrifugeDbUri.startsWith("http://") || centrifugeDbUri.startsWith("ftp://") ) {
    remoteCentrifugeDb = true
    Channel.from( "${ centrifugeDbUri }" )
           .set{ ch_centrifugeDbUri }
} else {
    println("Don't understand whether ${params.db_uri} is local or remote")
    helpMessage()
    exit 1
}

// Human genome fasta (set local or remote).
if ( humanGenomeFastaUri.startsWith("file://") ) {
    remoteHumanGenomeFasta = false
    Channel.fromPath( "${humanGenomeFastaUri}".replace("file:\\/\\/", "") )
           .set{ ch_humanGenomeUri }
} else if ( humanGenomeFastaUri.startsWith("http://") || humanGenomeFastaUri.startsWith("http://") || humanGenomeFastaUri.startsWith("ftp://") ) {
    remoteHumanGenomeFasta = true
    Channel.from( "${ humanGenomeFastaUri }" )
           .set{ ch_humanGenomeUri }
} else {
    println("Don't understand whether ${params.hg_uri} is local or remote")
    helpMessage()
    exit 1
}

// Dummy value channel to start krona taxonomy process
Channel.from( "taxonomy" )
       .set{ ch_kronaDummy }


if (remoteHumanGenomeFasta) {
    process PREPAREDBREMOTE_MINIMAP {

        tag { custom_runName }
    
        cpus 4

        input:
        val db_uri from ch_humanGenomeUri

        output:
        file("reference.idx") into ch_prepareDb_humanDepletion

        script:
        if ( db_uri.endsWith(".gz") )
            """
            curl -fsSL "${db_uri}" | gzip -d > reference.fna
            minimap2 -t ${task.cpus} -d reference.idx reference.fna
            """
        else 
            """
            curl -fsSL "${db_uri}"  > reference.fna
            minimap2 -t ${task.cpus} -d reference.idx reference.fna
            """
    }
} else {
    process PREPAREDBLOCAL_MINIMAP {

        tag { custom_runName }
    
        cpus 4

        input:
        file(db_uri) from ch_humanGenomeUri

        output:
        file("reference.idx") into ch_prepareDb_humanDepletion

        script:
        if ( db_uri.endsWith(".gz") )
            """
            zcat ${db_uri} > reference.fna
            minimap2 -t ${task.cpus} -d reference.idx reference.fna
            """
        else 
            """
            minimap2 -t ${task.cpus} -d reference.idx ${db_uri}
            """
    }
}


process DEPLETEHUMAN_MINIMAP {

    tag { custom_runName }
    
    cpus 8

    input:
    set sample_id, file(forward), file(reverse), file(hg_reference) from ch_inputReads.combine(ch_prepareDb_humanDepletion)

    output:
    set sample_id, file("${sample_id}.R1.fq.gz"), file("${sample_id}.R2.fq.gz") into ch_depleteHuman_trimReads

    script:
    """
    minimap2 -t ${task.cpus} -ax sr ${hg_reference} ${forward} ${reverse} | samtools view -b -f 13 - | samtools fastq -N -1 ${sample_id}.R1.fq.gz -2 ${sample_id}.R2.fq.gz -
    """
}

process TRIMREADS_TRIMGALORE {
    tag { sample_id }

    publishDir "${outputDir}/trimmed_reads", pattern: '*_val_{1,2}.fq.gz', mode: 'copy'
    publishDir "${outputDir}/fastqc", pattern: '*_fastqc.{zip,html}', mode: 'copy'
    publishDir "${outputDir}/trim_galore" , pattern: '*_trimming_report.txt', mode: 'copy'

    cpus 2

    input: 
    set sample_id, file(forward), file(reverse) from ch_depleteHuman_trimReads
 
    output:
    set sample_id, file("*_val_1.fq.gz"), file("*_val_2.fq.gz") optional true into ch_trimReads_insertSize, ch_trimReads_sampleComposition
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



if (remoteCentrifugeDb) {
    process PREPAREDBREMOTE_CENTRIFUGE {

        tag { custom_runName }
    
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
} else {
    process PREPAREDLOCAL_CENTRIFUGE {

        tag { custom_runName }
    
        cpus 1

        input:
        file(db_archive) from ch_centrifugeDbUri

        output:
        file("*.cf") into ch_prepareDb_sampleComposition

        script:
        """
        tar -xzf ${db_archive}
        find . -name "*.cf" -exec mv {} . \\;
        """
    }
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

process PREPAREDB_KRONA {

    tag { custom_runName }

    cpus 1

    input:
    val taxonomy_dir from ch_kronaDummy

    output:
    file("${taxonomy_dir}") into ch_prepareDb_sampleCompSummary

    script:
    """
    ktUpdateTaxonomy.sh ${taxonomy_dir} 
    """
}

process SAMPLECOMPSUMMARY_KRONA {

    tag { custom_runName }

    publishDir "${outputDir}/", mode: 'copy'

    input:
    file("centrifuge_reports/*") from ch_sampleComposition_compositionSummary.collect()
    file(krona_taxonomy) from ch_prepareDb_sampleCompSummary

    output:
    file "${custom_runName}_krona.html"

    script:
    """
    ktImportTaxonomy -tax ${krona_taxonomy}  -o ${custom_runName}_krona.html -m 5 -s 7 centrifuge_reports/*
    """
}
