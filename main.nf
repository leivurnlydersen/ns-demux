#!/usr/bin/env nextflow

params.rundir = null
params.outdir = null
params.samplesheet = null
params.whitelist = null
params.help = false

helpMessage = """
    Input:
    rundir:         Path to FASTQ run directory.
    outdir:         Where to store output data.
    samplesheet:    Path to sample sheet.
    """.stripIndent()

if (params.help) {
    log.info helpMessage
    exit 0
}

assert params.rundir != null, 'Input parameter "rundir" cannot be unassigned.'
assert params.outdir != null, 'Input parameter "outdir" cannot be unassigned.'
assert params.samplesheet != null, 'Input parameter "samplesheet" cannot be unassigned.'
assert params.whitelist != null, 'Input parameter "whitelist" cannot be unassigned.'

rundir = file(params.rundir)
outdir = file(params.outdir)
samplesheet = file(params.samplesheet)
whitelist = file(params.whitelist)

println "D E M U X    L I N K S E Q   "
println "================================="
println "rundir              : ${rundir}"
println "outdir              : ${outdir}"
println "samplesheet         : ${samplesheet}"
println "================================="
println "Command line        : ${workflow.commandLine}"
println "Profile             : ${workflow.profile}"
println "Project dir         : ${workflow.projectDir}"
println "Launch dir          : ${workflow.launchDir}"
println "Work dir            : ${workflow.workDir}"
println "Container engine    : ${workflow.containerEngine}"
println "================================="
println "Project             : $workflow.projectDir"
println "Git info            : $workflow.repository - $workflow.revision [$workflow.commitId]"
println "Cmd line            : $workflow.commandLine"
println "Manifest version    : $workflow.manifest.version"
println "================================="


// Call bcl2fastq, performing simultaneous basecalling and demultiplexing.
// --use-bases-mask will use RunInfo.xml (in the run directory) to determine the length of read 1 and 2
// and of the index.
// Adapter sequences (read 1 and read2) should be contained in the sample sheet.
process bcl2fastq {
    publishDir "$outdir", mode: 'copy', pattern: '.command.log', saveAs: {filename -> 'bcl2fastq/log.log'}
    publishDir "$outdir", mode: 'copy', pattern: 'outs/Stats/Stats.json', saveAs: {filename -> 'bcl2fastq/Stats.json'}

    output:
    file "outs/*fastq.gz" into fastq_trim_adapters_ch
    file '.command.log'
    file 'outs/Stats/Stats.json'

    script:
    if(task.cpus > 20) {
        p_threads = task.cpus - 8
        w_threads = 4
        r_threads = 4
    } else if(task.cpus > 10) {
        p_threads = task.cpus - 2
        w_threads = 1
        r_threads = 1
    } else if(task.cpus > 3) {
        p_threads = task.cpus - 2
        w_threads = 1
        r_threads = 1
    } else {
        p_threads = 1
        w_threads = 1
        r_threads = 1
    }
    """
    bcl2fastq \
        -R $rundir \
        -o outs \
        --interop-dir interop \
        --sample-sheet $samplesheet \
        --use-bases-mask Y*,I*,Y* \
        --minimum-trimmed-read-length 8 \
        --mask-short-adapter-reads 8 \
        --ignore-missing-positions \
        --ignore-missing-filter \
        --ignore-missing-bcls \
        -p $p_threads -r $r_threads -w $w_threads
    """
}

// Prepare the data channel for the trim_adapters process. We require (sample, lane, FASTQ read 1, FASTQ read 2) tuples.

// Get (key, FASTQ files) tuples, where key is a (sample, lane) tuple.
// First flatten the channel because each instance of process "bcl2fastq" outputs a tuple.
// Then map the channel to (key, FASTQ path) tuples. This channel has one record per file.
// Then group the records by key to obtain records with (key, FASTQ list) tuples, where
// "key" is unchanged and "FASTQ list" is a list of FASTQ files (read 1 and 2, specifically).
fastq_trim_adapters_ch.flatten()
    .map { file ->
        // Note that it is not necessary to extract all these fields, it's done for clarity.
        def sample = file.name.toString().split('_')[0]
        def number = file.name.toString().split('_')[1]
        def lane = file.name.toString().split('_')[2]
        def read = file.name.toString().split('_')[3]
        def key = tuple(sample, lane)
        return tuple(key, file)}.groupTuple().set { fastq_trim_adapters_ch }

// Convert the channel records from (key, FASTQ list) tuples to (key, read 1, read 2) tuples.
fastq_trim_adapters_ch.map { it ->
    key = it[0]
    fastq1 = it[1][0]
    fastq2 = it[1][1]
    // Get the "R1" or "R2" string from the first FASTQ.
    fastq1_read = fastq1.baseName.toString().split('_')[3]
    // Assign read 1 and 2 to fastq 1 and 2.
    if(fastq1_read == "R1") {
        read1 = fastq1
        read2 = fastq2
    } else {
        read1 = fastq2
        read2 = fastq1
    }
    return tuple(key, read1, read2)}.into {
        fastq_trim_adapters_ch; temp_ch}

// Parses samplesheet and saves adapter in a FASTA file.
process extract_adapter {
    label 'small_mem'

    output:
    file "adapter.fasta" into adapter_fasta_ch

    script:
    """
    samplesheet_extract_adapter.py $samplesheet > adapter.fasta
    """
}

// Combine the FASTQs with the adapter FASTA file to get (key, read1 FASTQ, read2 FASTQ, adapter FASTA) tuples.
trim_adapters_data_ch = fastq_trim_adapters_ch.combine(adapter_fasta_ch)

// Trim adapters.
process trim_adapters {
    label 'small_mem'

    publishDir "$outdir/$sample/logs/adapter_trim", mode: 'copy', pattern: 'bbduk.log', saveAs: { filename -> "${lane}.log" }

    input:
    set key, file(read1), file(read2), file(adapter_fasta) from trim_adapters_data_ch

    output:
    set key, file("*R1*adapter_trimmed.fastq.gz"), file("*R2*adapter_trimmed.fastq.gz") into fastq_bctrim_ch
    file 'bbduk.log'

    when:
    key[0] != "Undetermined"

    script:
    sample = key[0]
    lane = key[1]
    """
    # Trim adapters from 3' end (ktrim=r) with up to 2 mismatches (hdist=2).
    # k-mer size 21, and 11 at the end of the read (mink=11).
    # Use pair overlap detection (tbo), and trim both reads to the same length (tpe).
    bbduk.sh -Xmx${task.memory.toGiga()}g in1=$read1 in2=$read2 out1=$sample\\_$lane\\_R1\\_adapter_trimmed.fastq.gz out2=$sample\\_$lane\\_R2\\_adapter_trimmed.fastq.gz ref=$adapter_fasta ktrim=r k=21 mink=11 hdist=2 tbo tpe 2> bbduk.log
    """
}

// Trim 10x barcode from read 2.
// The barcode is taken from the first 16 bases of read 1.
// If the barcode does not match any in the list of known barcodes (whitelist), we do not trim.
// FIXME: trimR2bc.py reads FASTQ records in chunks. If these chunks are large, this process may require a non-trivial amount of memory.
process bctrim {
    label 'small_mem'

    publishDir "$outdir/$sample/logs/bctrim", mode: 'copy', pattern: 'bctrim_stats.log', saveAs: { filename -> "${lane}.log" }

    input:
    set key, file(read1), file(read2) from fastq_bctrim_ch

    output:
    set key, file("*R1*bctrimmed.fastq.gz"), file("*R2*bctrimmed.fastq.gz") into fastq_polyg_trim_ch, fastq_check_bctrim_sync_ch
    file 'bctrim_stats.log'

    script:
    sample = key[0]
    lane = key[1]
    """
    trimR2bc.py $read1 $read2 $whitelist $sample\\_$lane\\_R2\\_bctrimmed.fastq 1> bctrim_stats.log
    gzip -k $sample\\_$lane\\_R2\\_bctrimmed.fastq
    # Even though we did not change R1, we rename it before outputting it, as having differently named R1 and R2 can cause confusion.
    cp $read1 $sample\\_$lane\\_R1\\_bctrimmed.fastq.gz
    """
}

//// With this process I can unsynchronize the reads to check if the "check_sync" process works.
//process mess_up_sync_test {
//    input:
//    set key, file(fastqs) from fastq_check_sync_ch
//
//    output:
//    set key, file("*unsynced.fastq.gz"), file(read2) into fastq_check_sync_unsynced_ch
//
//    script:
//    sample = key[0]
//    lane = key[1]
//    read1 = fastqs[0]
//    read2 = fastqs[1]
//    """
//    # Mess up read 1.
//    # Take the first read and place it at the end of the file.
//    zcat $read1 > temp.fastq
//    tail -n 4 temp.fastq > $sample\\_$lane\\_R1\\_unsynced.fastq
//    head -n -4 temp.fastq >> $sample\\_$lane\\_R1\\_unsynced.fastq
//    gzip -k $sample\\_$lane\\_R1\\_unsynced.fastq
//    # Can use this if I want to remove lines:
//    #zcat $read1 | head -n -3 > $sample\\_$lane\\_R1\\_unsynced.fastq
//    """
//}

// Check that the read 1 and 2 are synchronized. If they are not, this process will throw an error
// and the pipeline will exit.
process check_bctrim_sync {
    label 'small_mem'

    input:
    set key, file(read1), file(read2) from fastq_check_bctrim_sync_ch

    script:
    """
    # Check if reads are synchronized.
    reformat.sh -Xmx${task.memory.toGiga()}g in=$read1 in2=$read2 vpair
    """
}

// Trim poly-G tail.
// FIXME: how much memory does fastp actually need?
process polyG_trim {
    publishDir "$outdir/$sample/logs/polyG_trim", mode: 'copy', pattern: 'polyG_trim.log', saveAs: { filename -> "${lane}.log" }

    input:
    set key, file(read1), file(read2) from fastq_polyg_trim_ch

    output:
    set key, file("*R1*polyGtrimmed.fastq.gz"), file("*R2*polyGtrimmed.fastq.gz") into fastq_qtrim_r1_ch, fastq_qtrim_r2_ch
    file 'polyG_trim.log'

    script:
    sample = key[0]
    lane = key[1]
    """
    # Trim poly G of reads
    # Q: Disable quality filter, -L: Disable length filter, -A: Disable adapter filtering, -g: Enable polyG trim
    fastp -i $read1 -I $read2 -o $sample\\_$lane\\_R1\\_polyGtrimmed.fastq.gz -O $sample\\_$lane\\_R2\\_polyGtrimmed.fastq.gz -Q -L -A -g -h "polyG_trim_log.html" -j "polyG_trim_log.json" 2> polyG_trim.log
    # NOTE: is does not seem like HTML and JSON reports (-h and -j) can be disabled.
    """
}

// Do quality trimming and minimum length filtering (read 1).
// FIXME: how much memory does sickle actually need?
process quality_trim_read1 {
    publishDir "$outdir/$sample/logs/quality_trim", mode: 'copy', pattern: 'sickle.log', saveAs: { filename -> "${lane}_R1.log" }

    input:
    set key, file(read1), file(read2) from fastq_qtrim_r1_ch

    output:
    set key, file("*R1*qtrimmed.fastq.gz") into fastq_qtrimmed_r1_ch
    file 'sickle.log'

    script:
    sample = key[0]
    lane = key[1]
    """
    # We don't trim from 5' end, because this would trim the barcode (-x).
    sickle se -f $read1 -t sanger -g -o $sample\\_$lane\\_R1\\_qtrimmed.fastq.gz -x -q 20 -l 58 1> sickle.log
    """
}

// Do quality trimming and minimum length filtering (read 2).
process quality_trim_read2 {
    publishDir "$outdir/$sample/logs/quality_trim", mode: 'copy', pattern: 'sickle.log', saveAs: { filename -> "${lane}_R2.log" }

    input:
    set key, file(read1), file(read2) from fastq_qtrim_r2_ch

    output:
    set key, file("*R2*qtrimmed.fastq.gz") into fastq_qtrimmed_r2_ch
    file 'sickle.log'

    script:
    sample = key[0]
    lane = key[1]
    """
    sickle se -f $read2 -t sanger -g -o $sample\\_$lane\\_R2\\_qtrimmed.fastq.gz -q 20 -l 35 1> sickle.log
    """
}

// Join the (key, read1) with the (key, read2) channels to obtain a (key, read1, read2) channel.
fastq_qtrimmed_ch = fastq_qtrimmed_r1_ch.join(fastq_qtrimmed_r2_ch)

// Synchronize reads, if the reads got out of order.
process sync_qtrim_reads {
    publishDir "$outdir/$sample/fastqs", mode: 'copy', pattern: "*.fastq.gz"
    publishDir "$outdir/$sample/logs/sync_reads", mode: 'copy', pattern: '.command.log', saveAs: { filename -> "${lane}.log" }

    input:
    set key, file(read1), file(read2) from fastq_qtrimmed_ch

    output:
    set key, file("*R1*.fastq.gz"), file("*R2*.fastq.gz") into fastq_qc_ch
    file '.command.log'

    script:
    sample = key[0]
    lane = key[1]
    """
    # We use ziplevel=1 to get fast but low-level compression.
    # NOTE: singletons.fastq.gz should be empty.
    repair.sh -Xmx${task.memory.toGiga()}g ziplevel=1 in1=$read1 in2=$read2 out1=$sample\\_$lane\\_R1.fastq.gz out2=$sample\\_$lane\\_R2.fastq.gz outs=singletons.fastq.gz repair
    """
}

// Run FastQC for QC metrics of raw data.
// NOTE: FastQC claims 250 MB of memory for every thread that is allocated to it.
process fastqc_analysis {
    memory { 250.MB * task.cpus }

    publishDir "$outdir/$sample/fastqc", mode: 'copy', pattern: '{*.zip,*.html}',
        saveAs: {filename -> filename.indexOf('.zip') > 0 ? "zips/$filename" : "$filename"}
    publishDir "$outdir/$sample/fastqc", mode: 'copy', pattern: '.command.log',
        saveAs: {filename -> 'fastqc.log'}

    input:
    set key, file(read1), file(read2) from fastq_qc_ch

    output:
    set sample, file('*.{zip,html}') into fastqc_report_ch
    set sample, file('.command.log') into fastqc_stdout_ch
    val 'done' into status_ch  // Tells the MultiQC process that the pipeline is finished.

    script:
    //fastq_list = (fastqs as List).join(' ')
    sample = key[0]
    lane = key[1] // Unused.
    """
    # We unset the DISPLAY variable to avoid having FastQC try to open the GUI.
    unset DISPLAY
    mkdir tmp
    fastqc -q --dir tmp --threads ${task.cpus} --outdir . $read1 $read2
    """
}

process multiqc {
    publishDir "$outdir/multiqc", mode: 'copy', overwrite: true

    input:
    val status from status_ch.collect()

    output:
    file "multiqc_report.html" into multiqc_report_ch
    file "multiqc_data" into multiqc_data_ch

    script:
    """
    multiqc -f $outdir
    """
}

workflow.onComplete {
    log.info "D E M U X    L I N K S E Q   "
    log.info "================================="
    log.info "rundir              : ${rundir}"
    log.info "outdir              : ${outdir}"
    log.info "samplesheet         : ${samplesheet}"
    log.info "================================="
    log.info "Command line        : ${workflow.commandLine}"
    log.info "Profile             : ${workflow.profile}"
    log.info "Project dir         : ${workflow.projectDir}"
    log.info "Launch dir          : ${workflow.launchDir}"
    log.info "Work dir            : ${workflow.workDir}"
    log.info "Container engine    : ${workflow.containerEngine}"
    log.info "================================="
    log.info "Project             : $workflow.projectDir"
    log.info "Git info            : $workflow.repository - $workflow.revision [$workflow.commitId]"
    log.info "Cmd line            : $workflow.commandLine"
    log.info "Manifest version    : $workflow.manifest.version"
    log.info "================================="
    log.info "Completed at        : ${workflow.complete}"
    log.info "Duration            : ${workflow.duration}"
    log.info "Success             : ${workflow.success}"
    log.info "Exit status         : ${workflow.exitStatus}"
    log.info "Error report        : ${(workflow.errorReport ?: '-')}"
    log.info "================================="
}

