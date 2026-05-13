/*
=================================================================================================================
    Input Check Subworkflow Modules
=================================================================================================================
*/

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'
include { LANE_MERGE        } from '../../modules/local/lane_merge'

/*
========================================================================================================
    Run Input Check Subworkflow
========================================================================================================
*/

workflow INPUT_CHECK {
    take:
    samplesheet

    main:
    ch_versions = channel.empty()

    channel.fromPath(samplesheet)
        .splitCsv( header:false, sep:',', skip:1 )
        .map { row -> stage_fastq(row) }
        .set{ precheck_reads }

    LANE_MERGE(precheck_reads)

    emit:
    reads    =   LANE_MERGE.out.reads
    versions =   ch_versions
}

def stage_fastq(ArrayList row) {
    def meta        = [:]
    meta.id         = row[0]
    meta.single_end = false
    def array       = []
    def filesarray  = []

    row.drop(1).eachWithIndex { value, idx ->
        def i = idx + 1
        if(row[i] == "")
        {
        } else if (!file(row[i]).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read $i FastQ file does not exist!\n${row[i]}"
        } else
        {
            filesarray.add(file(row[i]))
        }
    }

    if(filesarray.size() == 1)
    {
        meta.single_end = true
    } else if( (filesarray.size() % 2) != 0)
    {
        exit 1, "ERROR: Please check input samplesheet -> Number of samples is not an even number or 1.\n$row"
    }

    array = [ meta, filesarray]
    return array
}
