/*
====================================================================================================
    Nextclade Dataset and Analysis Subworkflow Modules
====================================================================================================
*/

include { UNTAR as UNTAR_NEXTCLADE_DB        } from '../../modules/nf-core/untar/main'
include { NEXTCLADE_DATASETGET               } from '../../modules/local/nextclade_datasetget.nf'
include { NEXTCLADE_RUN                      } from '../../modules/local/nextclade_run.nf'
include { NEXTCLADE_PARSER                   } from '../../modules/local/nextclade_parser.nf'
include { NEXTCLADE_REPORT                   } from '../../modules/local/nextclade_report.nf'
/*
====================================================================================================
    Run Nextclade Dataset and Analysis Subworkflow
====================================================================================================
*/

workflow NEXTCLADE_DATASET_AND_ANALYSIS {
    take:
    dataset
    HA

    main:
    ch_versions         = Channel.empty()
    ch_nextclade_report = Channel.empty()
    ch_aligned_fasta    = Channel.empty()

    if (params.skip_nextclade) return

    dataset_keyed = dataset.map { meta, ds_file ->
        tuple(meta.id, meta, ds_file)
    }

    ha_keyed = HA.map { meta, ha_fasta ->
        tuple(meta.id, meta, ha_fasta)
    }

    joined = dataset_keyed
        .join(ha_keyed)
        .map { id, meta1, ds_file, meta2, ha_fasta ->
            tuple(meta1, ds_file, ha_fasta)
        }

    NEXTCLADE_DATASETGET(joined.map { meta, ds_file, ha_fasta -> tuple(meta, ds_file) })
    ch_versions = ch_versions.mix(NEXTCLADE_DATASETGET.out.versions)

    dataset_2 = NEXTCLADE_DATASETGET.out.dataset_2

    dataset2_keyed = dataset_2.map { meta, ds2 -> tuple(meta.id, meta, ds2) }
    ha2_keyed      = joined.map { meta, ds_file, ha_fa -> tuple(meta.id, meta, ha_fa) }

    run_joined = dataset2_keyed
        .join(ha2_keyed)
        .map { id, meta1, ds2, meta2, ha_fa ->
            tuple(meta1, ds2, ha_fa)
        }

    NEXTCLADE_RUN(
        run_joined.map { meta, ds2, ha_fa -> tuple(meta, ds2) },
        run_joined.map { meta, ds2, ha_fa -> tuple(meta, ha_fa) }
    )

    ch_aligned_fasta  = ch_aligned_fasta.mix(NEXTCLADE_RUN.out.fasta_aligned)
    ch_nextclade_report = NEXTCLADE_RUN.out.csv

    NEXTCLADE_PARSER( NEXTCLADE_RUN.out.parser_input.filter { meta, f -> f } )
    parser_tsv_files = NEXTCLADE_PARSER.out.nextclade_parser_tsv

    parser_tsv_best = parser_tsv_files
        .groupTuple(by: 0)
        .map { meta, files ->
            def chosen = files.sort { a, b -> a.size() <=> b.size() ?: a.name <=> b.name }.last()
            tuple(meta, chosen)
        }

    ch_combined_parser_tsv_results = parser_tsv_best
        .map { meta, tsv -> tsv.text }
        .collect()
        .map { texts ->
            def header = null
            def rows = []
            texts.each { txt ->
                def lines = txt?.readLines()?.findAll { it?.trim() }
                if (!lines) return
                header = header ?: lines[0]
                rows.addAll(lines.drop(1))
            }
            if (header == null) return ""
            ([header] + rows.unique()).join("\n") + "\n"
        }

    NEXTCLADE_REPORT(ch_combined_parser_tsv_results)

    emit:
    fasta_aligned        = NEXTCLADE_RUN.out.fasta_aligned
    parser_input         = NEXTCLADE_RUN.out.parser_input
    nextclade_report     = ch_nextclade_report
    nextclade_report_tsv = NEXTCLADE_REPORT.out.nextclade_report_tsv
    versions             = ch_versions
}
