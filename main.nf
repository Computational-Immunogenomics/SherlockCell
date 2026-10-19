#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Include modules
include { COUNT_CELLS } from './modules/count_cells'
include { SCF } from './modules/SCF'
include { SWIFTCNV_ANNOTS } from './modules/create_swiftCNV_annots'
include { SWIFTCNV } from './modules/swiftCNV'
include { MALIGNANT_CLASSIF } from './modules/malignant_classif'
include { REPORT } from './modules/multiqc'

workflow {

    main:

    println "Running SherlockCell 🕵️ \n"

    ch_data_dirs = channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true, sep: '\t')
        .map { row -> 
            def adata_path = file(row.adata_path, checkIfExists: true)
            def base_path  = row.outdir ?: adata_path.parent

            def meta = [
                dataset               : row.dataset,
                out_dir               : "${base_path}/reannot_results",
                adata_path            : adata_path.toString(),
                cell_origin           : row.cell_origin,
                sample_key            : row.sample_key ?: params.sample_key,
                cell_type_key         : row.cell_type_key ?: params.cell_type_key,
                sample_type_key       : row.sample_type_key ?: params.sample_type_key,
                embedding_key         : row.embedding_key ?: params.embedding_key,
                exclude_from_reference: row.exclude_from_reference ?: null,
                counts_in_X           : (row.counts_in_X != null && row.counts_in_X != '') 
                                        ? row.counts_in_X.toBoolean() 
                                        : params.swiftCNV.counts_in_X
            ]

            return [ meta, adata_path ]

        }

    def geneAnnots = params.swiftCNV?.gene_annots
    if( !geneAnnots )
        error("Missing required parameter: params.swiftCNV.gene_annots. Please specify a GTF file to use.")

    gene_annots_file = file(geneAnnots, checkIfExists: true)
    
    COUNT_CELLS(ch_data_dirs)

    ch_cell_counts = COUNT_CELLS.out
        .map { meta, count -> [ meta, count.toString().trim().toInteger() ] }

    ch_SCF = ch_data_dirs.join(ch_cell_counts)
    SCF(ch_SCF)

    ch_swiftCNV_annots = ch_data_dirs
        .map { meta, _adata -> meta }
        .join(ch_cell_counts)
        .join(SCF.out.anndata)
        .join(SCF.out.scf_predictions)

    SWIFTCNV_ANNOTS(ch_swiftCNV_annots)


    ch_swiftCNV_input = ch_data_dirs
        .join(SWIFTCNV_ANNOTS.out.cell_annots)
        .join(ch_cell_counts)

    SWIFTCNV(
        ch_swiftCNV_input,
        gene_annots_file,
        params.swiftCNV.HMM,
        params.swiftCNV.plot_cnv,
        params.swiftCNV.sex_chr,
        params.swiftCNV.cutoff
    )

    ch_malig_input = ch_data_dirs
        .map { meta, _adata -> meta }
        .join(ch_cell_counts)
        .join(SCF.out.anndata)
        .join(SWIFTCNV.out.cnv_scores)
        .join(SWIFTCNV.out.gene_order_swiftCNV)
        .join(SWIFTCNV.out.cell_annots_swiftCNV)

    MALIGNANT_CLASSIF(ch_malig_input)

    // Group all plots across all datasets into a single collection
    ch_multiqc_files = SWIFTCNV_ANNOTS.out.query_pcts_plot
        .mix(
            MALIGNANT_CLASSIF.out.umap_classif,
            MALIGNANT_CLASSIF.out.cell_metadata,
            MALIGNANT_CLASSIF.out.cell_metadata_2,
            MALIGNANT_CLASSIF.out.cnv_summary_heatmap,
            MALIGNANT_CLASSIF.out.density_plot
        )
        .groupTuple()


    ch_logs = SCF.out.log_scf.map { meta, log_file -> [ step: "SCF", dataset: meta.id, file: log_file ] }
        .mix(
            SWIFTCNV_ANNOTS.out.log_annots.map { meta, log_file -> [ step: "SWIFTCNV_ANNOTS", dataset: meta.id, file: log_file ] },
            SWIFTCNV.out.log_swiftcnv.map { meta, log_file -> [ step: "SWIFTCNV", dataset: meta.id, file: log_file ] },
            MALIGNANT_CLASSIF.out.log_maligclf.map { meta, log_file -> [ step: "MALIGNANT_CLASSIF", dataset: meta.id, file: log_file ] }
        )
        .toList()

    ch_log_info  = ch_logs.map { list -> [ list.collect { [ step: it.step, dataset: it.dataset, filename: it.file.name ] } ] }
    ch_log_files = ch_logs.map { list -> [ list.collect { it.file } ] }
        

    ch_report_input = ch_multiqc_files
            .combine(ch_log_info)
            .combine(ch_log_files)

    REPORT(
        ch_report_input,
        file(params.multiqc_config, checkIfExists: true),
        file("${projectDir}/assets/custom.css", checkIfExists: true),
        file("${projectDir}/docs/images/SherlockCell_full.png", checkIfExists: true)
    )


    onComplete:

    if (workflow.success) {
        
        println "\n🎉 Pipeline finished successfully! \n"
    }

    onError:
    println "Error: Pipeline execution stopped with the following message: ${workflow.errorMessage}"

}