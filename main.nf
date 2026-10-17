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

    ch_data_dirs = channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true, sep: '\t')
        .map { row -> 
            def adata_path = file(row.adata_path, type: 'file')
            def base_path = row.outdir ?: adata_path.parent
            def out_dir   = "${base_path}/reannot_results"
            def cell_origin = row.cell_origin 
            def sample_key    = row.sample_key ?: params.sample_key
            def cell_type_key = row.cell_type_key ?: params.cell_type_key
            def sample_type_key = row.sample_type_key ?: params.sample_type_key
            def embedding_key = row.embedding_key ?: params.embedding_key
            def exclude_from_reference = row.exclude_from_reference ?: null
            def counts_in_X = (row.counts_in_X != null && row.counts_in_X != '') ? row.counts_in_X.toBoolean() : params.swiftCNV.counts_in_X

            return tuple(row.dataset, adata_path, out_dir, cell_origin, sample_key, cell_type_key, sample_type_key, embedding_key, exclude_from_reference, counts_in_X) 
        }

    def geneAnnots = params.swiftCNV?.gene_annots
    if( !geneAnnots )
        error("Missing required parameter: params.swiftCNV.gene_annots. Please specify a GTF file to use.")

    gene_annots_file = file(geneAnnots)

    ch_count_input = ch_data_dirs.map { dataset, adata_path, out_dir, cell_origin, sample_key, cell_type_key, sample_type_key, embedding_key, exclude_from_reference, counts_in_X -> 
        tuple(dataset, adata_path, out_dir) }
    
    COUNT_CELLS(ch_count_input)

    ch_SCF = ch_count_input
        .join(COUNT_CELLS.out.map { dataset, count -> tuple(dataset, count.trim().toInteger()) })

    SCF(ch_SCF)

    ch_swiftCNV_annots = ch_data_dirs.map {dataset, adata_path, out_dir, cell_origin, sample_key, cell_type_key, sample_type_key, embedding_key, exclude_from_reference, counts_in_X ->
        tuple(dataset, out_dir, cell_origin, sample_key, cell_type_key, sample_type_key, exclude_from_reference)}
        .join(COUNT_CELLS.out.map { dataset, count -> tuple(dataset, count.trim().toInteger()) })
        .join(SCF.out.anndata)
        .join(SCF.out.scf_predictions)

    SWIFTCNV_ANNOTS(ch_swiftCNV_annots)


    ch_swiftCNV_input = ch_data_dirs
        .join(SWIFTCNV_ANNOTS.out.cell_annots)
        .join(COUNT_CELLS.out.map { dataset, count -> tuple(dataset, count.trim().toInteger()) })


    SWIFTCNV(
        ch_swiftCNV_input,
        gene_annots_file, 
        params.swiftCNV.HMM,
        params.swiftCNV.plot_cnv,
        params.swiftCNV.sex_chr,
        params.swiftCNV.cutoff
    )

    ch_malig_input = ch_data_dirs.map {dataset, adata_path, out_dir, cell_origin, sample_key, cell_type_key, sample_type_key, embedding_key, exclude_from_reference, counts_in_X ->
        tuple(dataset, out_dir, cell_origin, sample_key, cell_type_key, sample_type_key, embedding_key, exclude_from_reference)}
        .join(SCF.out.anndata)
        .join(SWIFTCNV.out.cnv_scores)
        .join(SWIFTCNV.out.gene_order_swiftCNV)
        .join(SWIFTCNV.out.cell_annots_swiftCNV)
        .join(COUNT_CELLS.out.map { dataset, count -> tuple(dataset, count.trim().toInteger()) })

    
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

    ch_out_dir = ch_data_dirs.map { dataset, adata_path, out_dir, cell_origin, sample_key, cell_type_key, sample_type_key, embedding_key, exclude_from_reference, counts_in_X -> 
            tuple(dataset, out_dir) }

    ch_logs = SCF.out.log_scf.map { dataset, log_file -> [ step: "SCF", dataset: dataset, file: log_file ] }
        .mix(
            SWIFTCNV_ANNOTS.out.log_annots.map { dataset, log_file -> [ step: "SWIFTCNV_ANNOTS", dataset: dataset, file: log_file ] },
            SWIFTCNV.out.log_swiftcnv.map { dataset, log_file -> [ step: "SWIFTCNV", dataset: dataset, file: log_file ] },
            MALIGNANT_CLASSIF.out.log_maligclf.map { dataset, log_file -> [ step: "MALIGNANT_CLASSIF", dataset: dataset, file: log_file ] }
        )
        .toList()
        .map { logs -> [ logs ] }
        

    ch_multiqc = ch_multiqc_files
                    .join(ch_out_dir)
                    .combine(ch_logs)

    REPORT(
    ch_multiqc,
    file(params.multiqc_config),
    file("${projectDir}/assets/custom.css"),
    file("${projectDir}/docs/images/SherlockCell_full.png")
    )

}