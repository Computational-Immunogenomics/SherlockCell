#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Include modules
include { count_cells } from './modules/count_cells'
include { create_swiftCNV_annots } from './modules/create_swiftCNV_annots'
include { swiftCNV } from './modules/swiftCNV'
include { malignant_classif } from './modules/malignant_classif'

workflow {

    ch_data_dirs = channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true, sep: '\t')
        .map { row -> 
            def adata_path = file(row.adata_path, type: 'file')
            def base_path = params.outdir ?: adata_path.parent
            def out_dir   = "${base_path}/reannot_results"

            def ref_cells = row.reference_cells
            def cell_origin = row.cell_origin 
            def sample_key    = row.sample_key ?: params.sample_key
            def cell_type_key = row.cell_type_key ?: params.cell_type_key
            def annot_mode = row.annot_mode ?: params.annot_mode
        
            return tuple(row.dataset, adata_path, out_dir, ref_cells, cell_origin, sample_key, cell_type_key, annot_mode) 
        }

    def geneAnnots = params.swiftCNV?.gene_annots
    if( !geneAnnots )
        error("Missing required parameter: params.swiftCNV.gene_annots. Please specify a GTF file to use.")

    gene_annots_file = file(geneAnnots)

    ch_count_input = ch_data_dirs.map { dataset, adata_path, out_dir, ref_cells, cell_origin, sample_key, cell_type_key, annot_mode -> 
        tuple(dataset, adata_path) }
    
    count_cells(ch_count_input)

    ch_swiftCNV_annots = ch_count_input
                        .join(count_cells.out.map { dataset, count -> tuple(dataset, count.trim().toInteger()) })

    create_swiftCNV_annots(ch_swiftCNV_annots)


    ch_swiftCNV_input = ch_data_dirs
        .join(create_swiftCNV_annots.out.cell_annots)
        .join(count_cells.out.map { dataset, count -> tuple(dataset, count.trim().toInteger()) })


    swiftCNV(
        ch_swiftCNV_input,
        gene_annots_file, 
        params.swiftCNV.HMM,
        params.swiftCNV.plot_cnv
    )

    ch_malig_input = ch_data_dirs
        .join(swiftCNV.out.cnv_scores)
        .join(swiftCNV.out.gene_order_swiftCNV)
        .join(swiftCNV.out.cell_annots_swiftCNV)
        .join(count_cells.out.map { dataset, count -> tuple(dataset, count.trim().toInteger()) })

    
    malignant_classif(ch_malig_input)


}