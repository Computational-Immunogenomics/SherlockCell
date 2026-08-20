#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Include modules
include { count_cells } from './modules/count_cells'
include { SCF } from './modules/SCF'
include { create_swiftCNV_annots } from './modules/create_swiftCNV_annots'
include { swiftCNV } from './modules/swiftCNV'
include { malignant_classif } from './modules/malignant_classif'

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
        
            return tuple(row.dataset, adata_path, out_dir, cell_origin, sample_key, cell_type_key, sample_type_key) 
        }

    def geneAnnots = params.swiftCNV?.gene_annots
    if( !geneAnnots )
        error("Missing required parameter: params.swiftCNV.gene_annots. Please specify a GTF file to use.")

    gene_annots_file = file(geneAnnots)

    ch_count_input = ch_data_dirs.map { dataset, adata_path, out_dir, cell_origin, sample_key, cell_type_key, sample_type_key -> 
        tuple(dataset, adata_path, out_dir) }
    
    count_cells(ch_count_input)

    ch_SCF = ch_count_input
        .join(count_cells.out.map { dataset, count -> tuple(dataset, count.trim().toInteger()) })

    SCF(ch_SCF)

    ch_swiftCNV_annots = ch_data_dirs.map {dataset, adata_path, out_dir, cell_origin, sample_key, cell_type_key, sample_type_key ->
        tuple(dataset, out_dir, cell_origin, sample_key, cell_type_key, sample_type_key)}
        .join(count_cells.out.map { dataset, count -> tuple(dataset, count.trim().toInteger()) })
        .join(SCF.out.anndata)
        .join(SCF.out.scf_predictions)

    create_swiftCNV_annots(ch_swiftCNV_annots)


    ch_swiftCNV_input = ch_data_dirs
        .join(create_swiftCNV_annots.out.cell_annots)
        .join(count_cells.out.map { dataset, count -> tuple(dataset, count.trim().toInteger()) })


    swiftCNV(
        ch_swiftCNV_input,
        gene_annots_file, 
        params.swiftCNV.HMM,
        params.swiftCNV.plot_cnv,
        params.swiftCNV.sex_chr,
        params.cutoff
    )

    ch_malig_input = ch_data_dirs.map {dataset, adata_path, out_dir, cell_origin, sample_key, cell_type_key, sample_type_key ->
        tuple(dataset, out_dir, cell_origin, sample_key, cell_type_key, sample_type_key)}
        .join(SCF.out.anndata)
        .join(swiftCNV.out.cnv_scores)
        .join(swiftCNV.out.gene_order_swiftCNV)
        .join(swiftCNV.out.cell_annots_swiftCNV)
        .join(count_cells.out.map { dataset, count -> tuple(dataset, count.trim().toInteger()) })

    
    malignant_classif(ch_malig_input)


}