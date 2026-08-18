process count_cells {
    tag "$dataset"
    label "count_cells"

    input:
    tuple val(dataset), path(adata), val(out_dir)

    output:
    tuple val(dataset), stdout 

    script:
    """
    #!/usr/bin/env python3
    import sys
    import scanpy as sc

    try:
        adata = sc.read_h5ad('${adata}')
   
    except Exception as e:
        print(f"Error: ${dataset} is not a valid h5ad file! Details: {e}", file=sys.stderr)
        sys.exit(1)

    ncells = adata.n_obs
    print(ncells, end='')

    """
}