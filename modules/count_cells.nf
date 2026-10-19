process COUNT_CELLS {
    tag "${meta.dataset}"
    label "count_cells"

    input:
    tuple val(meta), path(adata)

    output:
    tuple val(meta), stdout

    script:
    """
    #!/usr/bin/env python3
    import sys
    import scanpy as sc

    try:
        adata = sc.read_h5ad('${adata}', backed='r')
   
    except Exception as e:
        print(f"Error: ${meta.dataset} is not a valid h5ad file! Details: {e}", file=sys.stderr)
        sys.exit(1)

    ncells = adata.n_obs
    print(ncells, end='')

    """
}