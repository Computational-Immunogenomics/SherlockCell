process count_cells {
    tag "$dataset"
    executor 'local' 

    input:
    tuple val(dataset), path(adata)

    output:
    tuple val(dataset), stdout 

    script:
    """
    #!/opt/env/bin/python
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