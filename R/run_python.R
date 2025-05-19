
prepare_python <- function(out_dir){

  #use_condaenv("py39", required = FALSE)
  #script_path <- "inst/run_GATMLP.py"
  reticulate::install_python(version = 3.9)
  virtualenv_create("env_SCTP")
  virtualenv_install("env_SCTP", packages = c("pytorch", "pytorch_geometric", "scikit-learn", "numpy", "scipy"))

  use_virtualenv("env_SCTP", required = FALSE)
    # The return value MUST be a pure R object, i.e., no reticulate
    # Python objects, no pointers to shared memory.


  #source_python(script_path)
  #res <- load_data(work_dir = out_dir)

  #return(list(beta=res,
  #            sc_data=sc_dataset@assays$RNA@data,
  #            st_data=st_dataset@assays$RNA@data) )
}









