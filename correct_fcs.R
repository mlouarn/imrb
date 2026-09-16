library(flowCore)

input_dir = "Files/LiveMono_md0.1_k30_out/" # with or without / at the end
output_dir = "FCS_corrected/" # with / at the end

list_filepath = list.files(input_dir, pattern = "fcs$", full.names = T)

fix_fcs = function(in_filepath, output_dir){
  cyto = read.FCS(in_filepath)
  out_filepath = paste0(output_dir, basename(in_filepath))
  write.FCS(cyto, out_filepath)
}

lapply(list_filepath, fix_fcs, output_dir)
