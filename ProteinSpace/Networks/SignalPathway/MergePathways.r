MergePathways <- function(file_dir, pattern = NULL) {
  # This function combine pathway files in a certain directory
  
  ## === input ===
  ### file_dir: pathway file directory (subnetworks)
  ### pattern: specific types of edges needed in the merged pathway; default as all edges
  
  ## == output ===
  ### merged pathway file
  
  # load files
  file_names = dir(file_dir,full.names = T)
  all_file = c()
  for (i in 1:length(file_names)) {
    tmp.file = read.delim(file_names[i],header = T,sep = ",")
    all_file = rbind(all_file, tmp.file[,1:6])
  }
  # generate adjacency matrix and get rid of the redundant linkages
  new_file = c()
  multi_linkage = c()
  unique_nodes = unique(c(as.character(as.matrix(all_file$node1)), as.character(as.matrix(all_file$node2))))
  adj.matrix = matrix(0, nrow = length(unique_nodes), ncol = length(unique_nodes),
                      dimnames = list(unique_nodes, unique_nodes))
  node1.list = as.character(as.matrix(all_file$node1))
  node2.list = as.character(as.matrix(all_file$node2))
  for (i in 1:dim(all_file)[1]) {
    current_state = adj.matrix[node1.list[i],node2.list[i]]
    if (current_state == 0) {
      adj.matrix[node1.list[i],node2.list[i]] = as.numeric(all_file$Subtype1[i])
      new_file = rbind(new_file, all_file[i,])
    } else {
      new_state = as.numeric(all_file$Subtype1[i])
      if (new_state != current_state) {
        if (is.null(dim(multi_linkage))) {
          multi_linkage = rbind(multi_linkage, c(node1.list[i], node2.list[i], new_state))
          new_file = rbind(new_file, all_file[i,])
        } else {
          flag = F
          for (j in 1:dim(multi_linkage)[1]) {
            if (multi_linkage[j,1] == node1.list[i] & node2.list[i] == multi_linkage[j,2] & new_state == multi_linkage[j,3]) {
              flag = T
              break
            }
          }
          if (!flag) {
            multi_linkage = rbind(multi_linkage, c(node1.list[i], node2.list[i], new_state))
            new_file = rbind(new_file, all_file[i,])
          }
        }
      }
    }
  }
  list(new_file,multi_linkage)
}


    

