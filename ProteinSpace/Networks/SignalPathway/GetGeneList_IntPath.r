GetGeneList_IntPath <- function(gene = NA, label = NA, direction = 1, steps = 20) {
  # This function gets certain downstream/upstream genes in given kegg pathways
  
  ## === input ===
  ### gene = NA: default as return all the genes
  ### label: list of pathway ids (some problem with the space!)
  ### direction: 1 for downstream; -1 for upstream
  ### step = 20: default as all the upstream/downstream genes; if specified, iterate given steps
  
  ## == output ===
  ### genelist: certain downstream/upstream genes in given kegg pathways
  
  load("./repertoire/IntPathRepertoire.RData")
  
  if (!is.na(label)) {
    filter = rep(F, dim(IntPathGenePairs)[1])
    for (i in 1:length(label)) {
      filter = filter | (as.character(IntPathGenePairs$label) == label[i])
    }
    filter = filter & !(as.character(IntPathGenePairs$reaction_type) == "GPrel") # here we don't consider binding as downstream
    trimmed.file = IntPathGenePairs[filter,]
  } else {
    trimmed.file = IntPathGenePairs
  }
  
  
  if (is.na(gene)) {
    result = unique(c(as.character(trimmed.file$node1), as.character(trimmed.file$node2)))
  } else {
    # generate matrix representing linkages
    node1.list = as.character(trimmed.file$node1)
    node2.list = as.character(trimmed.file$node2)
    all_nodes = unique(c(node1.list, node2.list))
    linkage.matrix = matrix(0, nrow = length(all_nodes), ncol = length(all_nodes), dimnames = list(all_nodes, all_nodes))
    for (j in 1:dim(trimmed.file)[1]) {
      linkage.matrix[node1.list[j], node2.list[j]] = 1
    }
    # get required genes
    if (direction == -1) {
      linkage.matrix = t(linkage.matrix)
    }
    result = c()
    current_number = 0
    for (cnt in 1:steps) {
      if (length(gene) > 1) {
        gene = colnames(linkage.matrix)[colSums(linkage.matrix[gene,]) > 0]
      } else {
        gene = colnames(linkage.matrix)[linkage.matrix[gene,] == 1]
      }
      result = c(result, gene)
      result = unique(result)
      tmp_cnt = length(result)
      if(current_number == tmp_cnt) {
        break
      } else {
        current_number = tmp_cnt
      }
    }
  }
  
  list(gene.list = result, steps =cnt)
}


