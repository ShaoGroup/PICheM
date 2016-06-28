# This function gets certain downstream/upstream genes in given kegg pathways

## === input ===
### file: pathway file name;
### id: list of kegg pathway ids (in the form of hsaXXXXXX)
### gene = NA: default as return all the genes
### direction: 1 for downstream; -1 for upstream
### step = 20: default as all the upstream/downstream genes; if specified, iterate given steps

## == output ===
### genelist: certain downstream/upstream genes in given kegg pathways

GetGeneList <- function(file, id, gene = NA, direction = 1, steps = 20) {
  edges.file = read.delim(file, header = T, sep = ",")
  filter = rep(F, dim(edges.file)[1])
  for (i in 1:length(id)) {
    filter = filter | (as.character(edges.file$Source) == id[i])
  }
  filter = filter & !(as.character(edges.file$Type) == "downstream")  # here we don't need "pathway nodes"
  filter = filter & !(as.character(edges.file$Subtype1) == "binding/association") # here we don't consider binding as downstream
  trimmed.file = edges.file[filter,]
  
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
  
  result
}


