


# Function blockSample
blockSample <- function(mat, el, el2, name) {
  block_rand <- c()
  names <- c()
  randnames <- c()
  #### Assumes el == info[rownames(binary_matrix), "site"]
  #### Assumes el2 == info[rownames(binary_matrix), "habitat"]
  
  for(i in 1:length(unique(el))) {
    for(j in 1:length(unique(el2))) {
      names <- c(names, name[which(el %in% unique(el)[i] & el2 %in% unique(el2)[j])])## Vector of IDs belonging to site[i] and habitat[j]
      randnames <- c(randnames, sample(name[which(el %in% unique(el)[i] & el2 %in% unique(el2)[j])]))## Shuffle those IDs
    }
  }
  if(is.vector(mat)) {
    block_rand <- mat[names, ]
    names(block_rand) <- randnames
    block_rand <- block_rand[names(mat), ]
  } else {
    block_rand <- mat[names, ]## binary_matrix with IDs belonging to site[i] and habitat[j]
    rownames(block_rand) <- randnames### Set the shuffled IDs as row names of the matrix above
    block_rand <- block_rand[rownames(mat), ]# Reorder the rows of block_rand based on rownames(mat)
    # Row order of block_rand matches mat
  }
  
  return(list(matrix = block_rand, rownames = randnames))# Contents of the matrix remain unchanged, row names are the shuffled IDs
}
# Generates a randomized matrix where IDs (items other than site and habitat -> host) are shuffled within each combination of site and habitat