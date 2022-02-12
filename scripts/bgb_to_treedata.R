split_n_sort<-function(vector, sep="-"){
  df=data.frame(x=vector)
  colmn <- paste0("col_", 1:(max(stringr::str_count(df$x, sep)) + 1))
  tmp=sort(colmn, decreasing = T)
  tmp2 <- df %>% tidyr::separate(x, sep = sep, convert=T, into = colmn, remove = FALSE) 
  for(i in tmp){
    tmp2=tmp2[order(tmp2[,i], na.last = F),]
  }
  order_vec<-rownames(tmp2)
  return(order_vec)
}

bgb_to_treedata <- function(results_path, geo_data_path, tree_path, area_names = NULL, nstates=NULL) {
  # load biogeobears results object
  load(results_path)
  # change data directories in results object
  res[["inputs"]]$geogfn <- geo_data_path
  res[["inputs"]]$trfn <- tree_path
  
  # read in tree
  tree <- RevGadgets::readTrees(paths = res[["inputs"]]$trfn)
  states <- res$inputs$all_geog_states_list_usually_inferred_from_areas_maxareas
  
  node_count=nrow(res$relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS)
  
  if (!is.null(nstates)){
    nstates=nstates
  }else{
    nstates=length(states)
  }
  
  column_count=(nstates*4)+3
  rev_data <- data.frame(matrix(nrow = node_count, ncol = column_count))
  colnames(rev_data) <- c(paste0("end_state_",1:nstates),
                          paste0("end_state_",1:nstates, "_pp"),
                          "end_state_other_pp",
                          paste0("start_state_",1:nstates),
                          paste0("start_state_",1:nstates, "_pp"),
                          "start_state_other_pp",
                          "node")
  
  # format df
  for (i in 1:node_count) {
    row <- res$ML_marginal_prob_each_state_at_branch_top_AT_node[i,]
    endstates<-order(row,decreasing = T)[1:nstates]
    endstates_pp <- row[order(row, decreasing=T)][1:nstates]
    if(nstates<length(row)){
      endstates_other_pp <- sum(row[order(row, decreasing=T)][nstates+1:length(row)], na.rm = T)
    }else{
      endstates_other_pp <- 0
    }
    
    row2 <- res$ML_marginal_prob_each_state_at_branch_bottom_below_node[i,]
    startstates<-order(row2,decreasing = T)[1:nstates]
    startstates_pp <- row2[order(row2, decreasing=T)][1:nstates]
    if(nstates<length(row2)){
      startstates_other_pp <- sum(row2[order(row2, decreasing=T)][nstates+1:length(row2)], na.rm = T)
    }else{
      startstates_other_pp <- 0
    }
    
    rev_data[i,]<-c(endstates, endstates_pp, endstates_other_pp, startstates, startstates_pp, startstates_other_pp, i)
  }
  
  # make better labels
  tipranges <- getranges_from_LagrangePHYLIP(res[["inputs"]]$geogfn)
  geo <- res$inputs$all_geog_states_list_usually_inferred_from_areas_maxareas
  geo_lab <- unlist(lapply(lapply(lapply(geo, function(x) x+1), as.character), paste0, collapse ="-"))
  label_dict <- data.frame(lab_num_short = 1:length(geo), lab_letters = geo_lab)
  tmp<-split_n_sort(label_dict$lab_letters)
  label_dict <- label_dict[tmp,]
  label_dict$col<-viridis::viridis(nrow(label_dict))
  
  # replace short number codes with labels 
  not_state_cols <- c(grep("_pp", colnames(rev_data)),
                      grep("node", colnames(rev_data)))
  state_cols <- c(1:ncol(rev_data))[!c(1:ncol(rev_data)) %in% not_state_cols]
  for (i in state_cols) { # loop through by column indices for the state columns
    col <- as.character(rev_data[,i])
    for (j in 1:length(col)) { # loop through each item in the column and replace with letters
      col[j] <- label_dict$lab_letters[which(label_dict$lab_num_short == col[j])]
    }
    rev_data[,i] <- col
  }
  
  #change "NA" to NA 
  rev_data %>%
    naniar::replace_with_na_all(condition = ~.x == "NA") -> rev_data
  # change any NAs in PP columns to 0 
  pp_cols <- grep("_pp", colnames(rev_data))
  for (p in pp_cols) { rev_data[ ,p][is.na(rev_data[ ,p])] <- 0 }
  
  # make treedata object (combine data and tree)
  tibble::as_tibble(tree[[1]][[1]]) %>%
    full_join(rev_data, by = 'node') %>%
    as.treedata() -> rev_treedata
  #add list of states
  attributes(rev_treedata)$state_labels <- as.character(na.omit(unique(unlist(rev_data[,state_cols]))))
  attributes(rev_treedata)$label_dict<-label_dict
  return(rev_treedata)
}