bgb_to_revgadgets <- function(results_path, geo_data_path, tree_path, area_names = NULL) {
  # load biogeobears results object
  load(results_path)
  # change data directories in results object
  res[["inputs"]]$geogfn <- geo_data_path
  res[["inputs"]]$trfn <- tree_path
  
  
  ##### Process data for plotting ##### 
  
  # read in tree separately
  
  tree <- RevGadgets::readTrees(paths = res[["inputs"]]$trfn)
  states <- res$inputs$all_geog_states_list_usually_inferred_from_areas_maxareas
  
  # create a dataframe with results in revgadgets compliant format
  rev_data <- data.frame(matrix(nrow = nrow(res$relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS),
                                ncol = 15))
  colnames(rev_data) <- c("end_state_1", "end_state_2", "end_state_3", 
                          "end_state_1_pp", "end_state_2_pp", "end_state_3_pp", 
                          "end_state_other_pp",
                          "start_state_1", "start_state_2", "start_state_3", 
                          "start_state_1_pp", "start_state_2_pp", "start_state_3_pp", 
                          "start_state_other_pp",
                          "node")
  # get end states
  for (i in 1:nrow(res$ML_marginal_prob_each_state_at_branch_top_AT_node)) {
    row <- res$ML_marginal_prob_each_state_at_branch_top_AT_node[i,]
    rev_data[i, 1] <- order(row,decreasing=T)[1]
    rev_data[i, 2] <- order(row,decreasing=T)[2]
    rev_data[i, 3] <- order(row,decreasing=T)[3]
    rev_data[i, 4] <- row[order(row,decreasing=T)[1]]
    rev_data[i, 5] <- row[order(row,decreasing=T)[2]]
    rev_data[i, 6] <- row[order(row,decreasing=T)[3]]
    rev_data[i, 7] <- sum(row[order(row,decreasing=T)[4:length(row)]]) 
  }
  # get start states
  for (i in 1:nrow(res$ML_marginal_prob_each_state_at_branch_bottom_below_node)) {
    row <- res$ML_marginal_prob_each_state_at_branch_bottom_below_node[i,]
    rev_data[i, 8] <- order(row,decreasing=T)[1]
    rev_data[i, 9] <- order(row,decreasing=T)[2]
    rev_data[i, 10] <- order(row,decreasing=T)[3]
    rev_data[i, 11] <- row[order(row,decreasing=T)[1]]
    rev_data[i, 12] <- row[order(row,decreasing=T)[2]]
    rev_data[i, 13] <- row[order(row,decreasing=T)[3]]
    rev_data[i, 14] <- sum(row[order(row,decreasing=T)[4:length(row)]]) 
  }
  rev_data$node <- 1:nrow(res$ML_marginal_prob_each_state_at_branch_bottom_below_node)
  
  # make better labels 
  tipranges <- getranges_from_LagrangePHYLIP(res[["inputs"]]$geogfn)
  
  geo <- res$inputs$all_geog_states_list_usually_inferred_from_areas_maxareas
  geo_num <- unlist(lapply(lapply(geo, as.character), paste0, collapse =""))
  if (is.null(area_names)) {
    area_names <- colnames(tipranges@df)
  }
  if (length(area_names) != length(colnames(tipranges@df))) stop("Number of specified area names is incorrect. Check your geo data file.")
  number_codes <- 0:(length(area_names)-1)
  geo_letters <- geo_num
  for (a in 1:length(area_names)) {
    geo_letters <- gsub(pattern = as.character(number_codes[a]), 
                        replacement = area_names[a], 
                        x = geo_letters)
  }

  label_dict <- data.frame(lab_num_short = 1:length(geo),
                           lab_num_long = geo_num,
                           lab_letters = geo_letters)
  
  # replace short number codes with letters 
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
  attributes(rev_treedata)$state_labels <- as.character(na.omit(unique(unlist(rev_data[,c(1:3, 8:10)]))))
  
  return(rev_treedata)
}



#==============================#
#### REQUIRED BGB FUNCTIONS ####
#==============================#

# setClass("tipranges", representation(df="data.frame"),
#          contains = "numeric"
# )
# 
# define_tipranges_object <- function(tmpdf=NULL)
# {
#   # Define the tipranges class;
#   # Now done separately
#   #tipranges_class = setClass("tipranges", representation(df="data.frame"))
#   
#   # Create an instance of the class
#   
#   # junk:
#   ## A class that extends the built-in data type "numeric"
#   # setClass("numWithId", representation(id = "character"),
#   #     contains = "numeric")
#   
#   
#   
#   # Make a simple default example df if needed
#   #	if (is.null(tmpdf))
#   #		{
#   # Default tip ranges data
#   areanames = c("A", "B", "C")
#   tipnames = c("tip1", "tip2", "tip3")
#   tip1_geog = c(1, 0, 0)
#   tip2_geog = c(0, 1, 0)
#   tip3_geog = c(0, 0, 1)
#   
#   # Make the temporary data frame
#   tmpdf2 = cbind(tip1_geog, tip2_geog, tip3_geog)
#   tmpdf2 = adf2(data.matrix(tmpdf2))
#   names(tmpdf2) = areanames
#   row.names(tmpdf2) = tipnames
#   #		}
#   
#   if (is.null(tmpdf))
#   {
#     tipranges_object = new("tipranges", df=tmpdf2)
#   } else {
#     tipranges_object = new("tipranges", df=tmpdf2)
#     tipranges_object@df = tmpdf
#   }
#   
#   # you can get the dataframe with
#   # tipranges_object@df
#   
#   
#   return(tipranges_object)
# }
# 
# adf <- function(x)
# {
#   return(as.data.frame(x, row.names=NULL, stringsAsFactors=FALSE))
# }
# 
# adf2 <- function(x)
# {
#   # Deals with the problem of repeated row names
#   rownames = 1:nrow(x)
#   return(as.data.frame(x, row.names=rownames, stringsAsFactors=FALSE))
# }
# 
# read_PHYLIP_data <- function(lgdata_fn="lagrange_area_data_file.data", regionnames=NULL)
# {
#   setup='
# 	lgdata_fn = "/Users/nickm/Desktop/__projects/_2011-07-15_Hannah_spider_fossils/_data/lagrange_for_nick2/palp_no_Lacun.data"
# 	'
#   
#   # Read the 1st line, split on whitespaces
#   firstline = scan(lgdata_fn, what="character", nlines=1)
#   
#   # Parse the firstline
#   ntips = as.numeric(firstline[1])
#   nareas = as.numeric(firstline[2])
#   
#   #######################################################
#   # If the length of the first line is > 2, parse that information 
#   #######################################################
#   if (length(firstline) > 2)
#   {
#     # Get region names	
#     regionnames = firstline[3:length(firstline)]
#     
#     # Remove "(" and ")"
#     regionnames = mapply(sub, pattern="\\(", replacement="", x=regionnames)
#     regionnames = mapply(sub, pattern="\\)", replacement="", x=regionnames)
#     names(regionnames) = NULL
#   }
#   
#   # Make the data.frame
#   tmpdf = matrix(data=NA, nrow=ntips, ncol=nareas)
#   
#   # Parse the remaining rows
#   tmplines = scan(lgdata_fn, what="character", sep="\t", skip=1)
#   tmplines = matrix(data=tmplines, ncol=2, byrow=TRUE)
#   tmplines
#   
#   tipnames = tmplines[,1]
#   areas_char = tmplines[,2]
#   
#   areas_char2 = unlist(mapply(strsplit, x=areas_char, split=""))
#   areas_char3 = matrix(data=areas_char2, nrow=ntips, byrow=TRUE)
#   
#   
#   # Store in a data.frame and return
#   tmpdf = areas_char3
#   tmpdf = adf(tmpdf)
#   
#   if ( is.null(regionnames) )
#   {
#     alphabet = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z")
#     regionnames	= alphabet[1:ncol(tmpdf)]
#     colnames(tmpdf) = regionnames
#   } else {
#     names(tmpdf) = regionnames
#   }
#   row.names(tmpdf) = tipnames
#   
#   return(tmpdf)	
# }
# 
# getranges_from_LagrangePHYLIP <- function(lgdata_fn="lagrange_area_data_file.data")
# {
#   # Make a tipranges instance
#   #tipranges_object = define_tipranges_object()
#   #tipranges_object@df = dfnums_to_numeric(tipranges_object@df)
#   
#   # Read the Lagrange geographic ranges data file
#   tmp_blah = read_PHYLIP_data(lgdata_fn)
#   tmp_input = adf2(data.matrix(tmp_blah))
#   #nums_as_char = as.numeric(unlist(tmp_input))
#   #tmpdf = adf2(matrix(data=nums_as_char, nrow=nrow(tmp_input), ncol=ncol(tmp_input), byrow=TRUE))
#   #names(tmpdf) = names(tmp_input)
#   #rownames(tmpdf) = rownames(tmp_input)
#   #sum(tmpdf)
#   #tmpdf = dfnums_to_numeric(tmp_input, printout=TRUE)
#   #cls.df(tmpdf)
#   
#   # Put into a tipranges object
#   tipranges_object = define_tipranges_object(tmpdf=tmp_input)
#   # you can get the dataframe with
#   # tipranges_object@df
#   
#   tipranges_object@df = adf2(data.matrix(tipranges_object@df))
#   
#   rownames(tipranges_object@df) = rownames(tmp_blah)
#   
#   return(tipranges_object)
# }
