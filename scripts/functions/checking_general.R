check_trueOrFalse <- function(arg){
  msg <- paste0("`", deparse(substitute(arg)), "` must be either TRUE or FALSE")
  if(is.list(unclass(arg))){stop(msg)}
  if(length(arg) != 1){stop(msg)}
  if(!is.logical(arg)){stop(msg)}
  if(is.na(arg)){stop(msg)}
}

check_replacementList <- function(data, object){
  if(is.null(object$cols)){stop()}
  if(is.list(unclass(object$cols))){stop()}
  if(!is.character(object$cols)){stop()}
  cols <- na.omit(object$cols)
  if(!all(cols %in% colnames(data))){stop()}

  if(is.null(object$old_new_pairs)){stop()}
  if(is.list(unclass(object$old_new_pairs))){stop()}
  if(!is.character(object$old_new_pairs)){stop()}
  pairs <- na.omit(object$old_new_pairs)
  if(length(pairs) < 1){stop()}
  pair_names <- names(pairs)
  pair_names <- pair_names[pair_names != ""]
  if(length(pair_names) != length(pairs)){stop()}
}
