firstContinuous <- function(x) {
	return( Position( function(p) p > 2,
										x,
										right = FALSE,
										nomatch = NA_integer )  
	)
}
