# Verifies that the input vectors are valid
check_vectors <- function(x, y, n){

	if(all(is.vector(x), is.vector(y), is.vector(n))){

		if(length(x) == length(y) & length(y) == length(n)){

			if(sum(is.na(c(x, y, n))) == 0){

				if(all(n >= y)){
					return(TRUE)
				}
				else
					stop("Successes (y) cannot exceed number of trials (n).")
			}
			else
				stop("Missing or NA values in data.")
		}
		else
			stop("Ensure all vectors are the same length.")
	}
	else
		stop("Inputs for x, y, and n must be of class 'vector'.")
}
