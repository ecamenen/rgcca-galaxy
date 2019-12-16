check_integer_func <- function(x,  type = "scalar", float = FALSE, min = 1){
    assign(x, check_integer(x, get(x, parent.frame()), type, float, min), 1)
}
