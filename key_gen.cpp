
/*

This file will:
    - generate the private key
        -> generate vectors of size n with entries between d and -d
        -> test how orthogonal they are by calculating the hadamard ratio
        -> set them up in n by n matrix V
    - generate the public key
        -> create an n by n unitary matrix U; det(U) = 1 or det(U) = -1
        -> "One way to create U is as a product of a large number of randomly chosen elementary matrice"
        -> compute W = UV
        -> the rows w_0, w_1, w_2, ... are now the public key 
*/









