# Minimal-Polynomial
one possible algorithm to calculate the minimal polynomial of a matrix

For a square matrix belongs to Mn(F) (F represents a field, a concept in abstract algebra), there must be some polynomials which belong to F[t] (and it's not "0") can make it be the zero matrix, and we call these polynomials "annihilator polynomials of square matrix A". This means that if you put A into these polynomials to replace the original undetermined element, you will get the zero matrix.

We can prove that there must be a polynomial which belongs to the set of annihilator polynomials of square matrix A has the lowest degree and the coefficient of its term which has the highest degree is 1, and we call this polynomial "the minimal polynomial of square matrix A".

Usually, it's not easy to calculate the minimal polynomial of a square matrix A. According to Cayley-Hamilton Theorem and its inferences, for a n-order square matrix A, the degree of its minimal polynomial must be equal or lower than n. To get it, you can calculate A^2, A^3......A^n, and find the linear correlated combination of <En, A, A^2.......> which has the lowest degree, or you can calculate the characteristic polynomial of A firstly, which is det(tEn - A), and combine its factors to check whether the new polynomial you get is an annihilator polynomial of A and has the lowest degree at the same time.

However, these two algorithms are too obscure ,and it's very hard to describe them accurately to make a steady algorithm for computers. To solve this problem, I came up with an algorithm which is easy to be described by programming languages (and I have made it come true by using Swift).

Let me use any 3-order square matrix to introduce this algorithm for you, and we call such a matrix "A".
