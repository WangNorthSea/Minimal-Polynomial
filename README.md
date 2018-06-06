# Minimal-Polynomial
one possible algorithm to calculate the minimal polynomial of a matrix

For a square matrix belongs to Mn(F) (F represents a field, a concept in abstract algebra), there must be some polynomials which belong to F[t] (and it's not "0") can make it be the zero matrix, and we call these polynomials "annihilator polynomials of square matrix A". This means that if you put A into these polynomials to replace the original undetermined element, you will get the zero matrix.

We can prove that there must be a polynomial which belongs to the set of annihilator polynomials of square matrix A has the lowest degree and the coefficient of its term which has the highest degree is 1, and we call this polynomial "the minimal polynomial of square matrix A".

Usually, it's not easy to calculate the minimal polynomial of a square matrix A. According to Cayley-Hamilton Theorem and its inferences, for a n-order square matrix A, the degree of its minimal polynomial must be equal or lower than n. To get it, you can calculate A^2, A^3......A^n, and find the linear correlated combination of <En, A, A^2.......> which has the lowest degree, or you can calculate the characteristic polynomial of A firstly, which is det(tEn - A), and combine its factors to check whether the new polynomial you get is an annihilator polynomial of A and has the lowest degree at the same time.

However, these two algorithms are too obscure ,and it's very hard to describe them accurately to make a steady algorithm for computers. To solve this problem, I came up with an algorithm which is easy to be described by programming languages (and I have made it come true by using Swift).

Let me use a random 3-order square matrix to introduce this algorithm for you, and we call such a matrix "A". At first, we know that its minimal polynomial must be in a form like this "μA(t) = x3t³ + x2t² + x1t + x0" (x3, x2, x1, x0 belong to F, and they won't all be zero), so what we are gonna do is to get the values of x3, x2, x1, x0. Put A into this polynomial, and because μA(A) equals to zero matrix, you can get a matrix which we can consider as an augmented matrix of some linear equations (limited by the format of this README file, I cannot show you this process by picture, but you can write it down on paper, it's not difficult). 

Then, we are gonna get the values of x3, x2, x1, x0 from this augmented matrix. However, this augmented matrix actually comes from homogeneous linear equations and its rank is smaller than the number of its undetermined elements, so in fact, the number of its possible solutions is infinite. Fortunately, this doesn't mean it's useless for us. We know that we can use some parameters like α1, α2, α3...... to represent all possible solutions of homogeneous linear equations which have infinite possible solutions. You give α1, α2, α3.... different numbers, and you'll get different solutions. So here comes the question, how to get such a general solution?

Firstly, as usual, we do Gauss Elimination to the augmented matrix we have. Then we will get an upper triangle matrix, and let's call it "Triangle". We can find some interesting things from "Triangle". Let's see it from its last row to its first row. First, if a row of it has only one non-zero element, the value of the undetermined element whose coefficient is this non-zero element will be zero. Second, if a row of it has more than one non-zero element, the number of parameters α1, α2, α3....... will be the number of non-zero elements(exclude which is the coefficient of undetermined element whose value is 0 according to "First") minus 1.

So far, we have got the number of parameters in the general solution of this augmented matrix. This means that we can use these parameters to represent the value of each undetermined element. The only difference between different undetermined elements is the coefficient of each parameter. We can use a two dimensional array to save the results we get during this process. The number of elements in this two dimensional array is the number of undetermined elements, and each element of this array saves the coefficients of parameters for each undetermined element. So each element of this array is a one dimensional array and the number of its elements is the number of parameters, which we have got. By using this two dimensional array, we can get the general solution of the augmented matrix I mentioned before.

The next step is very, very important. We know that the minimal polynomial has the lowest degree compared to annihilator polynomials of "A", and the coefficient of its term which has the highest degree is 1. We should keep this rule as basic principle to realize in the next step.

Now, let's go back to x3, x2, x1, x0. We have got the coefficient of each parameter for each x. In order to make the degree be the lowest, we assume that x3 = 0, and we consider its coefficients for parameters as a row vector and we do so to x2, x1, x0, too. If the rank of the matrix formed by the row vector from x3 and another row vector from x2 or x1 or x0 is 1, or in other words, is the same as the matrix formed by the row vector from x3 only, we will know that if x3 = 0, then x2 or x1 or x0 will be zero, too. In other words, if the row vector from x2 or x1 or x0 is linear correlated with the row vector from x3, then x2 or x1 or x0 will also be zero.

If x3 = 0 wouln't make x2, x1, x0 all have to be zero, then x3 has to be zero. After that, we keep the row vector from x3 in a new matrix, and in the situation x3 = 0, we assume that x2 = 0 to repeat this process until we find a x cannot be zero. If it cannot be zero, it has to be 1. Now we have got a new matrix whose last row is coefficients of parameters for the x which is 1 and other rows are coefficients of parameters for x which are 0. Then we consider the values of x we have known as a column vector, and put it into the new matrix as a new column. After this step, we will have another "augmented matrix" which has a unique solution and this unique solution contains the value of each parameter.

When we have got the value of each parameter, the last step will be using these values to calculate the value of each undetermined element, which is very easy. We can save the value of each x into a one dimensional array which has four elements: x3, x2, x1 and x0, and they are coefficients of the minimal polynomial of "A", which we have spent great efforts to get.
