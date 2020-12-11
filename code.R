#1. f <-function(v) c(v[1]^(3/4)*v[2]^(1/3)+v[1]-24, v[1]^(1/4)*v[2]^(2/3)+v[2]-24)
#Make a good initial guess
v0 = c(16,8); f(v0)
#Calculate the 2 x 2 Jacobian matrix
A <- jacobian(f, v0); A
#Invert the Jacobian
AInv =solve(A); AInv
#Use the same update formula as in the single-variable case
v1 <- v0 - AInv%*%f(v0); v1; f(v1)
#Repeat to improve the approximation
A <- jacobian(f, v1)
v2 <- v1 - solve(A)%*%f(v1); v2; f(v2)

A <- jacobian(f, v2)
v3 <- v2 - solve(A)%*%f(v2); v3; f(v3)

#2. library(numDeriv)    #for the grad() function
par(mar = c(2,2,1,1))  #maximize space for graphs
#Topic 1 - Single variable

#An easy case -- the square root of 2
#Define a function to use in an equation f(x) = 0
f <- function(x) x*(x^2 - 1)*(x^2 - 4) -1
#Graph the function on a interval where it has a root
curve(f(x), from = -2, to = 2.5)
abline(h=0, col = "green")


#Find a value that makes f close to zero  and evaluate f and f' there test -2
x0 <- -2; f(x0); grad(f, x0)
#Add the tangent line to the plot
curve(f(x0) + grad(f, x0)*(x-x0), col = "red", add = TRUE)
#Solve a linear approximation to find where the tangent line is zero
x1 <- x0 - f(x0)/grad(f, x0); abline(v = x1, col = "red"); f(x1)
#Repeat to improve the approximation
x2 <- x1 - f(x1)/grad(f, x1); abline(v = x2, col = "red", lty = 2); f(x2)
x3 <- x2 - f(x2)/grad(f, x2); abline(v = x3, col = "red", lty = 3); f(x3)
x3

#Find a value that makes f close to zero  and evaluate f and f' there test -1
x0 <- -1; f(x0); grad(f, x0)
#Add the tangent line to the plot
curve(f(x0) + grad(f, x0)*(x-x0), col = "red", add = TRUE)
#Solve a linear approximation to find where the tangent line is zero
x1 <- x0 - f(x0)/grad(f, x0); abline(v = x1, col = "red"); f(x1)
#Repeat to improve the approximation
x2 <- x1 - f(x1)/grad(f, x1); abline(v = x2, col = "red", lty = 2); f(x2)
x3 <- x2 - f(x2)/grad(f, x2); abline(v = x3, col = "red", lty = 3); f(x3)
x3

#Find a value that makes f close to zero  and evaluate f and f' there. test 0
x0 <- 0; f(x0); grad(f, x0)
#Add the tangent line to the plot
curve(f(x0) + grad(f, x0)*(x-x0), col = "red", add = TRUE)
#Solve a linear approximation to find where the tangent line is zero
x1 <- x0 - f(x0)/grad(f, x0); abline(v = x1, col = "red"); f(x1)
#Repeat to improve the approximation
x2 <- x1 - f(x1)/grad(f, x1); abline(v = x2, col = "red", lty = 2); f(x2)
x3 <- x2 - f(x2)/grad(f, x2); abline(v = x3, col = "red", lty = 3); f(x3)
x3

#Find a value that makes f close to zero  and evaluate f and f' there. Test 1
x0 <- 1; f(x0); grad(f, x0)
#Add the tangent line to the plot
curve(f(x0) + grad(f, x0)*(x-x0), col = "red", add = TRUE)
#Solve a linear approximation to find where the tangent line is zero
x1 <- x0 - f(x0)/grad(f, x0); abline(v = x1, col = "red"); f(x1)
#Repeat to improve the approximation
x2 <- x1 - f(x1)/grad(f, x1); abline(v = x2, col = "red", lty = 2); f(x2)
x3 <- x2 - f(x2)/grad(f, x2); abline(v = x3, col = "red", lty = 3); f(x3)
x3

#Find a value that makes f close to zero  and evaluate f and f' there. Now let's test near 2
x0 <- 2; f(x0); grad(f, x0)
#Add the tangent line to the plot
curve(f(x0) + grad(f, x0)*(x-x0), col = "red", add = TRUE)
#Solve a linear approximation to find where the tangent line is zero
x1 <- x0 - f(x0)/grad(f, x0); abline(v = x1, col = "red"); f(x1)
#Repeat to improve the approximation
x2 <- x1 - f(x1)/grad(f, x1); abline(v = x2, col = "red", lty = 2); f(x2)
x3 <- x2 - f(x2)/grad(f, x2); abline(v = x3, col = "red", lty = 3); f(x3)
x3
