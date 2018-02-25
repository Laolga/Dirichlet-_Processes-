MDP\_example
================
Olga Lazareva
12/9/2017

A discrimination problem
------------------------

First of all let's generate *X*<sub>*i*1</sub>, ...*X*<sub>*i**k*<sub>*i*</sub></sub> coming from *D*(*p*<sub>1</sub>)...*D*(*p*<sub>*k*</sub>) and *Y* ∼ *D*(*p*<sub>*j*</sub>) where 1 ≤ *j* ≤ *k*. We are assuming that *α*<sub>1</sub>, ..*α*<sub>*k*</sub> ∼ *P**o**i**s*(*λ*<sub>*k*</sub>) with the same support \[0,100\]. And we would like to minimize misclassification risk *r*<sub>*i*</sub> of *Y*<sub>*j*</sub> which is proportional to $P(Y|i) = \\prod\_{j=1}^n\\frac{\\alpha\_{i,j}^{'(k\_j)}}{M^{(n)}}$, where (here goes definition of alpha and so on)

``` r
k = 10 #number of samples X from different distributions we consider
lambda =seq(1,100,length.out = k)
M = 50
theta = seq(1,100,length= 100) #support
n = length(theta)
fun1 = function(n,lambda)(rpois(n,lambda))
sethuraman.cost <- function(fun,mu, M){
  n <- 5000
  y <- fun(n,mu)
  thet <- rbeta(n,shape1 = 1, shape2 = M)
  prob <- rep(0,n)
  prob[1] <- thet[1]
  prob[2:n] = sapply(2:n, function(i) thet[i] * prod(1 - thet[1:(i-1)]))
  return(list('alpha'=y,'pi'=prob))
}

X = matrix(NA, n,k)
j = sample(seq(2,k-2),1) #random choice of distribution for y
for (i in 1:k){
  obj= sethuraman.cost(fun1,mu = lambda[i],M = M)
  X[,i] <- sample(obj$alpha,size= n, prob=obj$pi,replace=T)
  if (i == j)(obj_y = obj)
}  
y <- sample(obj_y$alpha,size= n, prob=obj_y$pi,replace=T)
```

Let's compare disribution of *X*<sub>*j* − 1</sub>, *X*<sub>*j*</sub>, *X*<sub>*j* + 1</sub> and *Y*:

``` r
plot(density(X[,j]),col = "blue",lwd =2,
     xlim=range(min(X[,j-1]),max(X[,j+1])),
     ylim=range(0,max(c(max(density(X[,j-1])$y),max(density(X[,j])$y),max(density(X[,j+1])$y)))),
                main = "Samples distributions")
lines(density(X[,j-1]),col = "skyblue1",lwd =2)
lines(density(X[,j+1]),col = "skyblue1",lwd =2)
lines(density(y),lwd = 3)
legend("topright", 
       c("X(j)","X(j-1)","X(j+1)","Y"), col = c("blue", "skyblue1", "skyblue1","black"),
       lwd = 2)
```

![](dicriminative_problem_files/figure-markdown_github/unnamed-chunk-2-1.png)

Now we can compute risk *r*<sub>*i*</sub> for each of the distributions:

``` r
C = matrix(NA,n,k)
k_ = rep(NA,length.out = length(theta))
k_ = sapply(1:length(theta), function(i) sum(y==i)) #number of Y_i = i
for (i in 1:k){
  C[,i] = sapply(1:length(theta), function(m) sum(X[,i] == m))
}
alpha = matrix(NA,length(theta),k)
for (l in 1:k){
  alpha[,l] = dpois(theta,lambda[l])*M
  alpha[,l]+C[,l] #alpha'
}

poli = function(x,n){
  p = x
  for (i in 1:(n-1)){
    p = prod(p,x+i)
  }
  if (n ==0)(p = 1)
  return(p)
}
probs = matrix(NA,length(theta),length(lambda)) #risk
for (i in theta){
  for (l in 1:k){
    probs[i,l] = log(poli(alpha[i,l],k_[i])/poli(M,length(theta)))
  }
}
probs = sapply(1:length(lambda),function(i)prod(probs[,i]))
plot(probs,ylab = "log(risk)", xlab = "model",col = "blue", lwd = 4)
abline(v = j, lwd=3, lty=2)
points(x = which.min(probs),y = probs[which.min(probs)], col = "red", lwd = 3)
legend("topright", 
       c("risk_i","j","minimum risk"), col = c("blue", "black","red"),
       lwd = 2)
```

![](dicriminative_problem_files/figure-markdown_github/unnamed-chunk-3-1.png)

As we see, the minimum value of the risk correspondes to the case where *i* = *j* which means that the distribution of *Y* was classified correctly
