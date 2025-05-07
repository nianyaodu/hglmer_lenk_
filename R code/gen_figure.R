##### Code for generating figure and table of report ######

### Generate figure 1
library(mvtnorm)
n=100;m=10
psi = c(0.3,0.7)
thetas <- list(c(-10,7),c(5,5))
Lambdas <- list(matrix(c(25,9,9,4),nrow=2,byrow=T),
                matrix(c(9,-5,-5,5),nrow=2,byrow=T))
alpha=-1; tau2 = 4

Simdata <- list(length=n);Simind <- vector(length=n)
Simbeta <- matrix(NA_real_,nrow=n,ncol=2)

set.seed(3)
for(i in 1:n){
  Simind[i] <- ind <- sample(2,1,prob=psi)
  theta <- thetas[[ind]];
  Lambda <- Lambdas[[ind]];
  Simbeta[i,] <- beta <- rmvnorm(1,mean=theta,sigma=Lambda)
  X <- cbind(rep(1,m),rnorm(m,2,1))
  phi <- rnorm(1,alpha,sqrt(tau2))
  y <- X%*%t(beta) + rnorm(m,mean=0,sd=sqrt(exp(phi)))
  data <- cbind(y,X)
  Simdata[[i]] <- data
}

betas <- list()
for(k in 1:2){
  betas[[k]] <- Simbeta[Simind==k,]
}

table(Simind)

## generate figure

pdf(file = "fig_1.pdf",width = 5, height = 5)

par(mfrow=c(1,1))
plot(betas[[1]][,1],betas[[1]][,2],pch=1,col="red",xlim=c(-25,15),ylim=c(-5,15), xlab=expression(beta[0]), ylab=expression(beta[1]))
points(betas[[2]][,1],betas[[2]][,2],pch=2, col="blue")
points(thetas[[1]][1],thetas[[1]][2],pch=16, cex=2, col="red")
points(thetas[[2]][1],thetas[[2]][2],pch=17,cex=2,col="blue")
text(-10,7, "(-10,7)", pos=3, cex=1.5)
text(5,5, "(5,5)", pos=1, cex=1.5)

dev.off()



### Generate figure 2

load(file="hglmer_linear_K2_15chains.Rdata")

result = Result_fit_K2
# sample sets
well <- 2; switch <- 8

pdf(file = "fig_2.pdf",width = 8, height = 6)

par(mfrow=c(2,2),mar = c(2, 2, 2, 2))

### psi
plot(result[[well]]$psi[,1], type="l", col="black", 
     ylim=c(0,1), xlab="iterations",ylab="",main=expression(psi))
lines(result[[well]]$psi[,2], col="black")
abline(h=psi[1], col="black",lwd=2); abline(h=psi[2], col="black",lwd=2)

### theta[[1]]
plot(result[[well]]$theta[[1]][,1], type="l", col="black",
     ylim=c(-12,9),xlab="iterations",ylab="", main=expression(theta[1]))
lines(result[[well]]$theta[[1]][,2], col="black")
abline(h=thetas[[1]][1], col="black",lwd=2)
abline(h=thetas[[1]][2], col="black",lwd=2)

### theta[[2]]
plot(result[[well]]$theta[[2]][,1], type="l", col="black",
     ylim=c(-12,9),xlab="iterations",ylab="", main=expression(theta[2]))
lines(result[[well]]$theta[[2]][,2], col="black")
abline(h=thetas[[2]][1], col="black",lwd=2)
abline(h=thetas[[2]][2], col="black",lwd=2)

### alpha and tau2
plot(result[[well]]$alpha, type="l", ylim=c(-1.2,7),
     xlab="iterations",ylab="",main=expression(paste(alpha," & ",tau^2,sep="")))
abline(h=alpha, lwd=2)
lines(result[[well]]$tau2, type="l", ylim=c(0,180),
     xlab="iterations",ylab="")
abline(h=tau2, lwd=2)

dev.off()


mean(result[[well]]$alpha);sd()
mean(result[[well]]$tau2)

### Generate figure 3

pdf(file = "fig_3.pdf",width = 8, height = 6)

par(mfrow=c(2,2),mar = c(2, 2, 2, 2))

### psi
plot(result[[switch]]$psi[,1], type="l", col="black", 
     ylim=c(0,1), xlab="iterations",ylab="",main=expression(psi))
lines(result[[switch]]$psi[,2], col="black")
abline(h=psi[1], col="black",lwd=2); abline(h=psi[2], col="black",lwd=2)

### theta[[1]]
plot(result[[switch]]$theta[[1]][,1], type="l", col="black",
     ylim=c(-12,9),xlab="iterations",ylab="", main=expression(theta))
lines(result[[switch]]$theta[[1]][,2], col="black")
abline(h=thetas[[1]][1], col="black",lwd=2)
abline(h=thetas[[1]][2], col="black",lwd=2)

### theta[[2]]
plot(result[[switch]]$theta[[2]][,1], type="l", col="black",
     ylim=c(-12,9),xlab="iterations",ylab="", main=expression(theta))
lines(result[[switch]]$theta[[2]][,2], col="black")
abline(h=thetas[[2]][1], col="black",lwd=2)
abline(h=thetas[[2]][2], col="black",lwd=2)

### alpha and tau2
plot(result[[switch]]$alpha, type="l", ylim=c(-1.2,7),
     xlab="iterations",ylab="",main=expression(paste(alpha," & ",tau^2,sep="")))
abline(h=alpha, lwd=2)
lines(result[[switch]]$tau2, type="l", ylim=c(0,180),
     xlab="iterations",ylab="")
abline(h=tau2, lwd=2)


dev.off()


## Generate table 1

K = 2

## mean

### psi
psi_mean <- matrix(NA, nrow=15, ncol=K)
for (i in 1:15){
  psi.temp <- result[[i]]$psi
  psi_mean[i,] <- round(apply(psi.temp,2,mean),2)
}
true.psi <- matrix(psi,nrow=15,ncol=K,byrow=T)
psi_compare_mean <- cbind(psi_mean,true.psi)
colnames(psi_compare_mean) <- c(paste0("hat_psi.",seq(K)),paste0("psi.",seq(K)))


### theta
theta_mean <- theta_compare_mean <- list()
for (k in 1:K){
  theta_mean[[k]] <- matrix(NA, nrow=15, ncol=2)
}
for (k in 1:K){
  for (i in 1:15){
    theta.temp<- result[[i]]$theta[[k]]
    theta_mean[[k]][i,] <- round(apply(theta.temp,2,mean),2)
  }
}
for (k in 1:K){
  true.theta <- matrix(thetas[[k]],nrow=15,ncol=2,byrow=T)
  theta_compare_mean[[k]] <- cbind(theta_mean[[k]],true.theta)
  colnames(theta_compare_mean[[k]])<-
    c(paste0("hat_theta.",seq(2)),paste0("theta.",seq(2)))
}

### Lambdas
Lambda_mean <- Lambda_compare_mean <- list(list(),list())
for (k in 1:K){
  for (i in 1:15){
    Lambda_mean[[k]][[i]] <- matrix(0,nrow=2,ncol=2)
  }
}
for (k in 1:K){
  for (i in 1:15){
    Lambda_mean[[k]][[i]][1,1]<- round(mean(result[[i]]$Lambda[[k]][1,1,]),2)
    Lambda_mean[[k]][[i]][1,2]<- round(mean(result[[i]]$Lambda[[k]][1,2,]),2)
    Lambda_mean[[k]][[i]][2,1]<- round(mean(result[[i]]$Lambda[[k]][2,1,]),2)
    Lambda_mean[[k]][[i]][2,2]<- round(mean(result[[i]]$Lambda[[k]][2,2,]),2)
  }
}
for (k in 1:K){
  for (i in 1:15){
    Lambda_compare_mean[[k]][[i]] <- cbind(Lambda_mean[[k]][[i]],Lambdas[[k]])
    colnames(Lambda_compare_mean[[k]][[i]]) <- c("","hat_Lambda","","Lambda")
  }
}

### alpha
alpha_mean <- matrix(NA, nrow=15, ncol=1)
for (i in 1:15){
  alpha.temp <- result[[i]]$alpha
  alpha_mean[i] <- round(mean(alpha.temp),2)
}
alpha_compare_mean <- cbind(alpha_mean,rep(alpha,15))
colnames(alpha_compare_mean) <- c("hat_alpha","alpha")


### tau2
tau2_mean <- matrix(NA, nrow=15, ncol=1)
for (i in 1:15){
  tau2.temp <- result[[i]]$tau2
  tau2_mean[i] <- round(mean(tau2.temp),2)
}
tau2_compare_mean <- cbind(tau2_mean,rep(tau2,15))
colnames(tau2_compare_mean) <- c("hat_tau2","tau2")



# well estimated ---------------------------------------------

psi_compare_mean[well,]
alpha_compare_mean[well,]
tau2_compare_mean[well,]

theta_compare_mean[[1]][well,]
theta_compare_mean[[2]][well,]

Lambda_compare_mean[[1]][[well]]
Lambda_compare_mean[[2]][[well]]



## sd


### theta
theta_sd <- theta_compare_sd <- list()
for (k in 1:K){
  theta_sd[[k]] <- matrix(NA, nrow=15, ncol=2)
}
for (k in 1:K){
  for (i in 1:15){
    theta.temp<- result[[i]]$theta[[k]]
    theta_sd[[k]][i,] <- round(apply(theta.temp,2,sd),2)
  }
}
for (k in 1:K){
  true.theta <- matrix(thetas[[k]],nrow=15,ncol=2,byrow=T)
  theta_compare_sd[[k]] <- cbind(theta_sd[[k]],true.theta)
  colnames(theta_compare_sd[[k]])<-
    c(paste0("hat_theta.",seq(2)),paste0("theta.",seq(2)))
}

### Lambdas
Lambda_sd <- Lambda_compare_sd <- list(list(),list())
for (k in 1:K){
  for (i in 1:15){
    Lambda_sd[[k]][[i]] <- matrix(0,nrow=2,ncol=2)
  }
}
for (k in 1:K){
  for (i in 1:15){
    Lambda_sd[[k]][[i]][1,1]<- round(sd(result[[i]]$Lambda[[k]][1,1,]),2)
    Lambda_sd[[k]][[i]][1,2]<- round(sd(result[[i]]$Lambda[[k]][1,2,]),2)
    Lambda_sd[[k]][[i]][2,1]<- round(sd(result[[i]]$Lambda[[k]][2,1,]),2)
    Lambda_sd[[k]][[i]][2,2]<- round(sd(result[[i]]$Lambda[[k]][2,2,]),2)
  }
}
for (k in 1:K){
  for (i in 1:15){
    Lambda_compare_sd[[k]][[i]] <- cbind(Lambda_sd[[k]][[i]],Lambdas[[k]])
    colnames(Lambda_compare_sd[[k]][[i]]) <- c("","hat_Lambda","","Lambda")
  }
}

### alpha
alpha_sd <- matrix(NA, nrow=15, ncol=1)
for (i in 1:15){
  alpha.temp <- result[[i]]$alpha
  alpha_sd[i] <- round(sd(alpha.temp),2)
}
alpha_compare_sd <- cbind(alpha_sd,rep(alpha,15))
colnames(alpha_compare_sd) <- c("hat_alpha","alpha")


### tau2
tau2_sd <- matrix(NA, nrow=15, ncol=1)
for (i in 1:15){
  tau2.temp <- result[[i]]$tau2
  tau2_sd[i] <- round(sd(tau2.temp),2)
}
tau2_compare_sd <- cbind(tau2_sd,rep(tau2,15))
colnames(tau2_compare_sd) <- c("hat_tau2","tau2")



# well estimated ---------------------------------------------

alpha_compare_sd[well,]
tau2_compare_sd[well,]

theta_compare_sd[[1]][well,]
theta_compare_sd[[2]][well,]

Lambda_compare_sd[[1]][[well]]
Lambda_compare_sd[[2]][[well]]

