# Copyright (C) 2011 Pierrick Bruneau, see README for full notice

classicEM <- function(data, k, thres=0.1, maxit=NULL) {
	# initialize components to random position and data-scaled cov matrix
	n <- length(data[,1])
	d <- length(data[1,])
	
	maxs <- numeric(d)
	mins <- numeric(d)
	
	for(i in 1:d) {
		maxs[i] <- max(data[,i])
		mins[i] <- min(data[,i])
	}
	
	w <- numeric(k)
	mean <- list()
	cov <- list()
	

	
	sdvect <- numeric(d)
	for(i in 1:d) {
		sdvect[i] <- sd(data[,i])
	}
	
	for (i in 1:k) {
		w[i] <- 1/k
		mean[[i]] <- matrix(0, nrow=d, ncol=1)
		for(j in 1:d) {
			mean[[i]][j,1] <- runif(1, mins[j], maxs[j])
		}
	
		
		cov[[i]] <- diag(((maxs-mins)^2) / 64)
	}

	initmean <- mean
	initcov <- cov
	
	# start iterations
	# we detect covergence by computing argmax at each step, stop if no change.
	
	conv <- FALSE
	labels <- rep(1,n)
	resp <- matrix(0, nrow=n, ncol=k)
	lnresp <- matrix(0, nrow=n, ncol=k)
	oldresp <- resp
	
	nk <- list()
	agitation <- list()
	
	it <- 0
	rate <- 0

	likelihood <- -Inf
	#borneInf <- -Inf
	deltaLikelihood <- Inf

	noupdate <- numeric()
	
	while(!conv) {
		
		it <- it + 1
		print(paste("iteration", it, ", likelihood=", likelihood))

		
		oldlabs <- labels
		oldresp <- resp

		for(j in 1:k) {
			if(length(which(noupdate==j)) == 0) {
				resp[,j] <- w[j] * dmnorm(data, mean=as.numeric(mean[[j]]), varcov=cov[[j]])
				lnresp[,j] <- log(w[j]) + dmnorm(data, mean=as.numeric(mean[[j]]), varcov=cov[[j]], log=TRUE)
			} else {
				resp[,j] <- rep(0, n)
				lnresp[,j] <- rep(-Inf, n)
			}
		}
		
		for(i in 1:n) {
			if(sum(resp[i,]) > 0) {
				resp[i,] <- resp[i,] / sum(resp[i,])
				labels[i] <- which.max(resp[i,])
			} else {
				resp[i,] <- 0
				themax <- which.max(lnresp[i,])
				resp[i,themax] <- 1
				labels[i] <- themax
			}
		}
		
		
		

		# nouvelle règle : quand composante ajuste trop peu de données (<2), on cesse de la mettre à jour.
		for(i in 1:k) {
			if(length(which(noupdate==i)) == 0) {
				w[i] <- sum(resp[,i])
				mean[[i]] <- matrix(0, nrow=d, ncol=1)
			
				if(w[i] == 0) {
					print("null component detected")
					mean[[i]] <- initmean[[i]]
					cov[[i]] <- initcov[[i]]
					noupdate <- c(noupdate, i)
					likelihood <- -Inf
				} else {
					for(j in 1:n) {
						mean[[i]] <- mean[[i]] + resp[j,i] * as.matrix(data[j,])
					}
					mean[[i]] <- mean[[i]] / w[i]

					cov[[i]] <- matrix(0, nrow=d, ncol=d)
					for(j in 1:n) {
						cov[[i]] <- cov[[i]] + resp[j,i] * (as.matrix(data[j,]) - mean[[i]]) %*% t(as.matrix(data[j,]) - mean[[i]])
					}
					cov[[i]] <- cov[[i]] / w[i]
					
					if (w[i] < 2) {
						print("useless component detected")
						w[i] <- 0
						noupdate <- c(noupdate, i)
						likelihood <- -Inf					
					}
				}

			
				# control degenerate covariances
				# condition number : ratio of largest to smallest eigenvalue (http://tolstoy.newcastle.edu.au/R/e2/help/06/11/4890.html)
				# decreasing order
				vals <- eigen(cov[[i]], only.values=TRUE)$values
				cond <- vals[length(vals)] / vals[1]
				if(is.na(cond)) {
					print(cov[[i]])
					print(vals)
					print(w[i])
				}
				if(cond < 10^(-10)) {
					print("singular matrix detected : add diagonal noise and no more update.")
					w[i] <- 0
					#mean[[i]] <- initmean[[i]]
					cov[[i]] <- initcov[[i]]
					#cov[[i]] <- cov[[i]] + diag(d) * (10^(-10))
					noupdate <- c(noupdate, i)
					likelihood <- -Inf
				}
			}
					
		}
		
		
		w <- w / sum(w)
		
		# store nks and agitations
		if(it > 1) {
			nk[[it-1]] <- apply(resp, 2, sum)
			agitation[[it-1]] <- numeric()
			for(i in 1:k) {
				agitation[[it-1]][i] <- 0
				if(nk[[it-1]][i] > 0) {
					for(j in 1:n) {
						agitation[[it-1]][i] <- agitation[[it-1]][i] + abs(resp[j,i] - oldresp[j,i])
					}
					agitation[[it-1]][i] <- agitation[[it-1]][i] / nk[[it-1]][i]
				}
			}
		}
		
		
		newLikelihood <- getDataLikelihood(list(w=w, mean=mean, cov=cov), data)

		# Avoid Inf + Inf operation
		if((newLikelihood==-Inf) && (likelihood==-Inf)) {
			deltaLikelihood <- Inf
		} else {
			deltaLikelihood <- newLikelihood - likelihood
		}

		if(is.null(maxit)) {
			if(deltaLikelihood < thres) conv <- TRUE
		} else {
			if(it >= maxit) conv <- TRUE
		}
		likelihood <- newLikelihood
		
		
	}
	
	
	labels <- as.factor(labels)
	print(paste("final likelihood=", likelihood))
	
	for(i in 1:k) {
		mean[[i]] <- as.numeric(mean[[i]])
	}
	
	list(w=w, mean=mean, cov=cov, labels=labels, nk=nk, agitation=agitation)
		
}
			
			
			

			
			
			
		
		

	
	
