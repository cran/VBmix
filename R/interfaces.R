# Copyright (C) 2011 Pierrick Bruneau, see README for full notice

# interfaces to C entry points

varbayes <- function(data, ncomp, thres=0.1, maxit=NULL) {
	if(class(data) == "data.frame") data <- as.matrix(data)
	.Call("varbayes", data, ncomp, thres, maxit, PACKAGE="VBmix")
}

vbcomp <- function(models, ncomp, thres=0.1, maxit=NULL) {
	# in case it is not already, normalize weights of input model
	models$w <- models$w / sum(models$w)
	.Call("vbcomp", models, ncomp, thres, maxit, PACKAGE="VBmix")
}

vbconstr <- function(models, ncomp, thres=0.1, maxit=NULL) {
	models$w <- models$w / sum(models$w)
	.Call("vbconstr", models, ncomp, thres, maxit, rho=new.env(), PACKAGE="VBmix")
}

getLabels <- function(model, data) {
	# coerce means
	for(i in 1:length(model$w)) {
		model$mean[[i]] <- as.numeric(model$mean[[i]])
	}
	# coerce data
	data <- as.matrix(data)
	.Call("getLabels", model, data, PACKAGE="VBmix")
}


gmmkmsock <- function(models, names, ngroups, rho=new.env(), host="127.0.0.1") {
	.Call("gmmkmsock", models, names, ngroups, rho, host, PACKAGE="VBmix")
}


getCouple <- function(vec1, vec2) {
	.Call("couple", vec1, vec2, PACKAGE="VBmix")
}

getTimestamp <- function() {
	.Call("getTimestamp", PACKAGE="VBmix")
}

getElapsed <- function(stamp) {
	.Call("getElapsed", stamp, PACKAGE="VBmix")
}

readLabelFile <- function(name) {
	.Call("readLabelFile", name, PACKAGE="VBmix")
}

readPixmapFile <- function(name) {
	dat <- .Call("readPixmapFile", name, PACKAGE="VBmix")
	datmod <- list()
	print("post processing...")
	for(i in 1:length(dat)) {
		dat[[i]] <- (255-dat[[i]])/255
		dat[[i]] <- pixmapGrey(dat[[i]])
	}
	return(dat)
}

jsmc <- function(mod1, mod2, nsamp=5000) {
	.Call("jsmc", mod1, mod2, nsamp, PACKAGE="VBmix")
}

klmc <- function(mod1, mod2, nsamp=5000) {
	.Call("klmc", mod1, mod2, nsamp, PACKAGE="VBmix")
}


klut <- function(mod1, mod2) {
	.Call("klut", mod1, mod2, PACKAGE="VBmix")
}


jsut <- function(mod1, mod2) {
	.Call("jsut", mod1, mod2, PACKAGE="VBmix")
}



mppca <- function(data, ncomp, thres=0.1, maxit=NULL, qmax=NULL) {
	if(is.null(qmax)) qmax <- length(data[1,]) - 1
	.Call("mppca", data, ncomp, thres, maxit, qmax, PACKAGE="VBmix")
}


getResp <- function(data, model) {
	.Call("getResp", data, model, PACKAGE="VBmix")
}




mmppca <- function(mods, ncomp, thres=0.1, maxit=NULL) {
	mods$alpha <- (mods$alpha * 200000) / sum(mods$alpha)
	modred <- .Call("mmppca", mods, ncomp, thres, maxit, PACKAGE="VBmix")
	return(modred)
}


# return indexes associated to sorted values
sort_index <- function(vec, order=0) {
	# 0 = asc, 1 = desc
	.Call("sort_index", vec, order, PACKAGE="VBmix")
}

Rdct <- function(vect) {
	.Call("Rdct", vect, PACKAGE="VBmix")
	#.Call(Rdct, vect)
}

Rdct2D <- function(mat) {
	res <- .Call("Rdct2D", t(mat), PACKAGE="VBmix")
	return(t(res))
}

RinvDct2D <- function(mat) {
	.Call("RinvDct2D", t(mat), PACKAGE="VBmix")
}


extractSimpleModel <- function(model=model, labels=FALSE) {
	.Call("extractSimpleModel", model, labels, PACKAGE="VBmix")
}

dDirichlet <- function(alpha=0.1, x1, x2) {
	.Call("dDirichlet", alpha, x1, x2, PACKAGE="VBmix")
}

rDirichlet <- function(K, alpha=0.1) {
	.Call("rDirichlet", K, alpha, PACKAGE="VBmix")
}

gmmgen <- function(mod, nitem) {
	.Call("gmmgen", mod, nitem, PACKAGE="VBmix")
}

gmmdensity <- function(mod, data) {
	if(length(as.numeric(mod$mean[[1]])) > 1) {
		res <- .Call("gmmdensity", mod, data, PACKAGE="VBmix")
	} else {
		res <- rep(0, length(data))
		for(i in 1:length(mod$w)) {
			res <- res + mod$w[i] * dnorm(data, mean=mod$mean[[i]], sd=mod$cov[[i]])
		}
	}
	return(res)
}

mvngen <- function(mean, cov, nitem) {
	.Call("mvngen", mean, cov, nitem, PACKAGE="VBmix")
}

mvndensity <- function(mean, cov, data) {
	.Call("mvndensity", mean, cov, data, PACKAGE="VBmix")
}


multinomial <- function(weights, k) {
	.Call("multinomial", weights, k, PACKAGE="VBmix")
}


mixKnn <- function(data, labels, n=2, KLparam=500) {
	.Call("mixKnn", data, labels, n, KLparam, PACKAGE="VBmix")
}


mergeClassif <- function(data, labels, KLparam=500, rho=new.env()) {
	.Call("mergeClassif", data, labels, KLparam, rho, PACKAGE="VBmix")
}


constrClassif <- function(data, labels, KLparam=500, rho=new.env()) {
	.Call("constrClassif", data, labels, KLparam, rho, PACKAGE="VBmix")
}


sampleClassif <- function(data, labels, KLparam=500, rho=new.env()) {
	.Call("sampleClassif", data, labels, KLparam, rho, PACKAGE="VBmix")
}




