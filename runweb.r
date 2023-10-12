#### este programa se utiliza en GENESEARCH service 
#### http://bioinformatica.mty.itesm.mx:8080/genesearch

version <- "runweb.r v1.13 Sep/111/2015"

.errores <- list() #"**** Formato por líneas ****","Primero error y luego traceback (atrazado)."
agregaErrores <- function(...) { 
	if (length(.errores) < 5) {
		sink("errores.out")
		traceback() # traceback va atrazado
		cat(geterrmessage()) 
		sink()
		.errores[[length(.errores)+1]] <<- readLines("errores.out")
		unlink("errores.out")
	}
}

options(error=agregaErrores)

inputfile <- commandArgs(trailingOnly=TRUE)[1]
if (is.na(inputfile)) inputfile <- "lastinput.txt" 

#source("/Users/trevino/NetBeansProjects/AllScriptsApp/Scripts/R2LS.r")
source("loadR2LS.r")
if (inputfile != "lastinput.txt") file.copy(inputfile, "lastinput.txt", overwrite=TRUE)
xL <- readLines(inputfile)
#writeLines(xL, con="lastinput.txt")
len <- length(xL)
w <- grep("^data$", xL)
if (length(w) == 1) {
	l2 <- paste(xL[(w+1):len],collapse="\n")
	xL <- c(xL[1:w], l2)
	len <- length(xL)
}
if (len %% 2 == 1) {
	xL <- c(xL, "")
	len <- length(xL)
} 
params <- xL[seq(2, len, by=2)]
names(params) <- make.names(xL[seq(1, len-1, by=2)], unique=TRUE)

getParam <- function(pname, default=NULL) {
	xV <- params[pname]
	if (is.null(xV)  ||  is.na(xV)  || (is.character(xV) && length(xV) == 1 && nchar(xV) == 0))
		xV <- params[make.names(pname)]
	if (is.null(xV)  ||  is.na(xV) ||  (is.character(xV) && length(xV) == 1 && nchar(xV) == 0))
		xV <- default
	xV
}


outputfile <- function(n, ext=ifelse(isPDF,".pdf",".png")) {
	paste(gsub("\\..*", "", inputfile),n, ext, sep="")
}

as.HTML <- function(xdf, row.names=!is.null(rownames(xdf)), col.names=!is.null(colnames(xdf)), prefix="", posfix="") {
	#s[1] <- "<TABLE>"
	wordize <- function(s) gsub("([0-9])([A-Za-z])","\\1 \\2",gsub("([a-z])([A-Z0-9])","\\1 \\2",gsub("[_.]"," ",s)))
	m <- apply(xdf, 2, format)
	if (!is.matrix(m) && !is.data.frame(m)) {
		save(m, file="m.RData")
		m <- matrix(m, nrow=nrow(xdf), ncol=ncol(xdf))
	}
	colnames(m) <- if (!is.null(colnames(xdf))) colnames(xdf) else ""
	rownames(m) <- if (!is.null(rownames(xdf))) rownames(xdf) else as.character(1:nrow(m))
	if (row.names) {
		m <- cbind(rownames(m),m)
	}
	s <- character(nrow(m)+3)
	s[1] <- "<TABLE>"
	s[2] <- paste("<TR><TH>", paste(wordize(colnames(m)),collapse="</TH><TH>"),"</TD></TR>",sep="")
	s[-c(1,2,length(s))] <- apply(m,1,function(m) paste("<TR><TD>", prefix, paste(m,collapse=paste(posfix,"</TD><TD>",prefix,sep="")),posfix,"</TD></TR>",sep=""))
	#s[-c(1,2,length(s))] <- apply(m,1,function(m) paste("<TR><TD>", paste(m,collapse="</TD><TD>"),"</TD></TR>",sep=""))
	s[length(s)] <- "</TABLE>"
	return (if (col.names) s else s[-2])
}

concordance.cox <- function(cox, time, status, data, ranks=FALSE, prognostic.index=NULL) {

	if (TRUE) {
		o <- order(time)
		time <- time[o]
		status <- status[o]
		if (!is.null(data)) data <- data[,o,drop=FALSE]
		
		wst <- which(status == 1)
		ut <- sort(time[wst])
		ci <- matrix(NA, nrow=length(ut), ncol=length(time), dimnames=list(ut, time))
		if (is.null(prognostic.index)) {
			betas <- if (is.vector(cox)) cox else cox$coefficients
			if (any(is.na(betas))) betas[is.na(betas)] <- 0
			p.i <- if(!is.null(data)) betas %*% data else betas
		} else {
			p.i <- prognostic.index[o]
		}
		rpi <- rank(p.i)
		xn <- paste(1:length(time),":",time,"t",ifelse(status==1,"","+"),".",rpi, "r", sep="")
		colnames(ci) <- xn
		rownames(ci) <- xn[wst]
		#for (i in 1:(nrow(ci)-1)) {
		#	for (j in (wst[i]+1):ncol(ci)) {
		#		if (time[wst[i]] < time[j] || status[wst[i]] != status[j]) {
		#			ci[i,j] = if(ranks) rpi[j] - rpi[wst[i]] else (p.i[wst[i]] > p.i[j])
		#			# Evaluating the yield of medical tests - JAMA - Harrell
		#			# patient with higher prognostic index has higher survival time
		#		}
		#	}
		#}
		if (length(wst) >= 1) {
			for (i in 1:(nrow(ci)-1)) {
				wj <- (wst[i]+1):ncol(ci)
				wci <- which(time[wst[i]] < time[wj] | status[wst[i]] != status[wj])
				if (length(wci) > 0) {
					if (ranks) {
						ci[i, wj[wci]] <- rpi[wj[wci]] - rpi[wst[i]]
					} else {
						ci[i, wj[wci]] <- p.i[wst[i]] > p.i[wj[wci]]
					}
					# Evaluating the yield of medical tests - JAMA - Harrell
					# patient with higher prognostic index has higher survival time
				}
			}
		} else {
			ci[] <- 0
		}
	} else {
		ci <- matrix(0, nrow=3, ncol=3) ## it was used to test times
	}
	ci
}

concordance.index <- function(concordance.mat) {
	sum(concordance.mat, na.rm=TRUE) / sum(!is.na(concordance.mat))
}
somer.D.rank <- function(ci) {
	#ci is concordance index
	2*(ci-0.5)
}

concordance.risk <- function(cox, time, status, data, high.risk=time < median(time) || status == 1) {
	# This makes no sense because we don't know a priori which samples are high risk
	# Use only for special cases
	betas <- if (is.vector(cox)) cox else cox$coefficients
	o <- order(time)
	time <- time[o]
	status <- status[o]
	data <- data[,o,drop=FALSE]
	
	ut <- sort(time[status == 1])
	p.i <- betas %*% data
	rpi <- rank(p.i)

	(sum((rpi > median(rpi)) == high.risk[o]) / length(time))
}

getcutsfrompi <- function(p.i, risk.groups, all=FALSE) {
	brk <- seq(0,1,length.out=risk.groups+1)
	qcuts <- quantile(p.i, brk, na.rm=TRUE)
	if (length(qcuts) > length(unique(qcuts)) && length(qcuts) > 1) {
		opi <- sort(na.omit(p.i))
		km <- kmeans(opi, centers=risk.groups)
		qcuts <- c(min(opi), opi[diff(km$cluster) != 0], max(opi))
		brk <- c(0,which(diff(km$cluster) != 0)/length(opi),1)
		#se puede poner un algoritmo de balanced k-means 
		#o modificar el p.i artificialmente como abajo con xopi pero habria que modificar el p.i original
		#xopi <- opi+cumsum(rep(min(diff(unique(opi))),length(opi)))
		#e <- ecdf(opi)
		#for (q in 2:length(qcuts)) {
		#	if (qcuts[q] == qcuts[q-1]) {
		#	  w <- which(opi > qcuts[q-1])[1]
		#	  if (length(w) < 1) {
		#	      qcuts[q] <- max(max(p.i), qcuts+1,max(p.i))
		#	    } else {
		#	      qcuts[q] <- opi[w]
		#	      brk[q] <- e(opi[w])
		#	    }
		#	}
		#}
	}
	if (all) {
		return (list(cuts=qcuts,breaks=brk))
	} else {
		qcuts
	}
}


maximizeRiskGroups <- function(time, status, p.i, risk.groups, minpergroup=max(5,round(length(time)*(20/risk.groups)/100,0)), tr=NULL) {
	## el algoritmo es:
	## (1) Mover cada limitrofe entre el rango posible y calcular el p-value
	## (2) Establecer los nuevos rangos con el p-value menor de todos los limitrofes
	## (3) Volver a (1) si hubo cambios en los rangos
	if (is.null(tr)) tr <- 1:length(time)
	otime <- time
	ostatus <- status
	op.i <- p.i
	time <- time[tr]
	status <- status[tr]
	p.i <- p.i[tr]
	#breaks <- seq(0,1,length.out=risk.groups+1, na.rm=TRUE)
	lecuts <- getcutsfrompi(p.i, risk.groups, all=TRUE)
	qcuts <- lecuts$cuts #quantile(p.i, breaks)
	breaks <- lecuts$breaks
	ocuts <- pmax(0,round(breaks*length(p.i),0))
	if (length(unique(qcuts)) < length(qcuts)) {
		return (list(cluster.quantile=rep(1,length(p.i)), cluster.max=rep(1,length(p.i)), prisk=rep(1,length(p.i))))
	}
	pi.cluster  <- cut(p.i, breaks=qcuts, labels=FALSE, include.lowest=TRUE)
	pi.original <- pi.cluster
	minPrisk <- rep(1, length(pi.cluster))
	o <- order(p.i)
	#	prisks <- numeric(length(otime))
	#	prisks[] <- NA
	#	prisks[tr[o]] <- minPrisk
	#	return (list(cluster.quantile=pi.original, cluster.max=pi.cluster, prisk=prisks))
	steps <- max(1,length(pi.cluster) / (200*max(1,risk.groups-1)))
	if (length(pi.cluster) > minpergroup*risk.groups) {
		#ocuts contiene los rangos de pacientes ordenados por pi en clusters
		repeat {
			#cat(ocuts,"\n")
			prisk <- matrix(NA,ncol=risk.groups-1,nrow=length(pi.cluster))
			for (g in 1:(risk.groups-1)) {
				## cluster contiene los cluster actuales asignados segun los cortes en ocuts
				cluster <- numeric(length(pi.cluster))
				for (i in 1:risk.groups) {
					cluster[o[(ocuts[i]+1):ocuts[i+1]]] <- i
				}
				from <- (ocuts[g]+minpergroup)
				to <- (ocuts[g+2]-minpergroup)
				for (i in seq(from,to,length.out=(to-from+1)/steps)) { ### seq(from,to,lenght.out=(to-from)/steps) OR from:to
					kluster <- cluster
					kluster[o[from:i]] <- g
					if (i < to) kluster[o[(i+1):(to)]] <- g+1
					svd <- survdiff(Surv(time, status) ~ cluster, data=data.frame(time=time, status=status, cluster=kluster))
					prisk[i,g] <- 1-pchisq(svd$chisq,df=risk.groups-1)
				}
			}
			prm <- apply(prisk,2,min,na.rm=TRUE)
			w <- which.min(prm)
			if (length(w) > 1) w <- sample(w,1)
			minPrisk <- unlist(apply(prisk,1,min,na.rm=TRUE))
			#cat(minPrisk,"\n")
			#minPrisk <- prisk[,w]
			wna <- which(!is.finite(minPrisk))
			wnona <- which(is.finite(minPrisk))
			if (length(wna) > 0) minPrisk[wna] <- NA
			wg <- which.min(minPrisk)
			minPrisk[wna] <- minPrisk[sapply(wna, function(wi) if (any(wnona > wi)) wnona[wnona > wi][1] else max(wnona))]
			#cat(minPrisk,"\n")
			prevcuts <- ocuts
			ocuts[w+1] <- wg[1]-1
			if (all(ocuts == prevcuts) ||  risk.groups == 2) break()
		}
		pi.cluster <- numeric(length(pi.cluster))
		for (i in 1:risk.groups) {
			pi.cluster[o[(ocuts[i]+1):ocuts[i+1]]] <- i
		}
	}
	## revert to all data (not only training)
	pi.original <- cut(op.i, breaks=qcuts, labels=FALSE, include.lowest=TRUE)
	ncuts <- p.i[o[ocuts+1]]-abs(min(diff(sort(p.i))))/10
	ncuts[1] <- -Inf
	ncuts[length(ncuts)] <- +Inf
	if (length(ncuts) > length(unique(ncuts))) {
		for (q in 2:length(ncuts)) {
			if (ncuts[q] == ncuts[q-1]) {
			  w <- which(p.i > ncuts[q-1])[1]
			  if (length(w) < 1) {
			      ncuts[q] <- max(ncuts+1,max(p.i))
			    } else {
			      ncuts[q] <- p.i[w]
			    }
			}
		}
	}
	pi.cluster <- cut(op.i, breaks=ncuts, labels=FALSE, include.lowest=TRUE)
	prisks <- numeric(length(otime))
	prisks[] <- NA
	prisks[tr[o]] <- minPrisk
	return (list(cluster.quantile=pi.original, cluster.max=pi.cluster, prisk=prisks))
}



cox.log.likelihood <- function(time, status, betas, data) {
	## l(b) <- sum(statusi*(xi'*b-log(sum(exp(xj'*b)))
	o <- order(time)
	n <- length(o)
	time <- time[o]
	status <- status[o]
	data <- data[,o,drop=FALSE]
	t1 <- betas %*% data
	ll <- 0
	for (i in 1:n) {
		ll <- ll + status[i] * (t1[i] - log(sum(exp(betas %*% data[,i:n,drop=FALSE]))))
	}
	ll
}

cox.r2.max <- function(time, status, betas, data) {
	llb <- cox.log.likelihood(time, status, betas-betas, data)
	return (1-exp(2*llb/ncol(data)))
}

cox.r2 <- function(time, status, betas, data) {
	dev <- -2 * (cox.log.likelihood(time, status, betas, data) - 
				cox.log.likelihood(time, status, betas-betas, data))
	return (1-exp(dev/ncol(data)))
}

plot.raw.kaplan.per.strata <- function(time, status, strata, 
	col=1:nlevels(factor(strata))+1, conf.int=FALSE, main="",
	draw.main=FALSE) {
	strata <- factor(strata)
	nstrata <- nlevels(strata)
	pp = par(mar=c(5.1+max(0,nstrata-1), 4.1, 5.1, 2.1))
	on.exit(par(pp))
	surdif <- survdiff(Surv(time, status) ~ cluster, data=data.frame(time=time, status=status, cluster=strata))
	p <- 1-pchisq(surdif$chisq,df=nstrata-1)
	hr <- coxph(Surv(time, status) ~ cluster, data=data.frame(time=time, status=status, cluster=strata), method="breslow") 
	shr <- summary(hr)
	#print("******* plot.raw.kaplan.per.strata ******")
	#print(main)
	#print(hr)
	#print(shr)
	#print(str(shr))
	plot(survfit(Surv(time, status) ~ cluster, 
		data=data.frame(time=time, status=status, cluster=strata)), 
		col=col, mark.time=TRUE, pch=mark.char,
		lty=1, lwd=2, 
		conf.int=conf.int,
		main=paste(main,"\nLog-Rank Equal Curves p=",format(p, digits=4),
			ifelse(nstrata == 2, paste("\nRisk Groups Hazard Ratio = ",format(round(shr$conf.int[1,1],2),digits=4)," (conf. int. ",format(round(shr$conf.int[1,3],2),digits=4)," ~ ",format(round(shr$conf.int[1,4],2),digits=4),"), p=",format(shr$coefficients[1,5], digits=4), sep=""), ""), sep=""),
		cex.main=1.25)
	toview = 1:nstrata
	if (draw.main) {
		lines(survfit(Surv(time, status) ~ 1), mark.time=TRUE, col=16,
		 	conf.int=conf.int, lty=1, lwd=1, pch=mark.char)
		toview <- c(toview,nstrata+1)
	}
	xt <- axTicks(1)
	xt <- sort(c(xt, xt[-length(xt)]+diff(xt)/2))
	censxrisk <- numeric(nstrata)
	for (k in 1:nstrata) {
		kl <- levels(strata)[k]
		w <- which(strata == kl)
		censxrisk[k] <- sum(status[w] == 0)
		axis(1, xt, labels=sapply(xt, function(x) sum(time[w] >= x)), tick=FALSE, line=k, col.axis=col[k])
	}
	legend("topright", 
		paste(c(levels(strata),"(All)"),":",c(table(strata),length(time)),", +:",c(censxrisk,sum(status==0)),sep="")[toview], 
		col=c(col[1:nstrata],8)[toview], 
		ncol=1, lwd=2, bg=0)
}

plot.kaplan <- function(data, time, status, model, tr=1:ncol(data), risk.groups=2, 
	main="", tr.betas=tr, col=1:risk.groups+1, roundcoef=3, symbols=NULL, 
	betas=NULL, betas.p=rep(1, length(betas)), draw.text=TRUE, draw.main=TRUE, 
	needMaximizeGroups=FALSE, univariateBetas=FALSE,
	col2=1:risk.groups+1, draw=TRUE) {
	#col[1:risk.groups] <- col[risk.groups:1]
	library(survival)
	ci <- NA
	betas.supplied <- ! is.null(betas)
	coxp <- betas.p
	estimated <- FALSE
	cox <- NULL
	try (
		cox <- coxph(Surv(time[tr.betas], status[tr.betas]) ~ .,data.frame(t(data[model,tr.betas,drop=FALSE])), method="breslow")
		)
	if (is.null(betas)) {
		estimated <- is.null(cox)
		if (estimated  ||  univariateBetas) {
			betas <- c()
			coxp <- c()
			for (i in 1:length(model)) {
				cox <- coxph(Surv(time[tr.betas], status[tr.betas]) ~ x,data.frame(x=data[model[i],tr.betas]), method="breslow")
				betas <- c(betas, cox$coefficients)
				coxp <- c(coxp, summary(cox)$coefficients[,5])
			}
		} else {
			betas <- cox$coefficients
			coxp <- summary(cox)$coefficients[,5]
		}
	}
	coxo <- cox
	scoxo <- summary(coxo)
	if (any(is.na(betas))) 
		betas[is.na(betas)] <- 0
    ci <- concordance.index(concordance.cox(betas, time[tr], status[tr], data[model, tr, drop=FALSE]))
	p.i <- betas %*% data[model, , drop=FALSE]

	maxRisk <- toupper(getParam("max_risk", "YES")) %in% c("1","YES","TRUE")
	if (maxRisk || needMaximizeGroups) {
		mrg <- maximizeRiskGroups(time, status, p.i, risk.groups, tr=tr.betas)
		pi.cluster <- if (maxRisk) mrg$cluster.max else mrg$cluster.quantile
	} else {
		qcuts <- quantile(p.i, seq(0,1,length.out=risk.groups+1, na.rm=TRUE))
		pi.cluster <- cut(p.i, breaks=qcuts, labels=FALSE, include.lowest=TRUE)
		mrg <- NULL
	}

	### hazard-ratio estimation of the PI as predictor, THIS WILL ALWAYS BE 1
	### see Assessment of performance of survival prediction models for cancer prognosis - Chen - Chen - BMC Medical Research Methodology 2012.pdf
	### Instead, it is better to estimate the hazard ratio by the two-groups
	### That is, fitting a cox model to the cluster as predictor
	#hr <- coxph(Surv(time[tr.betas], status[tr.betas]) ~ PI, data.frame(PI=p.i[tr.betas]), method="breslow") 
	hr <- coxph(Surv(time, status) ~ cluster, data=data.frame(time=time[tr], status=status[tr], cluster=pi.cluster[tr]), method="breslow") 
	shr <- summary(hr)


	surdif <- survdiff(Surv(time, status) ~ cluster, data=data.frame(time=time[tr], status=status[tr], cluster=pi.cluster[tr]))
	p <- 1-pchisq(surdif$chisq,df=risk.groups-1)
	
	vars <- if (is.null(rownames(data))) model else rownames(data)[model]
	xm <- data.frame("Symbol:Id"=vars)
	xm$Beta <- round(betas,roundcoef)
	xm$pValue <- coxp
	xcol <- matrix(1, ncol=risk.groups*2+3, nrow=length(model)+1)
	xcol[-1,2:3] <- ifelse(coxp < 0.05, 16, ifelse(coxp < 0.1, 8, 1))
	xcol[ 1,2:3] <- 16
	res <- list()
	cox <- list()
	cox.ci <- list()
	cox.cirefit <- list()
	censxrisk <- rep(0,risk.groups)
	survcox <- list()
	xnames <- c()
	for (k in 1:risk.groups) {
	    xcol[1,2:3+k*2] <- col[k]
		w <- which(pi.cluster[tr] == k)
		censxrisk[k] <- sum(status[tr[w]] == 0)
		coxk <- NULL
		xnames <- c(xnames, paste(c("Beta","pValue"),k,sep=""))
		#if (draw.text) {
			try(coxk <- coxph(Surv(time[tr[w]], status[tr[w]]) ~ ., data.frame(t(data[model, tr[w], drop=FALSE])), method="breslow"))
		#}
		cox.ci[[k]] <- cox.cirefit[[k]] <- NA
		if(!is.null(coxk)) {
		  cox[[k]] <- coxk 
		  survcox[[k]] <- NULL
		  try(survcox[[k]] <- survfit(cox[[k]]))
		  res[[k]] <- data.frame(coef=round(cox[[k]]$coefficients,roundcoef), 
		  							p=summary(cox[[k]])$coefficients[,5])
		  xm <- cbind(xm, res[[k]])
		  xcol[-1,2:3+k*2] <- ifelse(res[[k]]$p < 0.05, col[k], ifelse(res[[k]]$p < 0.1, rgb2rgb(col2rgb(col[k])/2), 1))
		  ## con modelo del refit
		  try (cox.cirefit[[k]] <- concordance.index(concordance.cox(coxk, time[tr[w]], status[tr[w]], data[model, tr[w], drop=FALSE])) )
		} else {
		  cox[[k]] <- NA
		  xm <- cbind(xm, matrix(NA, ncol=2, nrow=length(model)))
		}
	  ## con modelo "original" (de todo el grupo)
	  try (cox.ci[[k]] <- concordance.index(concordance.cox(coxo, time[tr[w]], status[tr[w]], data[model, tr[w], drop=FALSE])) )
	}
	colnames(xm) <- c(colnames(xm)[1:3], xnames)
	
	#r2 <- scoxo$rsq["rsq"] 
	#r2max <- scoxo$rsq["maxrsq"]
	## if betas.supplied then tr is supposed to be test in this case, so r2 is about the correlation in test
	## these 2 functions made the estimation in any case
	r2 <- cox.r2(time[tr], status[tr], betas, data[model, tr, drop=FALSE])
	r2max <- cox.r2.max(time[tr], status[tr], betas, data[model, tr, drop=FALSE])

	pp = par(mar=c(5.1+max(0,risk.groups-2), 4.1, 5.1, 2.1))
	on.exit(par(pp))
	if (draw && length(unique(pi.cluster[tr])) > 1) plot(survfit(Surv(time, status) ~ cluster, 
			data=data.frame(time=time[tr], status=status[tr], cluster=pi.cluster[tr])), 
		col=col,
		lty=ifelse(col == col2, 1, 1:groups), lwd=2, 
		conf.int=FALSE,
		main=paste(main,"\nConcordance Index = ",round(ci*100,2),", Log-Rank Equal Curves p=",format(p, digits=4), 
			", R^2=",format(round(r2,3)),"/",format(round(r2max,3)),
			ifelse(risk.groups==2,paste("\nRisk Groups Hazard Ratio = ",format(round(shr$conf.int[1,1],2),digits=4)," (conf. int. ",format(round(shr$conf.int[1,3],2),digits=4)," ~ ",format(round(shr$conf.int[1,4],2),digits=4),"), p=",format(shr$coefficients[1,5], digits=4),sep=""),""),
			ifelse(estimated,"\n(Full Cox model invalid:using univariate Cox)",""),sep=""), 
		ylim=c(ifelse(draw.text, -0.15,0), 1.05),
		cex.main=1.25, mark.time=TRUE, pch=mark.char)

	#plot(0,0,type="n",
	#	main=paste(main,"\nci=",round(ci*100,2),"|p=",format(p),
	#		ifelse(estimated,"\n*** FULL COX MODEL INVALID ***",""),sep=""), 
	#	ylim=c(ifelse(draw.text, -0.15,0), 1),
	#	xlim=c(1,max(time[tr])))
	#for (k in 1:risk.groups) {
	#	ok <- which(pi.cluster == k)
	#	lines(survfit(Surv(time, status) ~ cluster, 
	#		data=data.frame(time=time[tr[ok]], status=status[tr[ok]], cluster=pi.cluster[ok])), 
	#		col=col[k], 
	#		lty=1, lwd=2, 
	#		conf.int=FALSE)
	#}
	toview <- 1:(risk.groups+ifelse(draw.main,1,0))
	#pic.cens <- paste(pi.cluster,status)[tr[status[tr]==0]]
	#table(pic.cens)
	if (draw) legend("topright", 
		paste(paste(c(table(pi.cluster[tr]),length(tr)),", +:",c(censxrisk,sum(status[tr]==0)),sep=""), ", CI=",round(c(unlist(cox.ci),ci)*100,1),if (draw.text) paste(", CIrefit=",round(c(unlist(cox.cirefit),ci)*100,1),sep="") else rep("",risk.groups+1),sep="")[toview], 
		col=c(col[1:risk.groups],8)[toview], 
		lty=ifelse(col == col2, 1, 1:groups), 
		ncol=round((risk.groups+ifelse(draw.main,1,0))/2.5,0), lwd=2, bg=0)
	pos <- c("bottomleft", "bottomright", "top", "right")
	xt <- axTicks(1)
	xt <- sort(c(xt, xt[-length(xt)]+diff(xt)/2))
	#xx=par("usr")[2]*min(1,(3+2*risk.groups)/9)
	if (!is.null(symbols)) {
		#xm <- cbind(data.frame(Symbol=symbols),xm)
		#xcol <- cbind(rep(1,nrow(xcol)),xcol)
	}
	if (draw.text) {
		plot.text(xm, col=xcol, add=TRUE, x0=0*par("usr")[1], xx=par("usr")[2],
		 	y0=min(1,nrow(xm)/20-0.1), yy=(-0.15+par("usr")[3])/2, 
		 	formatter=function(x) if (is.numeric(x)) format(x, digits=5) else format(x, digits=5),
		 	vlines=0, hlines=0)
	}
	#lines(survfit(maincox), mark.time=TRUE, col=8, conf.int=FALSE, lty=1, lwd=2)
	if (draw.main) {
		lines(survfit(Surv(time[tr], status[tr]) ~ 1), col=16,
		 	conf.int=FALSE, lty=1, lwd=1, mark.time=TRUE, pch=mark.char)
	}
	for (k in 1:risk.groups) {
		#if (!is.na(cox[[k]])) {
			#lines(survcox[[k]], mark.time=TRUE, col=k+1, conf.int=FALSE, lty=3, lwd=1)
		#}
		w <- which(pi.cluster[tr] == k)
		if (draw) axis(1, xt, labels=sapply(xt, function(x) sum(time[tr[w]] >= x)), tick=FALSE, line=k, col.axis=col[k])
	}

	result <- list(surdif, pi=p.i, pi.cluster=pi.cluster, betas=betas, fitting=xm, ticks=xt, 
		cox=coxo, coxes=cox, 
		hr=shr$conf.int[1,1], hrl=shr$conf.int[1,3], hru=shr$conf.int[1,4], hr.p=shr$coefficients[1,5],
		rsq=r2, rsqmax=r2max, ci=ci, log.rank=p)
	if (!is.null(mrg)) {
		result$p.risk <- mrg$prisk
		result$pi.cluster.original <- mrg$cluster.quantile		
	}
	result
}

plots_made <- c()

new_plot <- function(filename, ...) {
	close_plot()
	cat("Starting Plot:",filename,"\n")
	if (isPDF) {
		if (length(plots_made) == 0) {
			## filename <- paste(filename,".pdf",sep="")
			pdf(file=filename, width=ceiling(width / (72)), height=ceiling(height / (72)), bg=bg)
			plots_made <<- filename
		}
	} else {
		try(dev.off())
		png(filename=filename, width=width, height=height,bg=bg)
		plots_made <<- c(plots_made, filename)
	}
}

end_plots <- function(filename, ...) {
	close_plot()
	try(dev.off())
}

close_plot <- function(...) {
	if (length(plots_made) > 0) {
		cat("Ending Plot:",plots_made[length(plots_made)],"\n")
	} else {
		cat("No plot to end.\n")
	}
  if (dev.cur() > 1 & !isPDF) {
    try(mtext(paste(getParam("database",""), sep=""), cex=0.85, line=-1, side=1, adj=0, outer=TRUE, col=6))
    try(mtext(paste("(",version, ") ", format(Sys.time(), "%Y-%m-%d %H:%M"), sep=""), cex=0.666, line=-1, side=1, adj=1, outer=TRUE, col=6))
  }
}

################################## MAIN ##################################
xfunc <- getParam("function.", getParam("function", ""))
if (is.null(xfunc)) {
	write.table(file="salida.txt", "Null function")
	cat("salida.txt")
	q("no")
}

isPDF <- (getParam("pdf", "0") != "0")

outfile <- paste(gsub("\\..*", "", inputfile), ifelse(isPDF,".pdf",".png"), sep="")

mark.char <- getParam("char_mark","20")
if (is.finite(as.numeric(mark.char))) mark.char <- as.numeric(mark.char)

sep <- getParam("separator", "\t")
if (sep == "tab") sep = "\t"
xdata <- strsplit(getParam("data", ""), "\n|;;;")[[1]]

if (length(xdata) > 0) {
	xdata <- xdata[which(unlist(lapply(xdata, nchar)) > 10)]
	alldata <- sapply(xdata, strsplit, sep, fixed=TRUE, USE.NAMES=FALSE)
	#evitar cambios de tamaño de lineaas
	sizes <- unlist(lapply(alldata, length))
	ml <- max(sizes, na.rm=TRUE)
	if (any(sizes != ml)) {
		for (i in which(sizes != ml)) {
			alldata[[i]] <- c(alldata[[i]], rep("",ml-length(alldata[[i]])))
		}
	}
	alldata <- t(data.frame(alldata))
} else {
	alldata <- matrix(0, ncol=10000, nrow=10)
}
rownames(alldata) <- colnames(alldata) <- NULL
columns <- unlist(strsplit(getParam("columns", ""), ""))
data <- data.matrix(alldata[,which(columns=="D"),drop=FALSE])
mode(data) <- "numeric"


wsymb <- which(columns=="S")
symb <- character(nrow(data))
if (length(wsymb) > 0) {
	xs <- alldata[,wsymb,drop=FALSE]
	for (i in 1:ncol(xs)) {
		xs[,i] <- gsub(",.+","",xs[,i])
	}
	symb <- unlist(apply(xs,1,paste,collapse=":"))
}
wid <- which(columns=="I")
ids <- character(nrow(data))
if (length(wid) > 0) {
	ids <- gsub(",.+","",as.character(alldata[,wid[1],drop=FALSE]))
}
rownames(data) <- paste(symb,ids,sep=":")

width <- as.numeric(getParam("width",480))
height <- as.numeric(getParam("height",480))
bg <- getParam("bg","transparent")

if (getParam("standardize","0") != "0") {
	data <- t(scale(t(data), center=TRUE, scale=TRUE))
}

if (any(is.na(data))  &&  getParam("naimputation","none") != "none") {
	na <- getParam("naimputation","none")
	if (na == "KNN") {
		data <- na.impute(data, method="dynamic")
	} else if (na == "0") {
		data[is.na(data)] <- 0
	} else if (na == "row") {
		for (i in 1:nrow(data)) {
			if (any(is.na(data[i,,drop=FALSE]))) {
				data[i,is.na(data[i,]),drop=FALSE] <- mean(data[i,,drop=FALSE], na.rm=TRUE)
			}
		}
	} else if (na == "col") {
		for (i in 1:ncol(data)) {
			if (any(is.na(data[,i]))) {
				data[is.na(data[,i]), i, drop=FALSE] <- mean(data[,i,drop=FALSE], na.rm=TRUE)
			}
		}
	} else if (na == "mean") {
		data[is.na(data)] <- mean(data, na.rm=TRUE)
	} else if (na == "min") {
		data[is.na(data)] <- min(data, na.rm=TRUE)
	}
}


if (getParam("quantize","0") != "0") {
	quantos <- as.numeric(getParam("quantize",2))
	rg <- round(quantile(data,c(.025,.975)),0) ## range(data)
	if (rg[2] == rg[1]) rg[2] = rg[2] + 1
	brk <- seq(rg[1],rg[2],length.out=quantos+1)
	brk[1] <- min(data)
	brk[length(brk)] <- max(data)
	drn <- rownames(data)
	data <- t(apply(t(data), 2, function(x) cut(x, labels=FALSE, breaks=brk, include.lowest=TRUE)))
	if (any(is.na(data))) {
		data[is.na(data)] <- quantos
	}
	data <- matrix(seq(0,1,length.out=quantos)[data],ncol=ncol(data))
	rownames(data) <- drn
}

paramsplit <- strsplit(params, "\t")
inputParams <- params[which(unlist(lapply(paramsplit,length)) == 1)]
w <- which(names(inputParams) == "columns")
if (length(w) > 0) inputParams <- inputParams[-w]
logFile <- "/biomatecdatasets/datasets/runweb-log.txt"
inputParams <- c("Symbols"=paste(symb, collapse=";"), input=paste(rownames(data), collapse=";"), inputFile=inputfile, processTime=date(), inputParams)
if (file.exists(logFile)) {
	xf <- file(logFile, open="at")
	writeLines(substr(paste(names(inputParams),"=",inputParams),1,80), xf)
	close(xf)
}
if (isPDF) {
	new_plot(filename=outputfile("PAR"))
	plot.text(data.frame(Parameter=c(names(inputParams),"version"),Value=c(inputParams, version)), max.chars=80)
	new_plot(filename=outputfile("GENES"))
	plot.text(data.frame(Num=1:nrow(data),Symbol=symb,Id=ids), max.chars=80)
}

## build headers
headers <- NULL
for (i in 1:1000) {
	h <- getParam(paste("header",i,sep=""))
	if (length(h) == 0 || nchar(h) == 0) break;
	headers <- c(headers, h)
}
if (!is.null(headers)) {
	xh <- headers[which(nchar(headers) > 10)]
	allh <- sapply(xh, strsplit, sep, fixed=TRUE, USE.NAMES=FALSE)
	#evitar cambios de tamaño de lineaas
	sizes <- unlist(lapply(allh, length))
	ml <- max(sizes)
	if (any(sizes != ml)) {
		for (i in which(sizes != ml)) {
			allh[[i]] <- c(allh[[i]], rep("",ml-length(allh[[i]])))
		}
	}
	allh <- t(data.frame(allh, stringsAsFactors=FALSE))
	rownames(allh) <- colnames(allh) <- NULL
	headers <- t(allh)
	colnames(headers)[1:min(ncol(headers),ncol(data))] <- colnames(data)[1:min(ncol(headers),ncol(data))]
	w <- grep("#ORDINAL:|#CENSORED:|#CLASS:|#STRATA:",toupper(headers[1,]))
	if (length(w) == 0) {
		headers <- NULL
	} else {
		hc <- headers[,ncol(headers),drop=FALSE]
		colnames(data) <- make.names(gsub("#|DATA:","",toupper(hc[columns=="D"])))
		hc <- headers[1,w,drop=FALSE]
		headers <- headers[,w,drop=FALSE]
		colnames(headers) <- make.names(gsub("#","",hc))
		headers <- headers[columns=="D",,drop=FALSE]
	}
	#rownames(allh) <- make.names(gsub("#","",allh[,1]))
}

################ FUNCTIONS

if (xfunc == "heatmap") {
	xmargins <- as.numeric(unlist(strsplit(getParam("margins","5,7"),",")))
	xcolor <- -as.numeric(getParam("colors", 1))
	if (xcolor >= 0) xcolor = -1
	norm <- as.numeric(getParam("normalize", 1))
	xnames <- gsub("DATA\\.","",getParam("names",colnames(all))[columns=="D"])
	if (nrow(data) < 2) data <- rbind(data, data)
	if (ncol(data) < 2) data <- cbind(data, data)
	predictor <- getParam("predictor", NULL)
	if (is.null(predictor) || is.na(predictor)  || 
		is.null(getParam(predictor))  ||  is.na(getParam(predictor))) {
		predictor <- NULL
	} else {
		predictor <- data.frame(Predictor=unlist(strsplit(getParam(predictor),sep))[columns=="D"])
	}
	colnames(data) <- xnames
	#predictor <- NULL
	xcluster <- as.numeric(getParam("clusterized", 0))
	ncols <- as.numeric(getParam("colorlen", 24))
	new_plot(filename=outfile)
	plot.heatmap.info(data, scale=ifelse(norm==1, "row", ifelse(norm==2,"col","none")), col=redgreenblue(xcolor, len=ncols), margins=xmargins, ColSideFactors=predictor, ColvFactor=if (xcluster != 0 && !is.null(predictor)) 1 else NULL)
}

if (xfunc == "cna-survival-groups") {
	predictor = getParam("censored", "no predictor")
	if (length(predictor) == 1 && predictor == "no predictor") {
		new_plot(filename=outputfile("err"))
		plot.text(data.frame(Error=c("No predictor was specified.", "Please select a CENSORED option.")))
		close_plot()
		xfunc = "error"
	} else {
		if (getParam("riskgeneration","cox") == "cox") {
			xfunc = "cox-survival-groups"
		}
	}
}

if (xfunc == "cox-survival-groups") {
	predictor = getParam("censored", "no predictor")
	if (length(predictor) == 1 && predictor == "no predictor") {
		new_plot(filename=outputfile("err"))
		plot.text(data.frame(Error=c("No predictor was specified.", "Please select a CENSORED option.")))
		close_plot()
		xfunc = "error"
	}
}
if (xfunc == "cox-survival-groups") {
	xmargins <- as.numeric(unlist(strsplit(getParam("margins","5,20"),",")))
	groups <- as.numeric(getParam("groups", 2))
	xcolor <- rgb2rgb(col2rgb(strsplit(getParam("colors", paste(as.character(1:groups+1), collapse=",")),",")[[1]]))
	if (!any(is.na(as.numeric(xcolor)))) {
		xcolor <- as.numeric(xcolor)
	}
	xcolor <- rev(xcolor[1:groups])
	survcolor <- xcolor
	if (length(xcolor) < 2 || xcolor[1] == xcolor[2]) {
		xcolor <- gray(0:(groups-1) / (groups*2-2))
	}
	hmcol <- -as.numeric(getParam("hmcol", 1))
	#main.fit <- toupper(getParam("show_main","0")) %in% c("YES","1","TRUE")
	advanced <- toupper(getParam("advanced","0")) %in% c("YES","1","TRUE")
	fitting.info <- getParam("include_fitting_info","0") != "0"
	censored <- getParam(predictor,"")
	censored <- strsplit(censored, sep)[[1]][which(columns=="D")]
	time <- as.numeric(gsub("\\+| ","",censored))
	status <- rep(1, length(time))
	status[grep("+",censored,fixed=TRUE)] <- 0
	ok <- which(!is.na(time) & time >= 0)
	data <- data[,ok,drop=FALSE]
	time <- time[ok]
	status <- status[ok]
	if (!is.null(headers)) {
		headers <- headers[ok,,drop=FALSE]
	}

	ofactors <- getParam("otherfactors","")
	if (nchar(ofactors) > 0  &&  !is.null(headers)) {
		headnames <-  gsub("^[^.]+\\.","",colnames(headers))
		#Ejemplos:
		#N.STAGE
		#Agregar N.STAGE a los datos para la regresión, si es numérica la deja tal cual, si no, se convierte a factor y se hace la regresión
		#N.STAGE+T.STAGE
		#Agregar N.STAGE y T.STAGE a la regresion
		#N.STAGE+T.STAGE/I:I,IA,IB/II:II,IIA,IIB,IIC/III:III,3A,3B,4A,4B
		#Agregar N.STAGE tal cual y T.STAGE convertido a los grupos I, II, III
		#AGE/KID:<12/TEEN:<=18/ADULT:>18
		#NA are kept as NA or converted if NA is included in the list.
		#@T.STAGE/I:I,IA,IB/II:II,IIA,IIB,IIC/III:III,3A,3B,4A,4B
		#**** En este ejemplo, TSTAGE es transformado PERO SOLO EN LOS HEADERS esto se usa para cuando se hace STRATIFICATION
		of <- trim(unlist(strsplit(ofactors,"\\+")))
		for (i in 1:length(of)) {
			ofn <- trim(gsub("/.+$","",of[i]))
			isrep <- (length(grep("@",ofn)) > 0)
			ofn <- gsub("@","",ofn)
			facname <- ofn
			g <- grep("/",of[i])
			if (ofn %in% headnames) {
				v <- as.character(headers[,which(headnames == facname)])
				n <- as.numeric(v)
				nas <- sum(is.na(n))
				conv <- TRUE
				if (nas < length(v)*0.33) {
					# dejar numérico
					if (length(g) > 0) {
						# numérico a factor
						vv <- as.character(v)
						vv[] <- NA
						xvals <- unlist(strsplit(trim(substring(of[i],nchar(ofn)+2)),"/"))
						for (k in 1:length(xvals)) {
							xele <- unlist(strsplit(xvals[k],":"))
							tx <- paste("n",xele[2],sep="")
							xw <- which(eval(parse(text=tx)) & is.na(vv))
							if (length(xw) > 0) {
								vv[xw] <- xele[1]
							}
						}
						v <- vv
						g <- c()
					} else {
						# numérico as is
						v <- matrix(n, ncol=ncol(data), nrow=1)
						conv <- FALSE
					}
				}
				if (conv) {
					# caracter=factor
					if (length(g) > 0) {
						# convertir
						vv <- v
						vv[] <- NA
						xvals <- unlist(strsplit(trim(substring(of[i],nchar(ofn)+2)),"/"))
						for (k in 1:length(xvals)) {
							xele <- unlist(strsplit(xvals[k],":"))
							xv <- unlist(strsplit(xele[2],","))
							for (j in 1:length(xv)) {
								xw <- which(v == xv[j])
								vv[xw] <- xele[1]
							}
						}
						v <- vv
					} else {
						# as is
					}
					u <- unique(trim(v))
					vv <- matrix(0, ncol=ncol(data), nrow=length(u)) # dummy matrix
					nn <- c()
					for (k in 1:length(u)) {
						w <- which(v == u[k] | (is.na(v) & is.na(u[k])))
						vv[k,w] <- 1
						nn <- c(nn, paste(ofn,u[k],sep="."))
					}
					v <- vv
					ofn <- nn
					u <- apply(vv,1,function(x) length(unique(x)) == 1)
					if (any(u)) {
						vv <- vv[-which(u),]
						ofn <- ofn[-which(u)]
					}
				}
				if (!isrep) {
					rn <- rownames(data)
					data <- rbind(data, v)
					symb <- c(symb, ofn)
					ofn <- gsub("\\+","Pos",ofn)
					ofn <- gsub("\\-","Neg",ofn)
					rownames(data) <- c(rn, ofn)
				} else {
					if (conv) {
						ofn <- gsub(paste(facname,"\\.?",sep=""),"",ofn)
						v <- apply(v,2,function(x) ofn[which(x == 1)])
					}
					headers[,which(headnames == facname)] <- v
				}
			}
		}
	}

	strataPar <- getParam("data_groups","<SELECT>")
	maxRisk <- toupper(getParam("max_risk", "YES")) %in% c("1","YES","TRUE")
	if (strataPar != "<SELECT>") {
		separated <- as.numeric(getParam("separated","0"))
		strata <- strsplit(getParam(strataPar), sep)[[1]][which(columns=="D")]
		strata <- strata[ok]
		if (make.names(strataPar) %in% colnames(headers)  &&  nchar(ofactors) > 0) {
			# posible transformación, usar el dato de headers transformado
			strata <- headers[,which(make.names(strataPar) == colnames(headers))]
		}
		## 
	}
	tr.betas <- eval(parse(text=paste("c(",getParam("trainset", paste("1:",ncol(data),sep="")),")")))
	tr <- tr.betas
	test <- eval(parse(text=paste("c(",getParam("testset", paste("1:",ncol(data),sep="")),")")))
	train.is.dif <- (length(setdiff(tr,test)) != 0)
	weights <- as.numeric(unlist(strsplit(getParam("weights","0"),"[^0-9eE+.\\-]+")))
	useWeights <- any(weights != 0)
	if (train.is.dif || useWeights) {
		new_plot(filename=outputfile("ST"))
		plot.kaplan(data[,tr,drop=FALSE], time=time[tr], status=status[tr], model=1:nrow(data), risk.groups=groups, 
			col=xcolor, col2=survcolor,
			main=paste(paste(predictor,ifelse(train.is.dif,"(train only)",""),ifelse(useWeights,"(not using weights)","")),paste(symb,collapse="|"),sep="\n"), 
			symbols=symb,
			draw.text=fitting.info,
			draw.main=fitting.info)
	}
	outfile <- outputfile("A")
	new_plot(filename=outfile)
	km <- plot.kaplan(data, time=time, status=status, model=1:nrow(data), risk.groups=groups, 
		col=xcolor, col2=survcolor,
		main=paste(paste(predictor,ifelse(train.is.dif," (test only)",""),ifelse(useWeights," (Using weights)",""),sep=""),paste(symb,collapse="|"),sep="\n"), 
		symbols=symb, tr=test, tr.betas=tr.betas,
		betas=if (useWeights) rep(weights, length.out=nrow(data))[1:nrow(data)] else NULL,
		draw.text=fitting.info,
		draw.main=fitting.info,
		needMaximizeGroups=TRUE)

	fittingfile <- outputfile("S",ext=".html")
	sink(file=fittingfile, type="output")
	if (!isPDF) cat("<TABLE><TR><TD>")
	print(summary(km$cox))
	if (!isPDF) cat("</TD></TR></TABLE>")
	sink()
	#plots_made <- c(plots_made, fittingfile)

	allPI <- as.vector(km$betas %*% data[, , drop=FALSE])
	
	##
	##new_plot(filename=outputfile("HR"))
	##plot(rank(time), rank(-allPI), pch=ifelse(status!=0,1,3), col=xcolor[max(km$pi.cluster)-km$pi.cluster+1],
	##	xlab="Ranked Time", ylab="Inverse Rank: Prognostic Index (Risk Score)", 
	##	main="Risk Score vs Time (+=censored, o=no censored)\n.") #, type="n"
	##axis(3,at=quantile(rank(-allPI),c(.2,.8)),c("<-- Higher-Risk","Lower-Risk -->"))
	##axis(4,at=quantile(rank(time),c(.15,.85)),c("<-- Higher-Risk","Lower-Risk -->"))
	##abline(h=which(diff(km$pi.cluster[order(-allPI)]) != 0), lty=2)
	##abline(v=which(diff(km$pi.cluster[order(-allPI)]) != 0), lty=2)
	#text(rank(time), rank(-allPI), 1:length(time),  col=xcolor[ifelse(status!=0,1,2)])

	risks <- paste(c("Low","Low/Med","Medium","Med/High","High"),"Risk")
	if (groups == 2) {
		risks <- risks[c(1,5)]
	} else if (groups == 3) {
		risks <- risks[c(1,3,5)]
	} else if (groups == 4) {
		risks <- risks[-3]
	} else if (groups > 5) {
		risks <- 1:groups
	}
	if (is.null(colnames(data))) colnames(data) <- rep("",ncol(data))

	if (groups <= 5) {
		risks <- paste(groups:1,risks,sep="-")
	}
	hm <- as.numeric(getParam("heatmap",1))
	if (hm > 0) {
		norm <- as.numeric(getParam("normalize", 0))
		ncols <- as.numeric(getParam("hm_colorlen", 24))
		new_plot(filename=outputfile("H"))
		xdf <- data.frame(
				Risk=risks[km$pi.cluster],
				Censored=ifelse(status==0,"Yes","No"),
				Time=round(time,2),
				Prog.Idx=round(as.vector(km$pi),3)
				)
		forced = NULL
		if (strataPar != "<SELECT>") {
			xdf <- cbind(xdf, Class=strata)
		} else if (hm %in% c(6,7)) {
			hm <- 0
			forced <- "STRATA NOT SELECTED"
		}
		if (train.is.dif) {
			xdf <- cbind(xdf, Set=ifelse(1:ncol(data) %in% tr,"Train","Test"))
		}
		xdata <- data
		if (nrow(xdata) == 1) {
			xdata <- rbind(xdata, rep(mean(xdata), ncol(xdata)))
		}
		hmi <- plot.heatmap.info(xdata[,,drop=FALSE], 
			scale=ifelse(norm==1, "row", ifelse(norm==2,"col","none")), 
			ColSideFactors=xdf,
			ColSideFactors.type=c("i","i","h","l","i","i"),
			ColSideFactors.col=list(Risk=xcolor[groups:1], Censored=c(8,1)),
			ColvFactor=if (hm == 2 || hm == 3  ||  hm == 6) hm-1 else NULL,
			ColSideFactors.order.function = if (hm==8 || hm==9) ifelse(hm==8, order.by.ranks, order.by.correlations) else NULL,
			Colv=if (hm != 1 && hm != 4 && hm != 7) {if (hm == 5) order(data[as.numeric(getParam("heatmaprow",1)),,drop=FALSE]) else TRUE} else { order(if (hm==1) as.vector(km$pi) else if (hm==7) as.numeric(xdf$Class)+(km$pi/max(km$pi*2)) else time) },
			col=redgreenblue(hmcol, len=ncols),
			margins = xmargins 
			)
		if (!is.null(forced)) {
			title(sub=forced)
		}
		## dev.off()	
	}
	maxRisk <- toupper(getParam("max_risk", "YES")) %in% c("1","YES","TRUE")
	#if (maxRisk & groups == 2) {
		new_plot(filename=outputfile("MR"))
		o <- order(km$pi[tr])
		plot(1:length(tr), km$p.risk[tr[o]]+1e-20, log="y", xlab="Patients (Ordered by Prognostic Index)", ylab="p-Value",
			main=paste("Risk Group Optimization (+=censored, o=no censored)",ifelse(train.is.dif,"(train only)","")), 
			pch=ifelse(status[tr[o]]==1,1,3), 
			col=xcolor[1:groups][km$pi.cluster[tr[o]]], type="p")
		lines(1:length(tr), km$p.risk[tr[o]]+1e-20, col=8)
		#dc <- diff(km$pi.cluster.original[o])
		#cc <- which(dc != 0)
		#abline(v=cc, col=16)
		#abline(h=km$p.risk[cc[1]], col=16)
		#axis(cc, side=3, col=16, col.axis=16)
		#axis(km$p.risk[cc[1]], side=4, col=16, col.axis=16)
		#text(par("usr")[2], km$p.risk[cc[1]], round(km$p.risk[cc[1]],6), col=16, cex=1, adj=c(1,-0.25))
		#dc <- diff(km$pi.cluster[o])
		#cc <- which(dc != 0)
		#abline(v=cc, col=7)
		#abline(h=km$p.risk[cc[1]], col=7)
		#axis(cc, side=3, col=7, col.axis=7, line=-0.5)
		#axis(km$p.risk[cc[1]], side=4, col=7, col.axis=7, line=-0.5)
		#text(par("usr")[2], km$p.risk[cc[1]], round(km$p.risk[cc[1]],6), col=7, cex=1, adj=c(1,-0.25))

		#legend("topleft",legend=c("By Max Risk","By Size"), lty=c(1,1), col=c(7,16))

		## by prognostic index
		new_plot(filename=outputfile("MP"))
		o <- order(km$pi[tr])
		plot(1:length(o), km$pi[tr[o]], xlab="Patients (Ordered by Prognostic Index)", ylab="Prognostic Index (Risk Score)",
			main=paste("Risk Group Optimization (+=censored, o=no censored)",ifelse(train.is.dif,"(train only)","")), 
			pch=ifelse(status[tr[o]]==1,1,3), 
			col=xcolor[1:groups][km$pi.cluster[tr[o]]], type="h")
		points(1:length(o), km$pi[tr[o]], 
			pch=ifelse(status[tr[o]]==1,1,3), 
			col=xcolor[1:groups][km$pi.cluster[tr[o]]])
		#dc <- diff(km$pi.cluster.original[o])
		#cc <- which(dc != 0)
		#pim <- km$pi[o[cc]]
		#abline(v=cc+0.5, col=16)
		#abline(h=pim, col=16)
		#axis(cc+0.5, side=3, col=16, col.axis=16)
		#axis(pim, side=4, col=16, col.axis=16)
		#text(par("usr")[1], pim, round(pim,6), col=16, cex=1, adj=c(0,-0.25))

		#dc <- diff(km$pi.cluster[o])
		#cc <- which(dc != 0)
		#pim <- km$pi[o[cc]]
		#abline(v=cc+0.5, col=7)
		#abline(h=pim, col=7)
		#axis(cc+0.5, side=3, col=7, col.axis=7, line=-0.5)
		#axis(pim, side=4, col=7, col.axis=7, line=-0.5)
		#text(par("usr")[1], pim, round(pim,6), col=7, cex=1, adj=c(0,-0.25))

		#legend("topleft",legend=c("By Max Risk","By Size"), lty=c(1,1), col=c(7,16))
	#}
	attributes = toupper(getParam("attribute_plot", "YES")) %in% c("1","YES","TRUE")
	if (attributes) {
		new_plot(filename=outputfile("AP"))
		xdf <- data.frame(
				Risk=risks[km$pi.cluster],
				Censored=ifelse(status==0,"Yes","No"),
				Time=round(time,2),
				Prog.Idx=round(as.vector(km$pi),3)
				)
		#xdf <- cbind(xdf, classes[test,])
		if (strataPar != "<SELECT>") {
			xdf <- cbind(xdf, strata=strata)
			colnames(xdf)[ncol(xdf)] <- strataPar
		}
		if (train.is.dif) {
			xdf <- cbind(xdf, Set=ifelse(1:ncol(data) %in% tr,"Train","Test"))
		}
		if (!is.null(headers)) {
			xh <- headers[,,drop=FALSE]
			colnames(xh) <- gsub("^[^.]+\\.","",colnames(xh))
			xdf <- cbind(xdf,xh)
		}
		mx <- matrix(0,ncol=ncol(data),nrow=2)
		colnames(mx) <- colnames(data)
		#plot(0,0,main=paste(xmargins,collapse=","))
		if (1) {
		plot.heatmap.info(mx[,,drop=FALSE], 
			scale="none", 
			ColSideFactors=xdf,
			ColSideFactors.type=c("i","i","h","l","i",rep("?",ncol(xdf)-5)),
			ColSideFactors.col=list(Risk=xcolor[groups:1], Censored=c(8,1)),
			ColvFactor=NULL,
			dendrogram="none",
			Colv=if (hm <= 0) order(as.vector(km$pi)) else hmi$colInd, 
			Rowv=FALSE,
			main="Attributes Plot",
			key=FALSE,
			lhei=c(1,8),
			lwid=c(4,8),
			col=c("white","white"),
			ColSideFactors.height=rep(20/ncol(xdf),ncol(xdf)),
			labRow=rep("",nrow(mx)),
			margins=xmargins
			)
		}
		#title(sub="Attributes Plot")
	}
	if (advanced) {
		library(survivalROC)
		#tco <- 1:9/10 # c(.1,.25,.5,.75,.9) # 2:8/10
		#dt <- diff(range(time))
		#qcutoff <- min(time)+dt*tco # 
		qcutoff <- km$ticks[-1]
		new_plot(filename=outputfile("ROC"))
		plot(0,0,xlim=c(0,1),ylim=c(0,1),type="n",xlab="False Positives",ylab="True Positives",main="SurvivalROC using Prognostic Index (method=KM)")
		rocok <- if (ncol(data) > 500) seq(1,ncol(data),length.out=500) else 1:ncol(data)
		sroc <- list()
		for (i in 1:length(qcutoff)) {
			sroc[[i]] <- list(FP=0, TP=0, AUC=0)
			try(sroc[[i]] <- survivalROC(Stime=time[rocok],  
			    status=status[rocok],      
			    marker=km$pi[rocok],     
			    predict.time =  qcutoff[i], method="KM"))
			lines(sroc[[i]]$FP, sroc[[i]]$TP, col=i)
		}
		abline(0,1)
		legend("bottomright", legend=paste(
			#tco*100,"%, ",
			"time=",round(qcutoff,1),
			", AUC=",round(unlist(lapply(sroc, function(x) x$AUC)),3),
			sep=""), lty=1, col=1:length(qcutoff))

		new_plot(filename=outputfile("ROC2"))
		plot(0,0,xlim=c(0,1),ylim=c(0,1),type="n",xlab="False Positives",ylab="True Positives",main="SurvivalROC using Prognostic Index (method=NNE)")
		sroc <- list()
		for (i in 1:length(qcutoff)) {
			sroc[[i]] <- list(FP=0, TP=0, AUC=0)
			try(sroc[[i]] <- survivalROC(Stime=time[rocok],  
			    status=status[rocok],      
			    marker=km$pi[rocok],     
			    predict.time =  qcutoff[i], method="NNE", lambda=0.01)) # span = 0.25*length(time)^(-0.20)
			lines(sroc[[i]]$FP, sroc[[i]]$TP, col=i)
		}
		abline(0,1)
		legend("bottomright", legend=paste(
			#tco*100,"%, ",
			"t=",round(qcutoff,1),
			", AUC=",round(unlist(lapply(sroc, function(x) x$AUC)),3),
			sep=""), lty=1, col=1:length(qcutoff))
	}
	new_plot(filename=outputfile("BXP"))
	l <- list()	
	p <- numeric(nrow(data))
	for (i in 1:nrow(data)) {		
		for (k in 1:groups) {
			l[[length(l)+1]] <- as.vector(data[i,km$pi.cluster==k])
		}
		xdf <- data.frame(x=data[i,], cluster=km$pi.cluster)
		sm <- summary(lm(x~cluster, data=xdf))$fstatistic
		p[i] = pf(sm[1], sm[2], sm[3], lower.tail=FALSE)
	}
	#names(l) <- paste(rep(symb, each=groups),1:groups,sep=":")
	xn <- rep("",length(l))
	pos <- (1:nrow(data)-1)*(groups)+1
	xn[pos] <- symb
	xn[pos+1] <- paste("p=",format(p,digits=3),sep="")
	names(l) <- xn
	par(mar=c(8,5,8,2))
	mxl <- max(unlist(l), na.rm=TRUE)
	mnl <- min(unlist(l), na.rm=TRUE)
	boxplot(l, border=rep(xcolor[1:groups], nrow(data)), las=3, 
		main="", xaxt="n", ylab="Gene Expression Level", 
		ylim=c(mnl, mxl + 0.1 * max(0.1,range(c(mnl, mxl)))))
	title(main="Gene Expression By Risk Group", line=7)
	#axis(side=1,at=(1:nrow(data)-1)*(groups)+1+(groups-1)/2,labels=paste(symb,paste("p=",format(p,digits=3),sep=""),sep="\n"),las=3, tick=FALSE)
	axis(side=1,at=(1:nrow(data)-1)*(groups)+1+(groups-1)/2,labels=symb,las=3, tick=FALSE)
	axis(side=3,at=(1:nrow(data)-1)*(groups)+1+(groups-1)/2,labels=paste("p=",format(p,digits=3),sep=""),las=3, tick=FALSE)
	abline(v = 0.5+(0:nrow(data))*groups, col=8)
	legend("top", legend=risks, pch=22, col=xcolor[1:groups], text.col=xcolor[1:groups], ncol=min(5, groups))

	mDEG <- data.frame(Symbol=symb, p.DEG=p)


	if (strataPar != "<SELECT>") {
		sink(file=fittingfile, append=TRUE, type="output")
		xt <- table(strata)
		print("*** PERFORMING STRATIFICATION ***<BR>\n")
		print(data.frame(Stratas=names(xt), Occurrences=xt))
		sink()

		xu <- sort(unique(strata))
		new_plot(filename=outputfile(paste("ORGS",i,sep="")))
		try(
			plot.raw.kaplan.per.strata(time, status, strata, col=palette()[9:length(palette())],
			main=paste("Original Groups by ",strataPar," (No covariate fitting)",sep="")) 
			)
		if (separated) {
		  	for (i in 1:length(xu)) {
				w <- which(strata == xu[i])
				#w <- which(ok %in% w)
				if (length(w) >= 4*groups && nchar(xu[i]) > 0) {
					new_plot(filename=outputfile(paste("SZ",i,sep="")))
					try(
						skm <- plot.kaplan(data, time=time, status=status, 
						model=1:nrow(data), risk.groups=groups, 
						col=xcolor, col2=survcolor,
						main=paste("By",strataPar,"=",xu[i],paste("(",length(w),")",sep=""),"\n(Overall betas",ifelse(maxRisk,", and re-maximizing groups",""),")"), 
						symbols=symb, tr=w, tr.betas=tr, #betas=km$betas,
						draw.text=fitting.info,
						draw.main=fitting.info) 
						)
				}
			}			
		}
		if (!separated) {
			# this code needs to be cleaned since it also considers separated=TRUE
			oki <- c()
			okp <- c()
			okci <- c()
			okcol <- c()
			st <- table(strata[ok])
			stn <- sum(st >= 4*groups)
		  	for (i in 1:length(xu)) {
				w <- which(strata == xu[i])
				#w <- which(ok %in% w)
				if (length(w) >= 4*groups && nchar(xu[i]) > 0) {
					oki <- c(oki, i)
					p.i <- km$betas %*% data[, w, drop=FALSE]
					if (maxRisk) {						
						mrg <- maximizeRiskGroups(time[w], status[w], p.i, groups)
						pi.cluster <- if (maxRisk) mrg$cluster.max else mrg$cluster.quantile ## old code, only if will be used, new if (maxRisk) made to improve runtime
					} else {
						qcuts <- quantile(p.i, seq(0,1,length.out=groups+1, na.rm=TRUE))
						pi.cluster <- cut(p.i, breaks=qcuts, labels=FALSE, include.lowest=TRUE)
					}
					newplot <- (i==1  ||  separated)
					if (newplot) {
						new_plot(filename=outputfile(paste("S",i,sep="")))
						if (!separated) {
							pmar <- par("mar")
							pmar[1] <- pmar[1] + max(0,stn-3)
							par(mar=pmar)
						}
					}
					plotf <- if (newplot) plot else lines	
					surdif <- survdiff(Surv(time, status) ~ pi.cluster, data=data.frame(time=time[w], status=status[w], cluster=pi.cluster))
					p <- 1-pchisq(surdif$chisq,df=groups-1)
					okp <- c(okp, p)
					ci <- concordance.index(concordance.cox(km$betas, time[w], status[w], data[, w, drop=FALSE]))
					okci <- c(okci, ci)
					icol <- 8 + length(oki)
					okcol <- c(okcol, icol)
					plotf(survfit(Surv(time[w], status[w]) ~ pi.cluster), 
						col=icol, mark.time=TRUE, pch=mark.char,
						lty=1:groups, ylim=if (separated) c(0, 1.1) else c(-0.15, 1.15), conf.int=FALSE, lwd=ifelse(separated, 2, 3))
					xt <- axTicks(1)
					xt <- sort(c(xt, xt[-length(xt)]+diff(xt)/2))
					if (separated) {
						for (k in 1:groups) {
							axis(1, xt, labels=sapply(xt, function(x,times) sum(times >= x), times=time[w][pi.cluster == k]),
						 tick=FALSE, line=k, col.axis=icol)
						}
					} else {
						axis(1, xt, labels=sapply(xt, function(x) sum(time[w] >= x)),
						 tick=FALSE, line=i, col.axis=icol)
					}
					if (separated) {
						legend("top", legend=risks, col=icol, lty=1:groups, ncol=min(3,groups), cex=0.75)
						title(paste("By",strataPar,"=",xu[i],paste("(",length(w),")",sep=""),"\n(Using overall betas",ifelse(maxRisk,", and re-maximizing groups",""),")\nci=",round(ci*100,2),"|p=",p)) #(equal sized risk groups)\n
						## dev.off()
					}
				}
			}
			if (!separated) {
				legend("bottom", legend=paste(xu[oki],",p=",format(okp,digits=2),sep=""), col=okcol, lty=1, ncol=min(3,length(oki)/2), cex=0.75)
				legend("top", legend=paste(xu[oki],",ci=",format(okci,digits=2),sep=""), col=okcol, lty=1, ncol=min(3,length(oki)/2), cex=0.75)
				legend("topright", legend=risks, col=rep(1,groups), lty=1:groups, cex=0.75)
				title(paste("By ",strataPar, "\n(Using overall betas",ifelse(maxRisk,", and re-maximizing groups)",")"),sep=""))
				## dev.off()
			}
		}
		#if (separated) {
		  	for (i in 1:length(xu)) {
				w <- which(strata == xu[i])
				#w <- which(ok %in% w)
				if (length(w) >= 4*groups && nchar(xu[i]) > 0) {
					new_plot(filename=outputfile(paste("Z",i,sep="")))
					skm <- list(cox=NA)
					try(
					skm <- plot.kaplan(data, time=time, status=status, 
						model=1:nrow(data), risk.groups=groups, 
						col=xcolor, col2=survcolor,
						main=paste("By ",strataPar," = ",xu[i],paste(" (",length(w),")",sep=""),"\n(Betas re-estimated in this set",ifelse(maxRisk,", and re-maximizing groups",""),")",sep=""), 
						symbols=symb, tr=w, tr.betas=w,
						draw.text=fitting.info,
						draw.main=fitting.info)
						)
						sink(file=fittingfile, append=TRUE, type="output")
						if (!isPDF) cat("<TABLE><TR><TD>")
						cat("**** RE-ESTIMATION IN STRATA ",xu[i],"****<BR>\n")
						print(summary(skm$cox))
						if (!isPDF) cat("</TD></TR></TABLE>")
						sink()
				}
			}
		#}
	}

	if (getParam("network","none") != "none"  && getParam("network","none") != "") {
		netFile <- getParam("network","none")
		new_plot(filename=outputfile("Net1"))
		if (file.exists(netFile)) {
			netcnx <- read.delim(netFile, header=FALSE, as.is=TRUE)

			library(igraph)

			get.net <- function(cnx, ids, cytoscape=FALSE) {
				xids <- sort(unique(unlist(ids)))
				rows <- which(cnx[,1] %in% xids & cnx[,2] %in% xids)
				if (length(rows) == 0) return (NULL)
				if (cytoscape) {
					library(graph)
					g <- new ('graphNEL', edgemode='undirected')
					for (i in xids) g <- graph::addNode(i, g)
					for (i in rows) g <- graph::addEdge(cnx[i,1], cnx[i,2], g)
					g
				} else {
					n <- matrix(0, nrow=length(xids), ncol=length(xids))
					colnames(n) <- rownames(n) <- xids
					n[as.matrix(cnx[rows,1:2])] <- 1
					n[as.matrix(cnx[rows,2:1])] <- 1
					n
				}
			}

			net <- get.net(netcnx, symb)
			if (!is.null(net)) {
				xcol <- redgreenblue(c("green3","white","red"), len=nrow(nrow))#redgreenblue(c("green3","yellow","red"), len=nrow(n2))
				
				plot.igraph(graph.adjacency(net, mode="undirected"), 
					vertex.color="salmon",#xcol[rank(nodes.weight[rownames(n2)])], 
					vertex.size=max(min(20,20-5*nrow(net)/20),8), 
					main=paste("Network (", gsub("\\..+","",gsub(".+/","",netFile)), ")",sep=""),
					vertex.label.color="black", 
					vertex.label.dist=max(min(20,20-5*nrow(net)/20),8)*1.05/20,
					edge.width=1 #,edge.label=ltn3[ltn3 > 0]
					)
			} else {
				plot(0,0,main="No Connections Found")
			}
		} else {
			plot(0,0,main="No Net File")
		}
	}

	##################### BEGIN: WRITE TABLES 
	xdf <- data.frame(Sample=1:ncol(data),
					Name=colnames(data),
					Time=time,
					Status=status,
					Prognostic.Index=round(allPI,4),
					Risk.Rank=rank(-allPI), 
					Risk.Group=risks[km$pi.cluster])
	if (isPDF) {
		plot.text(km$fitting[,1:ifelse(advanced || fitting.info,ncol(km$fitting),3),drop=FALSE])
		plot.text(xdf[order(xdf$Risk.Rank)[1:min(100,nrow(xdf))],,drop=FALSE])
		xtdf <- data.frame(R.Model.Fitting.Output=readLines(fittingfile))
		xpages <- ceiling(nrow(xtdf)/100)
		for (p in 1:xpages) {
			fr <- (p-1)*100+1
			to <- min(p*100,nrow(xtdf))
			plot.text(xtdf[fr:to,,drop=FALSE], family="mono", max.chars=120)
		}
		
	} else {
		fitfile <- outputfile("F",ext=".html")
		tow <- as.HTML(km$fitting[,1:ifelse(advanced || fitting.info,ncol(km$fitting),3),drop=FALSE], row.names=FALSE)
		if (advanced || fitting.info) {
			tow <- c(as.HTML(data.frame("Fitting Information or Advanced Plots"=c("Please see the Tutorial for details.","Thanks!")), row.names=FALSE),tow)
		}
		tow <- c(tow, as.HTML(xdf[order(xdf$Risk.Rank),,drop=FALSE], row.names=FALSE))
		tow <- c(tow, as.HTML(data.frame(R.Model.Fitting.Output=gsub(" ","&nbsp;",readLines(fittingfile))), row.names=FALSE, prefix="<CODE>", posfix="</CODE>"))
		writeLines(tow, fitfile)
		plots_made <- c(plots_made[1], fitfile, plots_made[-1])
	}
	unlink(fittingfile)
	##################### END: WRITE TABLES 

	########################### Begin: Save Summary Table:
	resumen <- t(
			data.frame(Biomarker=paste(symb, collapse=", "),
					Database=getParam("database",""),
					Genes=length(symb),
					CI=km$ci,
					Log.Rank=km$log.rank,
					Hazard.Ratio=km$hr,
					C.I.Hazard.Ratio=paste("[", km$hrl, " - ", km$hru, "]", sep=""),
					pHR=km$hr.p,
					Correlation=km$rsq,
					Significant.Genes=sum(summary(km$cox)$coefficients[,5] <= 0.05),
					Marginal.Genes=sum(summary(km$cox)$coefficients[,5] <= 0.1),
					Cox.Interesting=paste(symb[summary(km$cox)$coefficients[,5] <= 0.1], collapse=", "),
					DEG=sum(p <= 0.05),
					DEG.Interesting=paste(symb[p <= 0.05], collapse=", "),
					Notes="Hazard Ratio was estimated by fitting a CoxPH using risk group as covariate.")
			)
	colnames(resumen) <- c("Summary.Result")
	if (isPDF) {
		plot.text(data.frame(Info=rownames(resumen), Summary.Result=resumen[,1]), max.chars=120)
		plot.text(mDEG, max.chars=120)
	} else {
		ff <- c(readLines(fitfile), as.HTML(resumen), as.HTML(cbind(mDEG,"."="")))
		writeLines(ff, fitfile)
	}
	########################### End: Save Summary Table:


}

if (xfunc == "correlation") {
	norm <- as.numeric(getParam("normalize", 1))
	if (norm != 0) {
		data <- apply(data, norm, scale)
		## data <- t(apply(t(data), 2, scale))
	}
	predictor <- getParam("predictor", NULL)
	if (is.null(predictor) || is.na(predictor)  || 
		is.null(getParam(predictor))  ||  is.na(getParam(predictor))) {
		kolor <- 1
	} else {
		klasses <- unlist(strsplit(getParam(predictor),sep))[columns=="D"]
		if (any(nchar(klasses) < 1)) {
			data <- data[, nchar(klasses) > 0, drop=FALSE]
			klasses <- klasses[nchar(klasses) > 0]
		}
		if (length(unique(klasses)) > length(palette()))
			kolor <- rainbow(length(unique(klasses))+2)[as.numeric(factor(klasses))]
		else
			kolor <- as.numeric(factor(klasses))
	}
	new_plot(filename=outfile)
	panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
	{
	    #usr <- par("usr"); on.exit(par(usr))	    
	    #par(usr = c(0, 1, 0, 1))
	    r <- abs(cor(x, y))
	    txt <- format(c(r, 0.123456789), digits=digits)[1]
	    txt <- paste(prefix, txt, sep="")
	    if(missing(cex.cor)) cex.cor <- 1
	    points(x,y,...)
	    text(par("usr")[1],par("usr")[4],txt,cex=cex.cor,adj=c(0,1))
	    #text(0.5, 0.5, txt, cex = cex.cor)
	}	
	panel.hist <- function(x, col="cyan", main="", ...)
	{
	    usr <- par("usr"); on.exit(par(usr))
	    par(usr = c(usr[1:2], 0, 1.5) )
	    h <- hist(x, plot = FALSE)
	    breaks <- h$breaks; nB <- length(breaks)
	    y <- h$counts; y <- y/max(y)
	    rect(breaks[-nB], 0, breaks[-1], y, col="cyan",  ...)
	}
	## from default panel.smooth
	panel.smooth1 <- function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
	    cex = 1, col.smooth = "red", span = 2/3, iter = 3, 
	    digits=2, prefix="", cex.cor, ...) 
	{
	    points(x, y, pch = ".", col = col, bg = bg, cex = cex * 2)
	    r <- abs(cor(x, y))
	    txt <- format(c(r, 0.123456789), digits=digits)[1]
	    txt <- paste(prefix, txt, sep="")
	    if(missing(cex.cor)) cex.cor <- 1.5
	    text(par("usr")[1],par("usr")[3],txt,cex=cex.cor,adj=c(-0.1,-0.2))
	    ok <- is.finite(x) & is.finite(y)
	    if (any(ok)) {
	        lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
	            col = col.smooth, lwd=1.5, ...)
	        if (length(unique(col)) > 1) {
	        	for (i in unique(col)) {
	        		ok <- is.finite(x) & is.finite(y) & col == i
			        lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
			            col = i, ...)
	        	}
	        }
       } 
	}
	panel.smooth2 <- function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
	    cex = 1, col.smooth = "red", span = 2/3, iter = 3, 
	    digits=2, prefix="", cex.cor, ...) 
	{
	    points(x, y, pch = pch, col = col, bg = bg, cex = cex)
	    r <- abs(cor(x, y))
	    txt <- format(c(r, 0.123456789), digits=digits)[1]
	    txt <- paste(prefix, txt, sep="")
	    if(missing(cex.cor)) cex.cor <- 1.25
	    text(par("usr")[1],par("usr")[4],txt,cex=cex.cor,adj=c(-0.1,1.1))
	    ok <- is.finite(x) & is.finite(y)
	    if (any(ok)) {
	        lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
	            col = col.smooth, lwd=1.5, ...)
	        if (length(unique(col)) > 1) {
	        	for (i in unique(col)) {
	        		ok <- is.finite(x) & is.finite(y) & col == i
			        lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
			            col = i, ...)
	        	}
	        }
       } 
	}
	pairs(t(data), gap=0, col=kolor, pch=20, 
		cex=min(max(1000/ncol(data),.25),2),
		main=paste(predictor,paste(symb,collapse="|"),sep="\n"), 
		lower.panel=panel.smooth2, col.smooth="burlywood3", prefix="r=",
		upper.panel=panel.smooth1, 
		diag.panel=panel.hist, cex.labels=1.25
		)
	if (length(unique(predictor)) <= length(palette())) {
	}
	new_plot(filename=outputfile("HC"))
	cl = round(sqrt(nrow(data)))
	hc <- plot.hclust.colored.clusters(as.dist(1-cor(t(data))), 
		k=cl, 
		silhouettes=TRUE,
		label.attr=list(lab.cex=1),
		sil.image=FALSE,
		main=paste("k=SQRT(",nrow(data),") ~ ",cl,sep=""))

	panel.text <- function (x, y, col = par("col"), bg = NA, cex.text=1, text.labels=rep("",length(x)), ...) 
	{
		points(x, y, pch=4, cex=1, col=col, bg=bg)
	    text(x, y, text.labels, col = col, bg = bg, cex = cex.text, adj=c(0.5, 0.5))
	}
	new_plot(filename=outputfile("MDS"))
	xcmd <- cmdscale(as.dist(1-cor(t(data))), k=cl)
	rownames(xcmd) <- colnames(xcmd) <- NULL
	colnames(xcmd) <- paste("MDS",1:ncol(xcmd))
	pairs(xcmd, 
		col=hc$cluster,
		gap=0, 
		upper.panel=panel.text,
		text.labels=gsub(":.+","",rownames(data)),
		main="MDS")

	tabF <- outputfile("ct",ext=".html")
	writeLines(as.HTML(cor(t(data))), tabF)
	plots_made <- c(plots_made,tabF)
}

if (xfunc == "classification") {
	norm <- as.numeric(getParam("normalize", 1))
	predictor <- getParam("predictor", NULL)
	strata <- getParam("strata", NULL)
	classname <- predictor
	if (is.null(predictor) || is.na(predictor)  || 
		is.null(getParam(predictor))  ||  is.na(getParam(predictor))) {
	} else {
		predictor <- unlist(strsplit(getParam(predictor),sep))[columns=="D"]
	}
	if (!is.null(strata)) {
		strata <- unlist(strsplit(getParam(strata),sep))[columns=="D"]		
	}
	if (norm) {
		data <- t(apply(t(data), 2, scale))
	}
	wd <- which(nchar(trim(predictor)) > 0)
	predictor <- predictor[wd]
	if (!is.null(strata)) {
		strata <- strata[wd]
		strats <- table(strata)
		strataOK <- names(strats)[strats > 2]
	}
	data <- data[,wd,drop=FALSE]
	library(galgo)
	if (file.exists("galgoS3.r")) source("galgoS3.r")
	train <- eval(parse(text=paste("c(",getParam("trainset", paste("1:",ncol(data),sep="")),")")))
	test <- eval(parse(text=paste("c(",getParam("testset", paste("1:",ncol(data),sep="")),")")))
	# c("knn","mlhd","svm","nearcent","rpart","nnet","ranforest"),

	doConfigBB <- function(predictor, data, info) {
		configBB.VarSel(data=data,
			classes=as.factor(predictor),
			classification.method=getParam("classifier","knn"),
			classification.rutines=ifelse(getParam("classifier","knn") %in% c("svm","rpart","nnet","ranforest"),"R","C"),
			populationSize=nrow(data),
			main=paste(info, "Class=",classname,sep=""),
			train=rep(1/3,20),
			scale=FALSE
			)		
	}

	bb <- doConfigBB(predictor, data, "") 
	if (any(nchar(c(getParam("trainset",""),getParam("testset",""))) > 1)) {
		bb$data$splitTrain <- list(train)
		bb$data$splitTest <- list(test)
	}
	#confusion, confusionpamr, confusionbox
	if (getParam("confusion","no") != "no") {
		new_plot(filename=outputfile("C"))
		plot(bb, type=getParam("confusion","no"), chromosomes=c(1:nrow(data)), set=eval(parse(text=paste("c(",getParam("confusion_set", paste("1:",ncol(data),sep="")),")"))))
		if (!is.null(strata)) {
			for (i in 1:length(strataOK)) {
				s <- strataOK[i]
				bb2 <- doConfigBB(predictor[strata==s], data[,strata==s,drop=FALSE], paste("Strata=",s,", ")) 
				new_plot(filename=outputfile(paste("C",i,sep="")))
				plot(bb2, type=getParam("confusion","no"), chromosomes=c(1:nrow(data)), set=eval(parse(text=paste("c(",getParam("confusion_set", paste("1:",ncol(data),sep="")),")"))))				
			}
		}
	}
	if (getParam("heatmap","1") == "1") {
		new_plot(filename=outputfile("H"))
		heatmapModels(bb, models=1:nrow(data), margins=c(5, 15))
		new_plot(filename=outputfile("H2"))
		plot.heatmap.info(data, ColSideFactors=data.frame(Class=predictor), ColvFactor=1, margins=c(5, 20))
	}
	if (getParam("PCA","0") == "1") {
		new_plot(filename=outputfile("B"))
		pcaModels(bb, models=1:nrow(data), show.loadings=FALSE)
		if (getParam("PCA","0") == "1") {
			new_plot(filename=outputfile("B2"))
			pcaModels(bb, models=1:nrow(data), show.loadings=TRUE,)
		}
	}
	#geneprofiles, genevalues, genevalueslines, sampleprofiles, genevaluesbox
	theouts <- c("D","E","F","G","I")
	thetype <- c("geneprofiles","gene_values", "genevalueslines", "sampleprofiles", "genevaluesbox")
	for (i in 1:length(thetype)) {
		if (getParam(thetype[i],"0") == "1") {
			mar <- par("mar") ## to avoid galgo bug
			par(mar=mar+c(2,0,0,0))
			new_plot(filename=outputfile(theouts[i]))
			plot(bb, y=1:nrow(data), type=gsub("_","",thetype[i]), chromosomes=list(1:nrow(data)))
		}
	}
}

if (xfunc == "gene-survival-permutations") {
	source("gene-survival-permutations.r")
}

if (xfunc == "cna-survival-groups") {
	source("CNAriskGroups.R")
}

if (toupper(getParam("zipscript", "NO")) %in% c("1","YES","TRUE","ON","CHECKED")) {
	zipFile <- outputfile("zip", ext=".zip")
	withinZip <- c(inputfile, "loadR2LS.r", "R2LS.r", "runweb.r")
	if (exists("zip") && is.function(zip)) {
		zip(zipFile, withinZip)
	} else {
		system(paste("zip",zipFile,paste(withinZip,collapse=" "),sep=" "))
	}
	plots_made <- c(plots_made, zipFile)
}

if (length(.errores) > 0) {
	new_plot(filename=outputfile("ERR"))
	names(.errores) <- NULL
	edf <- data.frame(.errores, stringsAsFactors=FALSE, check.names=FALSE)
	colnames(edf) <- rep("Error List",ncol(edf))
	plot.text(edf,max.chars=80)
}

end_plots()

if (file.exists(logFile)) {
	inputParams <- c(inputFile=inputfile, endTime=date())
	xf <- file(logFile, open="at")
	writeLines(substr(paste(names(inputParams),"=",inputParams),1,80), xf)
	close(xf)
}
cat(paste(plots_made, collapse=";"))
cat("\n")


# http://bioinformatica.mty.itesm.mx/sites/default/files/HPRD_PPI.txt
# http://bioinformatica.mty.itesm.mx/sites/default/files/PPI_HPRD_BIOGRID.txt
# /var/www/html/drupal/sites/default/files/HPRD_PPI.txt
# /var/www/html/drupal/sites/default/files/PPI_HPRD_BIOGRID.txt
