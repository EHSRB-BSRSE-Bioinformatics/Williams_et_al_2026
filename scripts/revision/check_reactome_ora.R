source("R/revision/reactome_ora.R")
r <- data.frame(log2FoldChange = c(1, 2, -3, 4, -1, 5, NA),
 padj = c(.01, .2, .01, .4, .01, .01, .01),
 row.names = c("ALPHA_p1", "ALPHA_p2", "BETA_p1", "GAMMA_p1", "DELTA_p1", "UNMAPPED_p1", "EPSILON_p1"))
map <- list(ALPHA = "1", BETA = "2", GAMMA = "2", DELTA = c("3", "4", "3"))
z <- historical_ora_input(list(DESeqRes = r), map)
stopifnot(setequal(z$selected, c("2", "3", "4")), setequal(z$background, as.character(1:4)))
# First tied probe wins; strict padj and fold-change boundaries; no final DEG use.
r <- data.frame(log2FoldChange = c(1, -1, log2(1.5), 1, 1), padj = c(.05, .01, .01, .049, NA),
 row.names = c("A_1", "A_2", "B_1", "C_1", "D_1"))
z <- historical_ora_input(list(DESeqRes = r, DESeqResDEGs = "must not read"), list(A="1", B="2", C="3", D="4"))
stopifnot(identical(z$selected, "3"), length(z$background) == 4L)
stopifnot(sum(c(.05-.Machine$double.eps,.05,.05+.Machine$double.eps,NA)<.05,na.rm=TRUE)==1L)
p <- list(small=as.character(1:14), lower=as.character(1:15), upper=as.character(1:500), large=as.character(1:501))
input <- list(selected=as.character(1:10), background=as.character(1:600))
x <- reactome_hypergeom(input, p)
y <- reactome_hypergeom(input, p, compiled=compile_reactome(p))
stopifnot(isTRUE(all.equal(x,y)), identical(x$terms$pathway, c("lower", "upper")),
 identical(x$terms$padj, p.adjust(x$terms$pvalue,"BH")),
 isTRUE(all.equal(x$terms$pvalue[1], phyper(9,15,486,10,lower.tail=FALSE))))
stopifnot(isTRUE(all.equal(x$terms$pvalue[1],exp(lchoose(15,10)-lchoose(501,10)))))
# A zero-overlap term must remain in BH: the hit term is raw-p positive but
# fails adjustment over the two-term family (a hit-only family would reject).
f <- list(hit=as.character(1:15),miss=as.character(16:30))
z <- reactome_hypergeom(list(selected=as.character(1:4),background=as.character(1:30)),f)
stopifnot(nrow(z$terms)==2L,z$terms$pvalue[1]<.05,z$terms$padj[1]>.05,
 z$terms$pvalue[2]==1,z$status=="none")
input$selected <- character()
stopifnot(reactome_hypergeom(input,p)$status == "insufficient_input",
 reactome_hypergeom(input,p)$n_significant == 0L)
ids <- paste0("TK6_3_3_",1:1000)
df <- trial_identity(ids); df$status <- "insufficient_input"
check_trial_coverage(df,ids)
stopifnot(mean(rep(FALSE,nrow(df))) == 0,
 inherits(try(check_trial_coverage(df[-1,],ids),silent=TRUE),"try-error"),
 inherits(try(trial_identity(c(ids,ids[1])),silent=TRUE),"try-error"))
cat("Mapping, boundaries, family sizes, direct/compiled equivalence, identity and denominator tests passed.\n")
