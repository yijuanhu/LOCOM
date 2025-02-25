n_confounder = 0

plot_size <- function(res, main=NULL, names=NULL, ylabel=NULL) {
    barplot(main=main, height=c(res$p_new_omni, res$p_new_wt1, res$p_new), 
            ylim=c(0, 0.1), border=F, 
            col=c("red", "black", "black"),
            names.arg=names, las=3, ylab=ylabel, 
            density=c(-1,-1,40))
    # box(lwd=1.5)
    abline(h=0.05, col="gray", lwd=1.5, lty=2)
}

names = c("Com-2seq (o)", "Com-2seq (e)", "Com-2seq (c)")

for (depth_fold in c(10, 1)) {
    for (disp1 in c(0.01, 0.001)) {
        pdf(paste("barplot_n100_C", n_confounder, "_disp", disp1, "_dfold", depth_fold, "_omni.pdf", sep=""), height=3, width=8)
        par(mfcol=c(1,3), pty="s", cex.lab=1.2, cex.main=1.2, mar=c(6, 4, 1, 2) + 0.1)
        for (causal_type in c(1, 2, 3)) {
            res=read.table(paste("summary_forR_C", n_confounder, "_cau", causal_type, "_disp", disp1, "_dfold", depth_fold, ".txt", sep=""), header=TRUE, as.is=TRUE)
            res=res[which(res$n_conf==n_confounder & res$beta==1.0),]
            main=paste("M", causal_type, sep="")
            ylabel = ifelse(causal_type==1, "Type I error of global test", "")
            
            if (depth_fold==10 & disp1==0.01) plot_size(res=res, main=main, names=names, ylabel=ylabel)
            if (depth_fold==1 & disp1==0.01) plot_size(res=res, main=main, names=names, ylabel=ylabel)
            if (depth_fold==10 & disp1==0.001) plot_size(res=res, main=main, names=names, ylabel=ylabel)
            if (depth_fold==1 & disp1==0.001) plot_size(res=res, main=main, names=names, ylabel=ylabel)
        }
        dev.off()
    }
}

