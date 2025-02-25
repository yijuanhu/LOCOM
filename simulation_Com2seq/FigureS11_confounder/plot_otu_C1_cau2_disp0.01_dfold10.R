n_confounder = 1
causal_type = 2

plot_global_otu <- function(res_sub, ylab1="", ylab2="", ylab3="") {
    
    sub1 <- c(2:7)
    beta_range1 = range(res_sub$beta[sub1])
    d1=(beta_range1[2]-beta_range1[1])*0.03
    plot(NA, main="", xlab="", ylab=ylab1, ylim=c(0.05, 1), xlim=c(beta_range1[1]-d1, beta_range1[2]+d1), bty="n")
    box(lwd=1.5)
    points(res_sub$beta[sub1], res_sub$p_comp[sub1], col="green3", type="o", lwd=1.5, lty=1, pch=17)
    #points(res_sub$beta[sub1], res_sub$p_compHM[sub1], col="green3", type="o", lwd=1.5, lty=2, pch=25)
    points(res_sub$beta[sub1], res_sub$p_pool[sub1], col="purple", type="o", lwd=1.5, lty=1, pch=18)
    points(res_sub$beta[sub1], res_sub$p_16S[sub1], col="blue", type="o", lwd=1.5, lty=1, pch=15)
    points(res_sub$beta[sub1], res_sub$p_shot[sub1], col="lightblue", type="o", lwd=1.5, lty=1, pch=19)
    points(res_sub$beta[sub1], res_sub$p_new[sub1], col="red", type="o", lwd=1.5, lty=2, pch=16)
    points(res_sub$beta[sub1], res_sub$p_new_wt1[sub1], col="red", type="o", lwd=1.5, lty=1, pch=16)
    points(res_sub$beta[sub1], res_sub$p_new_omni[sub1], col="black", type="o", lwd=1.5, lty=1, pch=4)
    
    sub2 <- c(3:8)
    beta_range1 = range(res$beta[sub2])
    d1=(beta_range1[2]-beta_range1[1])*0.03
    plot(NA, main="", xlab="", ylab=ylab2, ylim=c(0.1, 0.9), xlim=c(beta_range1[1]-d1, beta_range1[2]+d1), bty="n")
    box(lwd=1.5)
    points(res_sub$beta[sub2], res_sub$pcom_sen[sub2], col="green3", type="o", lwd=1.5, lty=1, pch=17)
    points(res_sub$beta[sub2], res_sub$pcomHM_sen[sub2], col="green3", type="o", lwd=1.5, lty=2, pch=25)
    points(res_sub$beta[sub2], res_sub$pool_sen[sub2], col="purple", type="o", lwd=1.5, lty=1, pch=18)
    points(res_sub$beta[sub2], res_sub$S16_sen[sub2], col="blue", type="o", lwd=1.5, lty=1, pch=15)
    points(res_sub$beta[sub2], res_sub$shot_sen[sub2], col="lightblue", type="o", lwd=1.5, lty=1, pch=19)
    points(res_sub$beta[sub2], res_sub$new_sen[sub2], col="red", type="o", lwd=1.5, lty=2, pch=16)
    points(res_sub$beta[sub2], res_sub$newwt1_sen[sub2], col="red", type="o", lwd=1.5, lty=1, pch=16)
    points(res_sub$beta[sub2], res_sub$newOmni_sen[sub2], col="black", type="o", lwd=1.5, lty=1, pch=4)
    
    plot(NA, main="", xlab="", ylab=ylab3, ylim=c(0, 1), xlim=c(beta_range1[1]-d1, beta_range1[2]+d1), bty="n")
    box(lwd=1.5)
    abline(h=0.2, col="gray", lwd=1.5, lty=2)
    points(res_sub$beta[sub2], res_sub$pcom_fdr[sub2], col="green3", type="o", lwd=1.5, lty=1, pch=17)
    points(res_sub$beta[sub2], res_sub$pcomHM_fdr[sub2], col="green3", type="o", lwd=1.5, lty=2, pch=25)
    points(res_sub$beta[sub2], res_sub$pool_fdr[sub2], col="purple", type="o", lwd=1.5, lty=1, pch=18)
    points(res_sub$beta[sub2], res_sub$S16_fdr[sub2], col="blue", type="o", lwd=1.5, lty=1, pch=15)
    points(res_sub$beta[sub2], res_sub$shot_fdr[sub2], col="lightblue", type="o", lwd=1.5, lty=1, pch=19)
    points(res_sub$beta[sub2], res_sub$new_fdr[sub2], col="red", type="o", lwd=1.5, lty=2, pch=16)
    points(res_sub$beta[sub2], res_sub$newwt1_fdr[sub2], col="red", type="o", lwd=1.5, lty=1, pch=16)
    points(res_sub$beta[sub2], res_sub$newOmni_fdr[sub2], col="black", type="o", lwd=1.5, lty=1, pch=4)
    if (causal_type==1) legend(x=beta_range1[1]-d1, y=1, 
                               c("New-omni", "New-we", "New-wc", "Com-count", "Com-p-C", "Com-p-HM", "16S", "SMS"), 
                               col=c("black", "red", "red", "purple", "green3", "green3", "blue", "lightblue"), 
                               pch=c(4, 16, 16, 18, 17, 25, 15, 19),
                               lty=c(1, 1, 2, 1, 1, 2, 1, 1),
                               lwd=c(1.5), cex=1.1, bty="n")

}


for (disp1 in c(0.01)) {
    for (depth_fold in c(10)) {

        res = read.table(paste("summary_forR_C", n_confounder, "_cau", causal_type, "_disp", disp1, "_dfold", depth_fold, ".txt", sep=""), header=TRUE, as.is=TRUE)
        res = res[order(res$n_conf, res$beta),]
        res_sub = res[which(res$n_conf==n_confounder),]
        pdf_file = paste("plot_n100_C", n_confounder, "_cau", causal_type, "_disp", disp1, "_dfold", depth_fold, ".pdf", sep="")
        
        
        pdf(pdf_file, height=8.5, width=3)
        par(mfcol=c(3,1), pty="s", cex.lab=1.5, mar=c(2.5, 4.5, 1, 0.5))
        if (causal_type==1) {
            plot_global_otu(res_sub=res_sub, ylab1="Power of global test", ylab2="Sensitivity", ylab3="FDR")
        } else {
            plot_global_otu(res_sub=res_sub)
        }
        dev.off()
        
        
    }
}

