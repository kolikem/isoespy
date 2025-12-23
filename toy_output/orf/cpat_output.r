load("/Users/kikemoto/Desktop/PROJECTS/IsoformVisualizer/Revision_1/isoespy1.1.1/test_data/Human_logitModel.RData")
test <- read.table(file="cpat_output.ORF_info.tsv",sep="\t",header=T)
test$Coding_prob <- predict(mylogit,newdata=test,type="response")
write.table(test, file="cpat_output.ORF_prob.tsv", quote=F, sep="\t",row.names=FALSE, col.names=TRUE)
