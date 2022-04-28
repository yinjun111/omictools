#####
#omictools
#####

#Test omictools workflow

#process
omictools rnaseq-process -c /data/jyin/Pipeline_Test/Testdata/rnaseq/Test_samples_anno_v1.txt -t Mouse.B38.Ensembl88 -r cluster

#merge
omictools rnaseq-merge -c /data/jyin/Pipeline_Test/Testdata/rnaseq/Test_samples_anno_v1.txt -t Mouse.B38.Ensembl88 -r cluster

#de
omictools rnaseq-de -c /data/jyin/Pipeline_Test/Testdata/rnaseq/Test_samples_anno_v1.txt -t Mouse.B38.Ensembl88 --treatment KO --reference WT -r local

#summary
omictools rnaseq-summary -c /data/jyin/Pipeline_Test/Testdata/rnaseq/Test_samples_anno_v1.txt -t Mouse.B38.Ensembl88 --run_gsea-gen cluster
