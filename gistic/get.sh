
cp ~/share/exigens/brain-mets/gistic/absolute/*.{sh,R,m} .
cp ~/share/exigens/brain-mets/gistic/absolute/*.cnv.seg .
cp ~/share/exigens/brain-mets/gistic/absolute/*.lgr.seg .
cp ~/share/exigens/brain-mets/gistic/absolute/*.cnvf.seg .
cp ~/share/exigens/tcga/tcga-luad/gistic/absolute/tcga-luad_pass_absolute.lgr.cntf.cenf.cnvf*.seg .

mkdir -p luad_absolute_bmet-only_cnvf
cp ~/share/exigens/brain-mets/gistic/absolute/luad_absolute_bmet-only_cnvf/{all_lesions.conf_95.txt,amp_genes.conf_95.txt,del_genes.conf_95.txt,scores.gistic} luad_absolute_bmet-only_cnvf/

mkdir -p luad_absolute_prim-only_cnvf
cp ~/share/exigens/brain-mets/gistic/absolute/luad_absolute_prim-only_cnvf/{all_lesions.conf_95.txt,amp_genes.conf_95.txt,del_genes.conf_95.txt,scores.gistic} luad_absolute_prim-only_cnvf/

mkdir -p stage3/luad_absolute_bmet-only_cnvf
cp ~/share/exigens/brain-mets/gistic/absolute/stage3/luad_absolute_bmet-only_cnvf/{all_lesions.conf_95.txt,amp_genes.conf_95.txt,del_genes.conf_95.txt,scores.gistic} stage3/luad_absolute_bmet-only_cnvf/
cp ~/share/exigens/brain-mets/gistic/absolute/stage3/*.{R,m} stage3/

mkdir -p wstage/luad_absolute_bmet-only_cnvf
cp ~/share/exigens/brain-mets/gistic/absolute/wstage/luad_absolute_bmet-only_cnvf/{all_lesions.conf_95.txt,amp_genes.conf_95.txt,del_genes.conf_95.txt,scores.gistic} wstage/luad_absolute_bmet-only_cnvf/
cp ~/share/exigens/brain-mets/gistic/absolute/wstage/*.{R,m} wstage/

mkdir -p wpurity/luad_absolute_bmet-only_cnvf
cp ~/share/exigens/brain-mets/gistic/absolute/wpurity/luad_absolute_bmet-only_cnvf/{all_lesions.conf_95.txt,amp_genes.conf_95.txt,del_genes.conf_95.txt,scores.gistic} wpurity/luad_absolute_bmet-only_cnvf/
cp ~/share/exigens/brain-mets/gistic/absolute/wpurity/*.{R,m} wpurity/

cp -r ~/share/exigens/brain-mets/gistic/absolute/plot .

mkdir -p gpldiff
cp ~/share/exigens/brain-mets/gistic/absolute/gpldiff/*.{R,tsv,pdf} gpldiff/
cp -r ~/share/exigens/brain-mets/gistic/absolute/gpldiff/stan gpldiff/

mkdir -p wpurity/gpldiff
cp ~/share/exigens/brain-mets/gistic/absolute/wpurity/gpldiff/*.{R,tsv} wpurity/gpldiff/

mkdir -p wstage/gpldiff
cp ~/share/exigens/brain-mets/gistic/absolute/wstage/gpldiff/*.{R,tsv} wstage/gpldiff/

mkdir -p stage3/gpldiff
cp ~/share/exigens/brain-mets/gistic/absolute/stage3/gpldiff/*.{R,tsv} stage3/gpldiff

