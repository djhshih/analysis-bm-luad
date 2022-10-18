library(io);
library(dplyr);

options(stringsAsFactors=FALSE);

clin <- qread("~/exigens/brain-mets/annot/patient-info_stage2.rds");
pheno.all <- qread("~/exigens/brain-mets/annot/sample-info_wes_stage3_pass_luad.tsv");
patients <- as.character(unique(pheno.all$clinical_id));
clin <- clin[match(patients, clin$clinical_id), ];

genes <- c("YAP1", "MYC", "CDKN2A");
directions <- c(1, 1, -1);

amp.cut <- 8;
del.cut <- 0.5;

pheno.all$pick[pheno.all$sample_id %in% c("LS-016-M1", "LS-016-P-2")] <- TRUE;

#pheno.sel <- filter(pheno.all, pick);
#pheno <- pheno.sel;

pheno <- pheno.all;

cn <- qread("~/exigens/brain-mets/comut/brain-mets_pass_luad_absolute-1-4_gene_corrected_CN.txt", type="tsv");
cn <- as.matrix(cn);
cn.sel <- cn[genes, ];
samples <- colnames(cn);

cna.samples <- mapply(
	function(gene, direction) {
		if (direction > 0) {
			samples[ cn[gene, ] > amp.cut ]
		} else {
			samples[ cn[gene, ] < del.cut ]
		}
	},
	genes, directions,
	SIMPLIFY=FALSE
);

qwrite(cna.samples, "cna-samples.rds");

cna.patients <- lapply(cna.samples,
	function(samples) {
		patients <- unique(pheno$clinical_id[match(samples, pheno$sample_id)]);
		y <- lapply(patients,
			function(patient) {
				as.character(pheno$sample_id[which(pheno$clinical_id == patient)])
			}
		);
		names(y) <- patients;
		y
	}
);

qwrite(cna.patients, "cna-patients.rds");

# samples from patients with any sample that harbours a CNA of target
cna.samples.any <- unname(unique(unlist(cna.patients)));
qwrite(cna.samples.any, "cna-samples-any.vtr");

