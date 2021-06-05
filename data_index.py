"""The `data_index` module contains remote data file information about CsTSI and multiple GS gene transcripts."""

CSTSI = "cstsi.fas"
CSGSI = "csgsi.fas"
VVGS = "vvgs.fas"
HSGS = "hsgs.fas"
CGGS = "cggs.fas"
PAGS = "pags.fas"
ZMGS = "zmgs.fas"
FNGS = "fngs.fas"
NTGS = "ntgs.fas"

CSTSI_PROTEIN = "cstsi.protein.fas"
CSGSI_PROTEIN = "csgsi.protein.fas"

DATA_URLS = {
    CSTSI: {
        "url": "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&id=1573295573&extrafeat=null&conwithfeat=on&hide-cdd=on",
    },
    CSGSI: {
        "url": "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&id=161789847&extrafeat=null&conwithfeat=on&hide-cdd=on",
    },
    VVGS: {
        "url": "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&id=1134897&extrafeat=null&conwithfeat=on&hide-cdd=on",
    },
    HSGS: {
        "url": "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&id=45331252&extrafeat=null&conwithfeat=on&hide-cdd=on",
    },
    CGGS: {
        "url": "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&id=57283126&extrafeat=null&conwithfeat=on&hide-cdd=on",
    },
    PAGS: {
        "url": "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&id=160785&extrafeat=null&conwithfeat=on&hide-cdd=on",
    },
    ZMGS: {
        "url": "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&id=434327&extrafeat=null&conwithfeat=on&hide-cdd=on",
    },
    FNGS: {
        "url": "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&id=16505733&extrafeat=null&conwithfeat=on&hide-cdd=on",
    },
    NTGS: {
        "url": "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&id=1419093&extrafeat=null&conwithfeat=on&hide-cdd=on",
    },
    CSTSI_PROTEIN: {
        "url": "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta_cds_aa&id=1573295573&extrafeat=null&conwithfeat=on&hide-cdd=on",
    },
    CSGSI_PROTEIN: {
        "url": "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta_cds_aa&id=161789847&extrafeat=null&conwithfeat=on&hide-cdd=on",
    },
}
