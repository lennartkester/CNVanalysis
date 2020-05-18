#!/bin/bash

#SBATCH -t 04:00:00
#SBATCH --mem=50G
#SBATCH -o CreateReadCountPanelOfNormals.out
#SBATCH -e CreateReadCountPanelOfNormals.err

cd CNV_before_20200320

module load Java/1.8.0_60
export LD_PRELOAD=/usr/lib64/libopenblas-r0.3.3.so
java -jar /hpc/pmc_gen/lkester/bin/gatk-4.1.5.0/gatk-package-4.1.5.0-local.jar CreateReadCountPanelOfNormals \
	-I  PMABM000BPH_PMCRZ158NAS_WXS.counts.hdf5  \
	-I  PMABM000BSM_PMCRZ852MJI_WXS.counts.hdf5  \
	-I  PMABM000BUY_PMCRZ449QJU_WXS.counts.hdf5  \
	-I  PMABM000CBZ_PMCRZ755GXO_WXS.counts.hdf5  \
	-I  PMABM000ABH_PMCRZ763DSH_WXS.counts.hdf5  \
	-I  PMABM000AZR_PMCRZ749VFD_WXS.counts.hdf5  \
	-I  PMABM000CEJ_PMCRZ675NRR_WXS.counts.hdf5  \
	-I  PMABM000CJM_PMCRZ547TLW_WXS.counts.hdf5  \
	-I  PMABM000ACW_PMCRZ795GNG_WXS.counts.hdf5  \
	-I  PMABM000CDH_PMCRZ356NCP_WXS.counts.hdf5  \
	-I  PMABM000BSK_PMCRZ424FYZ_WXS.counts.hdf5  \
	-I  PMABM000CHB_PMCRZ445ZFZ_WXS.counts.hdf5  \
	-I  PMABM000CJV_PMCRZ569XUT_WXS.counts.hdf5  \
	-I  PMABM000BSQ_PMCRZ141WOC_WXS.counts.hdf5  \
	-I  PMABM000BPV_PMCRZ007HML_WXS.counts.hdf5  \
	-I  PMABM000CHD_PMCRZ088QTX_WXS.counts.hdf5  \
	-I  PMABM000CLC_PMCRZ297BYU_WXS.counts.hdf5  \
	-I  PMABM000CLO_PMCRZ281NLW_WXS.counts.hdf5  \
	-I  PMABM000CLE_PMCRZ233PPH_WXS.counts.hdf5  \
	-I  PMABM000CLK_PMCRZ263AEG_WXS.counts.hdf5  \
	-I  PMABM000CLW_PMCRZ325FCI_WXS.counts.hdf5  \
	-I  PMABM000CMG_PMCRZ815LPJ_WXS.counts.hdf5  \
	-I  PMABM000CNA_PMCRZ356VBD_WXS.counts.hdf5  \
	-I  PMABM000CID_PMCRZ042EZR_WXS.counts.hdf5  \
	-I  PMABM000BYP_PMCRZ849SZM_WXS.counts.hdf5  \
	-I  PMABM000CED_PMCRZ942RUE_WXS.counts.hdf5  \
	-I  PMABM000CMI_PMCRZ118ZSU_WXS.counts.hdf5  \
	-I  PMABM000CMW_PMCRZ512UVC_WXS.counts.hdf5  \
	-I  PMABM000CNO_PMCRZ014WRV_WXS.counts.hdf5  \
	-I  PMABM000COQ_PMCRZ222QGW_WXS.counts.hdf5  \
	-I  PMABM000CRY_PMCRZ843VAW_WXS.counts.hdf5  \
	-I  PMABM000CSA_PMCRZ030QFZ_WXS.counts.hdf5  \
	-I  PMABM000CQC_PMCRZ269CQP_WXS.counts.hdf5  \
	-I  PMABM000CPW_PMCRZ019SYB_WXS.counts.hdf5  \
	-I  PMABM000COO_PMCRZ326HBW_WXS.counts.hdf5  \
	-I  PMABM000CNS_PMCRZ648GZF_WXS.counts.hdf5  \
	-I  PMABM000CPK_PMCRZ577GMU_WXS.counts.hdf5  \
	-I  PMABM000CQS_PMCRZ594GFB_WXS.counts.hdf5  \
	-I  PMABM000CRI_PMCRZ735LQU_WXS.counts.hdf5  \
	-I  PMABM000CTM_PMCRZ056RZY_WXS.counts.hdf5  \
	-I  PMABM000CTO_PMCRZ343QMT_WXS.counts.hdf5  \
	-I  PMABM000CSO_PMCRZ686IIT_WXS.counts.hdf5  \
	-I  PMABM000CPS_PMCRZ123MOW_WXS.counts.hdf5  \
	-I  PMABM000CUQ_PMCRZ851QGU_WXS.counts.hdf5  \
	-I  PMABM000CUO_PMCRZ198KLS_WXS.counts.hdf5  \
	-I  PMABM000CVC_PMCRZ092EYP_WXS.counts.hdf5  \
	-I  PMABM000CFM_PMCRZ121JOM_WXS.counts.hdf5  \
	-I  PMABM000COS_PMCRZ602HDZ_WXS.counts.hdf5  \
	-I  PMABM000BZP_PMCRZ178DET_WXS.counts.hdf5  \
	-I  PMABM000CXJ_PMCRZ964YLM_WXS.counts.hdf5  \
	-I  PMABM000CSQ_PMCRZ681ULP_WXS.counts.hdf5  \
	-I  PMABM000CVI_PMCRZ348LWG_WXS.counts.hdf5  \
	-I  PMABM000CVG_PMCRZ154WHT_WXS.counts.hdf5  \
	-I  PMABM000DEB_PMCRZ436GAH_WXS.counts.hdf5  \
	-I  PMABM000CKD_PMCRZ996YUP_WXS.counts.hdf5  \
	-I  PMABM000CWE_PMCRZ709EHX_WXS.counts.hdf5  \
	-I  PMABM000CPM_PMCRZ308TRH_WXS.counts.hdf5  \
	-I  PMABM000CVE_PMCRZ377OLZ_WXS.counts.hdf5  \
	-I  PMABM000CWC_PMCRZ092QJE_WXS.counts.hdf5  \
	-I  PMABM000CWU_PMCRZ593CMR_WXS.counts.hdf5  \
	-I  PMABM000CWZ_PMCRZ022DOS_WXS.counts.hdf5  \
	-I  PMABM000CZT_PMCRZ645EVS_WXS.counts.hdf5  \
	-I  PMABM000CZP_PMCRZ853RIG_WXS.counts.hdf5  \
	-I  PMABM000CXN_PMCRZ141TXI_WXS.counts.hdf5  \
	-I  PMABM000COU_PMCRZ522RFK_WXS.counts.hdf5  \
	-I  PMABM000CYB_PMCRZ689ZPR_WXS.counts.hdf5  \
	-I  PMABM000CYH_PMCRZ100SMD_WXS.counts.hdf5  \
	-I  PMABM000DAB_PMCRZ729GBQ_WXS.counts.hdf5  \
	-I  PMABM000DAJ_PMCRZ004IKD_WXS.counts.hdf5  \
	-I  PMABM000DAL_PMCRZ281LUI_WXS.counts.hdf5  \
	-I  PMABM000DBL_PMCRZ362LKK_WXS.counts.hdf5  \
	-I  PMABM000CYP_PMCRZ393ETF_WXS.counts.hdf5  \
	-I  PMABM000CYN_PMCRZ819TOF_WXS.counts.hdf5  \
	-I  PMABM000DAF_PMCRZ950ZJJ_WXS.counts.hdf5  \
	-I  PMABM000DAN_PMCRZ839YYH_WXS.counts.hdf5  \
	-I  PMABM000DED_PMCRZ333GHV_WXS.counts.hdf5  \
	-I  PMABM000DDI_PMCRZ525ELX_WXS.counts.hdf5  \
	-I  PMABM000DDZ_PMCRZ493CFK_WXS.counts.hdf5  \
	-I  PMABM000DFD_PMCRZ800VMF_WXS.counts.hdf5  \
	-I  PMABM000DAH_PMCRZ067SCE_WXS.counts.hdf5  \
	-I  PMABM000DBH_PMCRZ367SVS_WXS.counts.hdf5  \
	-I  PMABM000DCD_PMCRZ979YGB_WXS.counts.hdf5  \
	-I  PMABM000DNH_PMCRZ077EBM_WXS.counts.hdf5  \
	-I  PMABM000DPJ_PMCRZ879LMM_WXS.counts.hdf5  \
	-I  PMABM000DPH_PMCRZ539IBD_WXS.counts.hdf5  \
	-I  PMABM000DPP_PMCRZ101KQT_WXS.counts.hdf5  \
	-I  PMABM000DDT_PMCRZ405WQX_WXS.counts.hdf5  \
	-I  PMABM000DMX_PMCRZ099ZSP_WXS.counts.hdf5  \
	-I  PMABM000DNV_PMCRZ140RMX_WXS.counts.hdf5  \
	-I  PMABM000DOL_PMCRZ687VQV_WXS.counts.hdf5  \
	-I  PMABM000DPD_PMCRZ818RGV_WXS.counts.hdf5  \
	-I  PMRBM000AFX_PMCRZ443PUI_WXS.counts.hdf5  \
	-I  PMABM000DGJ_PMCRZ241MLV_WXS.counts.hdf5  \
	-I  PMABM000DPB_PMCRZ484LRQ_WXS.counts.hdf5  \
	-I  PMABM000DTG_PMCRZ108UMM_WXS.counts.hdf5  \
	-I  PMABM000DQZ_PMCRZ354HJD_WXS.counts.hdf5  \
	-I  PMABM000DQR_PMCRZ368DMG_WXS.counts.hdf5  \
	-I  PMABM000DQN_PMCRZ450XAU_WXS.counts.hdf5  \
	-I  PMABM000DJN_PMCRZ004SFH_WXS.counts.hdf5  \
	-I  PMABM000DIF_PMCRZ023OIL_WXS.counts.hdf5  \
	-I  PMABM000DNL_PMCRZ258ESE_WXS.counts.hdf5  \
	-I  PMABM000DZI_PMCRZ965NKR_WXS.counts.hdf5  \
	-I  PMABM000CXB_PMCRZ196DVH_WXS.counts.hdf5  \
	-I  PMABM000DWA_PMCRZ452TUK_WXS.counts.hdf5  \
	-I  PMABM000DZW_PMCRZ716LSL_WXS.counts.hdf5  \
	-I  PMABM000EPQ_PMCRZ743ZTV_WXS.counts.hdf5  \
	-I  PMABM000EMG_PMCRZ231KTX_WXS.counts.hdf5  \
	-I  PMABM000CBP_PMCRZ220NXR_WXS.counts.hdf5  \
	-I  PMABM000DEF_PMCRZ007XQX_WXS.counts.hdf5  \
	-I  PMABM000EQE_PMCRZ787WQO_WXS.counts.hdf5  \
	-I  PMABM000EYC_PMCRZ507BJB_WXS.counts.hdf5  \
	-I  PMABM000FDK_PMCRZ576TGZ_WXS.counts.hdf5  \
	-I  PMABM000FCT_PMCRZ419ILU_WXS.counts.hdf5  \
	-I  PMABM000FDU_PMCRZ206FJU_WXS.counts.hdf5  \
	-I  PMABM000FFR_PMCRZ206EHX_WXS.counts.hdf5  \
	-I  PMABM000FGW_PMCRZ504PHS_WXS.counts.hdf5  \
	-I  PMABM000FNM_PMCRZ670XQT_WXS.counts.hdf5  \
	-I  PMABM000FNA_PMCRZ003OKM_WXS.counts.hdf5  \
	-I  PMABM000FFP_PMCRZ495KJC_WXS.counts.hdf5  \
	-I  PMABM000DKT_PMCRZ069GHW_WXS.counts.hdf5  \
	-I  PMABM000FQS_PMCRZ724XCV_WXS.counts.hdf5  \
	-I  PMABM000FRJ_PMCRZ162NAT_WXS.counts.hdf5  \
	-I  PMABM000FTO_PMCRZ710QHZ_WXS.counts.hdf5  \
	-I  PMABM000FTG_PMCRZ590NQO_WXS.counts.hdf5  \
	-I  PMABM000FVZ_PMCRZ295VUX_WXS.counts.hdf5  \
	-I  PMABM000FXF_PMCRZ804UXD_WXS.counts.hdf5  \
	-I  PMABM000FUL_PMCRZ238QIN_WXS.counts.hdf5  \
	-I  PMABM000FYP_PMCRZ764OYX_WXS.counts.hdf5  \
	-I  PMABM000FWB_PMCRZ540VYQ_WXS.counts.hdf5  \
	-I  PMABM000FVY_PMCRZ438OHE_WXS.counts.hdf5  \
	-I  PMABM000FKB_PMCRZ049JWE_WXS.counts.hdf5  \
	-I  PMABM000FPV_PMCRZ700NUI_WXS.counts.hdf5  \
	-I  PMABM000FYM_PMCRZ494BVB_WXS.counts.hdf5  \
	-I  PMABM000FZT_PMCRZ741INS_WXS.counts.hdf5  \
	-I  PMABM000FZW_PMCRZ453AUE_WXS.counts.hdf5  \
	-I  PMABM000GBT_PMCRZ622PSW_WXS.counts.hdf5  \
	-I  PMABM000GER_PMCRZ864JXO_WXS.counts.hdf5  \
	-I  PMABM000GGF_PMCRZ761QSP_WXS.counts.hdf5  \
	-I  PMABM000GLE_PMCRZ229HTC_WXS.counts.hdf5  \
	-I  PMABM000GIY_PMCRZ671ISA_WXS.counts.hdf5  \
	-I  PMABM000GGG_PMCRZ099MSG_WXS.counts.hdf5  \
	-I  PMABM000FMS_PMCRZ845SCE_WXS.counts.hdf5  \
	-I  PMABM000GGB_PMCRZ952QFH_WXS.counts.hdf5  \
	--annotated-intervals /hpc/pmc_gen/lkester/CNV_calling/reference/MedExome_hg38_capture_targets_padded250_annotated.interval_list \
	--minimum-interval-median-percentile 10.0 \
	-O 20200428_PoN_padded250_mimp10.hdf5

