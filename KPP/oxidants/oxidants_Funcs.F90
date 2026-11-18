!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: oxidants_Funcs
!
! !DESCRIPTION: Module containing routines for passing data from GEOS-Chem
!  to the prescribed oxidants mechanism solver code.
!\\
!\\
! !INTERFACE:
!
MODULE oxidants_Funcs
!
! !USES:
!
  USE gckpp_Precision
  USE gckpp_Parameters
  USE gckpp_Global
  USE Precision_Mod,   ONLY : fp
  USE rateLawUtilFuncs

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: oxidants_ConvertKgToMolecCm3
  PUBLIC :: oxidants_ComputeRateConstants
  PUBLIC :: oxidants_ConvertMolecCm3ToKg
  PUBLIC :: oxidants_InitKPPFuncs
  PUBLIC :: oxidants_CleanupKPPFuncs
  PUBLIC :: oxidants_Get_COfromCH4_Flux
  PUBLIC :: oxidants_Get_COfromNMVOC_Flux
  PUBLIC :: oxidants_Get_CH4byOH_Loss
  PUBLIC :: oxidants_Get_CH4byCl_Loss
  PUBLIC :: oxidants_Get_CH4byO1D_Loss
  PUBLIC :: oxidants_Get_CH4bySoil_Loss
  PUBLIC :: oxidants_Get_CH4SoilUptakeRate
  PUBLIC :: oxidants_Get_CObyOH_Loss
  PUBLIC :: oxidants_Get_TaggedCH4Loss
  PUBLIC :: GC_OHCO
!
! !PUBLIC TYPES:
!
  ! Species ID flags
  INTEGER  :: id_CH4
  INTEGER  :: id_CO
  INTEGER, PARAMETER :: N_BASE_TAGGED_CH4 = 39
  INTEGER, PARAMETER :: N_TAGGED_CH4 = 389

  CHARACTER(LEN=11), PARAMETER :: TAGGED_CH4_NAMES(N_TAGGED_CH4) = (/ CHARACTER(LEN=11) :: &
       'CH4_OIL', 'CH4_GAS', 'CH4_COL', &
       'CH4_LIV', 'CH4_LDF', 'CH4_WST', &
       'CH4_RIC', 'CH4_OTA', 'CH4_BBN', &
       'CH4_WTL', 'CH4_SEE', 'CH4_LAK', &
       'CH4_TER', 'CH4_RES', 'CH4_USA', &
       'CH4_CAN', 'CH4_MEX', 'CH4_CAC', &
       'CH4_BRA', 'CH4_NSA', 'CH4_SCN', &
       'CH4_WEU', 'CH4_EEU', 'CH4_RUS', &
       'CH4_MDE', 'CH4_NAF', 'CH4_WAF', &
       'CH4_CAF', 'CH4_EAF', 'CH4_SIO', &
       'CH4_SAC', 'CH4_CAS', 'CH4_IND', &
       'CH4_SAS', 'CH4_CHN', 'CH4_EAS', &
       'CH4_SEA', 'CH4_OCE', 'CH4_ROW', &
       'CH4_USA_OIL', 'CH4_USA_GAS', 'CH4_USA_COL', &
       'CH4_USA_LIV', 'CH4_USA_LDF', 'CH4_USA_WST', &
       'CH4_USA_RIC', 'CH4_USA_OTA', 'CH4_USA_BBN', &
       'CH4_USA_WTL', 'CH4_USA_SEE', 'CH4_USA_LAK', &
       'CH4_USA_TER', 'CH4_USA_RES', 'CH4_CAN_OIL', &
       'CH4_CAN_GAS', 'CH4_CAN_COL', 'CH4_CAN_LIV', &
       'CH4_CAN_LDF', 'CH4_CAN_WST', 'CH4_CAN_RIC', &
       'CH4_CAN_OTA', 'CH4_CAN_BBN', 'CH4_CAN_WTL', &
       'CH4_CAN_SEE', 'CH4_CAN_LAK', 'CH4_CAN_TER', &
       'CH4_CAN_RES', 'CH4_MEX_OIL', 'CH4_MEX_GAS', &
       'CH4_MEX_COL', 'CH4_MEX_LIV', 'CH4_MEX_LDF', &
       'CH4_MEX_WST', 'CH4_MEX_RIC', 'CH4_MEX_OTA', &
       'CH4_MEX_BBN', 'CH4_MEX_WTL', 'CH4_MEX_SEE', &
       'CH4_MEX_LAK', 'CH4_MEX_TER', 'CH4_MEX_RES', &
       'CH4_CAC_OIL', 'CH4_CAC_GAS', 'CH4_CAC_COL', &
       'CH4_CAC_LIV', 'CH4_CAC_LDF', 'CH4_CAC_WST', &
       'CH4_CAC_RIC', 'CH4_CAC_OTA', 'CH4_CAC_BBN', &
       'CH4_CAC_WTL', 'CH4_CAC_SEE', 'CH4_CAC_LAK', &
       'CH4_CAC_TER', 'CH4_CAC_RES', 'CH4_BRA_OIL', &
       'CH4_BRA_GAS', 'CH4_BRA_COL', 'CH4_BRA_LIV', &
       'CH4_BRA_LDF', 'CH4_BRA_WST', 'CH4_BRA_RIC', &
       'CH4_BRA_OTA', 'CH4_BRA_BBN', 'CH4_BRA_WTL', &
       'CH4_BRA_SEE', 'CH4_BRA_LAK', 'CH4_BRA_TER', &
       'CH4_BRA_RES', 'CH4_NSA_OIL', 'CH4_NSA_GAS', &
       'CH4_NSA_COL', 'CH4_NSA_LIV', 'CH4_NSA_LDF', &
       'CH4_NSA_WST', 'CH4_NSA_RIC', 'CH4_NSA_OTA', &
       'CH4_NSA_BBN', 'CH4_NSA_WTL', 'CH4_NSA_SEE', &
       'CH4_NSA_LAK', 'CH4_NSA_TER', 'CH4_NSA_RES', &
       'CH4_SCN_OIL', 'CH4_SCN_GAS', 'CH4_SCN_COL', &
       'CH4_SCN_LIV', 'CH4_SCN_LDF', 'CH4_SCN_WST', &
       'CH4_SCN_RIC', 'CH4_SCN_OTA', 'CH4_SCN_BBN', &
       'CH4_SCN_WTL', 'CH4_SCN_SEE', 'CH4_SCN_LAK', &
       'CH4_SCN_TER', 'CH4_SCN_RES', 'CH4_WEU_OIL', &
       'CH4_WEU_GAS', 'CH4_WEU_COL', 'CH4_WEU_LIV', &
       'CH4_WEU_LDF', 'CH4_WEU_WST', 'CH4_WEU_RIC', &
       'CH4_WEU_OTA', 'CH4_WEU_BBN', 'CH4_WEU_WTL', &
       'CH4_WEU_SEE', 'CH4_WEU_LAK', 'CH4_WEU_TER', &
       'CH4_WEU_RES', 'CH4_EEU_OIL', 'CH4_EEU_GAS', &
       'CH4_EEU_COL', 'CH4_EEU_LIV', 'CH4_EEU_LDF', &
       'CH4_EEU_WST', 'CH4_EEU_RIC', 'CH4_EEU_OTA', &
       'CH4_EEU_BBN', 'CH4_EEU_WTL', 'CH4_EEU_SEE', &
       'CH4_EEU_LAK', 'CH4_EEU_TER', 'CH4_EEU_RES', &
       'CH4_RUS_OIL', 'CH4_RUS_GAS', 'CH4_RUS_COL', &
       'CH4_RUS_LIV', 'CH4_RUS_LDF', 'CH4_RUS_WST', &
       'CH4_RUS_RIC', 'CH4_RUS_OTA', 'CH4_RUS_BBN', &
       'CH4_RUS_WTL', 'CH4_RUS_SEE', 'CH4_RUS_LAK', &
       'CH4_RUS_TER', 'CH4_RUS_RES', 'CH4_MDE_OIL', &
       'CH4_MDE_GAS', 'CH4_MDE_COL', 'CH4_MDE_LIV', &
       'CH4_MDE_LDF', 'CH4_MDE_WST', 'CH4_MDE_RIC', &
       'CH4_MDE_OTA', 'CH4_MDE_BBN', 'CH4_MDE_WTL', &
       'CH4_MDE_SEE', 'CH4_MDE_LAK', 'CH4_MDE_TER', &
       'CH4_MDE_RES', 'CH4_NAF_OIL', 'CH4_NAF_GAS', &
       'CH4_NAF_COL', 'CH4_NAF_LIV', 'CH4_NAF_LDF', &
       'CH4_NAF_WST', 'CH4_NAF_RIC', 'CH4_NAF_OTA', &
       'CH4_NAF_BBN', 'CH4_NAF_WTL', 'CH4_NAF_SEE', &
       'CH4_NAF_LAK', 'CH4_NAF_TER', 'CH4_NAF_RES', &
       'CH4_WAF_OIL', 'CH4_WAF_GAS', 'CH4_WAF_COL', &
       'CH4_WAF_LIV', 'CH4_WAF_LDF', 'CH4_WAF_WST', &
       'CH4_WAF_RIC', 'CH4_WAF_OTA', 'CH4_WAF_BBN', &
       'CH4_WAF_WTL', 'CH4_WAF_SEE', 'CH4_WAF_LAK', &
       'CH4_WAF_TER', 'CH4_WAF_RES', 'CH4_CAF_OIL', &
       'CH4_CAF_GAS', 'CH4_CAF_COL', 'CH4_CAF_LIV', &
       'CH4_CAF_LDF', 'CH4_CAF_WST', 'CH4_CAF_RIC', &
       'CH4_CAF_OTA', 'CH4_CAF_BBN', 'CH4_CAF_WTL', &
       'CH4_CAF_SEE', 'CH4_CAF_LAK', 'CH4_CAF_TER', &
       'CH4_CAF_RES', 'CH4_EAF_OIL', 'CH4_EAF_GAS', &
       'CH4_EAF_COL', 'CH4_EAF_LIV', 'CH4_EAF_LDF', &
       'CH4_EAF_WST', 'CH4_EAF_RIC', 'CH4_EAF_OTA', &
       'CH4_EAF_BBN', 'CH4_EAF_WTL', 'CH4_EAF_SEE', &
       'CH4_EAF_LAK', 'CH4_EAF_TER', 'CH4_EAF_RES', &
       'CH4_SIO_OIL', 'CH4_SIO_GAS', 'CH4_SIO_COL', &
       'CH4_SIO_LIV', 'CH4_SIO_LDF', 'CH4_SIO_WST', &
       'CH4_SIO_RIC', 'CH4_SIO_OTA', 'CH4_SIO_BBN', &
       'CH4_SIO_WTL', 'CH4_SIO_SEE', 'CH4_SIO_LAK', &
       'CH4_SIO_TER', 'CH4_SIO_RES', 'CH4_SAC_OIL', &
       'CH4_SAC_GAS', 'CH4_SAC_COL', 'CH4_SAC_LIV', &
       'CH4_SAC_LDF', 'CH4_SAC_WST', 'CH4_SAC_RIC', &
       'CH4_SAC_OTA', 'CH4_SAC_BBN', 'CH4_SAC_WTL', &
       'CH4_SAC_SEE', 'CH4_SAC_LAK', 'CH4_SAC_TER', &
       'CH4_SAC_RES', 'CH4_CAS_OIL', 'CH4_CAS_GAS', &
       'CH4_CAS_COL', 'CH4_CAS_LIV', 'CH4_CAS_LDF', &
       'CH4_CAS_WST', 'CH4_CAS_RIC', 'CH4_CAS_OTA', &
       'CH4_CAS_BBN', 'CH4_CAS_WTL', 'CH4_CAS_SEE', &
       'CH4_CAS_LAK', 'CH4_CAS_TER', 'CH4_CAS_RES', &
       'CH4_IND_OIL', 'CH4_IND_GAS', 'CH4_IND_COL', &
       'CH4_IND_LIV', 'CH4_IND_LDF', 'CH4_IND_WST', &
       'CH4_IND_RIC', 'CH4_IND_OTA', 'CH4_IND_BBN', &
       'CH4_IND_WTL', 'CH4_IND_SEE', 'CH4_IND_LAK', &
       'CH4_IND_TER', 'CH4_IND_RES', 'CH4_SAS_OIL', &
       'CH4_SAS_GAS', 'CH4_SAS_COL', 'CH4_SAS_LIV', &
       'CH4_SAS_LDF', 'CH4_SAS_WST', 'CH4_SAS_RIC', &
       'CH4_SAS_OTA', 'CH4_SAS_BBN', 'CH4_SAS_WTL', &
       'CH4_SAS_SEE', 'CH4_SAS_LAK', 'CH4_SAS_TER', &
       'CH4_SAS_RES', 'CH4_CHN_OIL', 'CH4_CHN_GAS', &
       'CH4_CHN_COL', 'CH4_CHN_LIV', 'CH4_CHN_LDF', &
       'CH4_CHN_WST', 'CH4_CHN_RIC', 'CH4_CHN_OTA', &
       'CH4_CHN_BBN', 'CH4_CHN_WTL', 'CH4_CHN_SEE', &
       'CH4_CHN_LAK', 'CH4_CHN_TER', 'CH4_CHN_RES', &
       'CH4_EAS_OIL', 'CH4_EAS_GAS', 'CH4_EAS_COL', &
       'CH4_EAS_LIV', 'CH4_EAS_LDF', 'CH4_EAS_WST', &
       'CH4_EAS_RIC', 'CH4_EAS_OTA', 'CH4_EAS_BBN', &
       'CH4_EAS_WTL', 'CH4_EAS_SEE', 'CH4_EAS_LAK', &
       'CH4_EAS_TER', 'CH4_EAS_RES', 'CH4_SEA_OIL', &
       'CH4_SEA_GAS', 'CH4_SEA_COL', 'CH4_SEA_LIV', &
       'CH4_SEA_LDF', 'CH4_SEA_WST', 'CH4_SEA_RIC', &
       'CH4_SEA_OTA', 'CH4_SEA_BBN', 'CH4_SEA_WTL', &
       'CH4_SEA_SEE', 'CH4_SEA_LAK', 'CH4_SEA_TER', &
       'CH4_SEA_RES', 'CH4_OCE_OIL', 'CH4_OCE_GAS', &
       'CH4_OCE_COL', 'CH4_OCE_LIV', 'CH4_OCE_LDF', &
       'CH4_OCE_WST', 'CH4_OCE_RIC', 'CH4_OCE_OTA', &
       'CH4_OCE_BBN', 'CH4_OCE_WTL', 'CH4_OCE_SEE', &
       'CH4_OCE_LAK', 'CH4_OCE_TER', 'CH4_OCE_RES', &
       'CH4_ROW_OIL', 'CH4_ROW_GAS', 'CH4_ROW_COL', &
       'CH4_ROW_LIV', 'CH4_ROW_LDF', 'CH4_ROW_WST', &
       'CH4_ROW_RIC', 'CH4_ROW_OTA', 'CH4_ROW_BBN', &
       'CH4_ROW_WTL', 'CH4_ROW_SEE', 'CH4_ROW_LAK', &
       'CH4_ROW_TER', 'CH4_ROW_RES' /)
  INTEGER :: id_CH4_tagged(N_TAGGED_CH4)

  INTEGER, PARAMETER :: ind_CH4_tagged_base(N_BASE_TAGGED_CH4) = (/        &
       ind_CH4_OIL, ind_CH4_GAS, ind_CH4_COL, &
       ind_CH4_LIV, ind_CH4_LDF, ind_CH4_WST, &
       ind_CH4_RIC, ind_CH4_OTA, ind_CH4_BBN, &
       ind_CH4_WTL, ind_CH4_SEE, ind_CH4_LAK, &
       ind_CH4_TER, ind_CH4_RES, &
       ind_CH4_USA, ind_CH4_CAN, ind_CH4_MEX, &
       ind_CH4_CAC, ind_CH4_BRA, ind_CH4_NSA, &
       ind_CH4_SCN, ind_CH4_WEU, ind_CH4_EEU, &
       ind_CH4_RUS, ind_CH4_MDE, ind_CH4_NAF, &
       ind_CH4_WAF, ind_CH4_CAF, ind_CH4_EAF, &
       ind_CH4_SIO, ind_CH4_SAC, ind_CH4_CAS, &
       ind_CH4_IND, ind_CH4_SAS, ind_CH4_CHN, &
       ind_CH4_EAS, ind_CH4_SEA, ind_CH4_OCE, &
       ind_CH4_ROW /)
  INTEGER, PARAMETER :: ind_LCH4_tagged_byOH(N_BASE_TAGGED_CH4) = (/        &
       ind_LCH4_OILbyOH, ind_LCH4_GASbyOH, ind_LCH4_COLbyOH, &
       ind_LCH4_LIVbyOH, ind_LCH4_LDFbyOH, ind_LCH4_WSTbyOH, &
       ind_LCH4_RICbyOH, ind_LCH4_OTAbyOH, ind_LCH4_BBNbyOH, &
       ind_LCH4_WTLbyOH, ind_LCH4_SEEbyOH, ind_LCH4_LAKbyOH, &
       ind_LCH4_TERbyOH, ind_LCH4_RESbyOH, &
       ind_LCH4_USAbyOH, ind_LCH4_CANbyOH, ind_LCH4_MEXbyOH, &
       ind_LCH4_CACbyOH, ind_LCH4_BRAbyOH, ind_LCH4_NSAbyOH, &
       ind_LCH4_SCNbyOH, ind_LCH4_WEUbyOH, ind_LCH4_EEUbyOH, &
       ind_LCH4_RUSbyOH, ind_LCH4_MDEbyOH, ind_LCH4_NAFbyOH, &
       ind_LCH4_WAFbyOH, ind_LCH4_CAFbyOH, ind_LCH4_EAFbyOH, &
       ind_LCH4_SIObyOH, ind_LCH4_SACbyOH, ind_LCH4_CASbyOH, &
       ind_LCH4_INDbyOH, ind_LCH4_SASbyOH, ind_LCH4_CHNbyOH, &
       ind_LCH4_EASbyOH, ind_LCH4_SEAbyOH, ind_LCH4_OCEbyOH, &
       ind_LCH4_ROWbyOH /)
  INTEGER, PARAMETER :: ind_LCH4_tagged_byCl(N_BASE_TAGGED_CH4) = (/        &
       ind_LCH4_OILbyCl, ind_LCH4_GASbyCl, ind_LCH4_COLbyCl, &
       ind_LCH4_LIVbyCl, ind_LCH4_LDFbyCl, ind_LCH4_WSTbyCl, &
       ind_LCH4_RICbyCl, ind_LCH4_OTAbyCl, ind_LCH4_BBNbyCl, &
       ind_LCH4_WTLbyCl, ind_LCH4_SEEbyCl, ind_LCH4_LAKbyCl, &
       ind_LCH4_TERbyCl, ind_LCH4_RESbyCl, &
       ind_LCH4_USAbyCl, ind_LCH4_CANbyCl, ind_LCH4_MEXbyCl, &
       ind_LCH4_CACbyCl, ind_LCH4_BRAbyCl, ind_LCH4_NSAbyCl, &
       ind_LCH4_SCNbyCl, ind_LCH4_WEUbyCl, ind_LCH4_EEUbyCl, &
       ind_LCH4_RUSbyCl, ind_LCH4_MDEbyCl, ind_LCH4_NAFbyCl, &
       ind_LCH4_WAFbyCl, ind_LCH4_CAFbyCl, ind_LCH4_EAFbyCl, &
       ind_LCH4_SIObyCl, ind_LCH4_SACbyCl, ind_LCH4_CASbyCl, &
       ind_LCH4_INDbyCl, ind_LCH4_SASbyCl, ind_LCH4_CHNbyCl, &
       ind_LCH4_EASbyCl, ind_LCH4_SEAbyCl, ind_LCH4_OCEbyCl, &
       ind_LCH4_ROWbyCl /)
  INTEGER, PARAMETER :: ind_LCH4_tagged_byO1D(N_BASE_TAGGED_CH4) = (/       &
       ind_LCH4_OILbyO1D, ind_LCH4_GASbyO1D, ind_LCH4_COLbyO1D, &
       ind_LCH4_LIVbyO1D, ind_LCH4_LDFbyO1D, ind_LCH4_WSTbyO1D, &
       ind_LCH4_RICbyO1D, ind_LCH4_OTAbyO1D, ind_LCH4_BBNbyO1D, &
       ind_LCH4_WTLbyO1D, ind_LCH4_SEEbyO1D, ind_LCH4_LAKbyO1D, &
       ind_LCH4_TERbyO1D, ind_LCH4_RESbyO1D, &
       ind_LCH4_USAbyO1D, ind_LCH4_CANbyO1D, ind_LCH4_MEXbyO1D, &
       ind_LCH4_CACbyO1D, ind_LCH4_BRAbyO1D, ind_LCH4_NSAbyO1D, &
       ind_LCH4_SCNbyO1D, ind_LCH4_WEUbyO1D, ind_LCH4_EEUbyO1D, &
       ind_LCH4_RUSbyO1D, ind_LCH4_MDEbyO1D, ind_LCH4_NAFbyO1D, &
       ind_LCH4_WAFbyO1D, ind_LCH4_CAFbyO1D, ind_LCH4_EAFbyO1D, &
       ind_LCH4_SIObyO1D, ind_LCH4_SACbyO1D, ind_LCH4_CASbyO1D, &
       ind_LCH4_INDbyO1D, ind_LCH4_SASbyO1D, ind_LCH4_CHNbyO1D, &
       ind_LCH4_EASbyO1D, ind_LCH4_SEAbyO1D, ind_LCH4_OCEbyO1D, &
       ind_LCH4_ROWbyO1D /)
  INTEGER, PARAMETER :: ind_LCH4_tagged_bySoil(N_BASE_TAGGED_CH4) = (/      &
       ind_LCH4_OILbySoil, ind_LCH4_GASbySoil, ind_LCH4_COLbySoil, &
       ind_LCH4_LIVbySoil, ind_LCH4_LDFbySoil, ind_LCH4_WSTbySoil, &
       ind_LCH4_RICbySoil, ind_LCH4_OTAbySoil, ind_LCH4_BBNbySoil, &
       ind_LCH4_WTLbySoil, ind_LCH4_SEEbySoil, ind_LCH4_LAKbySoil, &
       ind_LCH4_TERbySoil, ind_LCH4_RESbySoil, &
       ind_LCH4_USAbySoil, ind_LCH4_CANbySoil, ind_LCH4_MEXbySoil, &
       ind_LCH4_CACbySoil, ind_LCH4_BRAbySoil, ind_LCH4_NSAbySoil, &
       ind_LCH4_SCNbySoil, ind_LCH4_WEUbySoil, ind_LCH4_EEUbySoil, &
       ind_LCH4_RUSbySoil, ind_LCH4_MDEbySoil, ind_LCH4_NAFbySoil, &
       ind_LCH4_WAFbySoil, ind_LCH4_CAFbySoil, ind_LCH4_EAFbySoil, &
       ind_LCH4_SIObySoil, ind_LCH4_SACbySoil, ind_LCH4_CASbySoil, &
       ind_LCH4_INDbySoil, ind_LCH4_SASbySoil, ind_LCH4_CHNbySoil, &
       ind_LCH4_EASbySoil, ind_LCH4_SEAbySoil, ind_LCH4_OCEbySoil, &
       ind_LCH4_ROWbySoil /)

#ifdef JACOBIAN
  ! Jacobian CH4 tracers
  INTEGER              :: numJacobianTracers
  INTEGER, ALLOCATABLE :: JacobianIDs(:)
#endif

  ! Kg species / molec species
  REAL(fp) :: xnumol_CH4
  REAL(fp) :: xnumol_CO

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: oxidants_get_taggedch4_kpp_id
!
! !DESCRIPTION: Returns the KPP species index for a tagged CH4 tracer.
!  Region-by-sector tracers are generated in contiguous groups of five.
!\
!\
! !INTERFACE:
!
  PURE INTEGER FUNCTION oxidants_Get_TaggedCH4KppId( tag ) RESULT ( KppId )
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN) :: tag
!EOP
!------------------------------------------------------------------------------
!BOC

    KppId = 0
    IF ( tag >= 1 .and. tag <= N_BASE_TAGGED_CH4 ) THEN
       KppId = ind_CH4_tagged_base(tag)
    ELSE IF ( tag > N_BASE_TAGGED_CH4 .and. tag <= N_TAGGED_CH4 ) THEN
       KppId = ind_CH4_USA_OIL + ( 5 * ( tag - N_BASE_TAGGED_CH4 - 1 ) )
    ENDIF

  END FUNCTION oxidants_Get_TaggedCH4KppId
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: oxidants_ConvertKgToMolecCm3
!
! !DESCRIPTION: Converts species from kg to molec/cm3 and stores into
!  the "C" concentration array used by the KPP-generated solver code.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE oxidants_ConvertKgToMolecCm3( I, J, L, State_Chm, State_Met )
!
! !USES:
!
    USE Species_Mod,   ONLY : SpcConc
    USE State_Chm_Mod, ONLY : ChmState
    USE State_Met_Mod, ONLY : MetState
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN) :: I, J, L     ! Grid box indices
    TYPE(ChmState), INTENT(IN) :: State_Chm   ! Chemistry State object
    TYPE(MetState), INTENT(IN) :: State_Met   ! Meterorology State object
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                :: N
    REAL(fp)               :: airvol_cm3

    ! Pointers
    TYPE(SpcConc), POINTER :: Spc(:)

    !========================================================================
    ! oxidants_ConvertKgToMolecCm3 begins here!
    !========================================================================

    ! Point to species array
    Spc => State_Chm%Species

    ! Grid box volume [cm3]
    airvol_cm3 = State_Met%AIRVOL(I,J,L) * 1.0e+6_fp

    ! Convert species kg to molec/cm3 and store in the C array
    ! (skip if species are not present)
    IF ( id_CH4 > 0 ) THEN
       C(ind_CH4) = Spc(id_CH4)%Conc(I,J,L) * xnumol_CH4 / airvol_cm3
    ENDIF

    IF ( id_CO > 0 ) THEN
       C(ind_CO)  = Spc(id_CO)%Conc(I,J,L)  * xnumol_CO  / airvol_cm3
    ENDIF

    DO N = 1, N_TAGGED_CH4
       IF ( id_CH4_tagged(N) > 0 ) THEN
          C(oxidants_Get_TaggedCH4KppId(N)) =                               &
               Spc(id_CH4_tagged(N))%Conc(I,J,L) * xnumol_CH4 / airvol_cm3
       ENDIF
    ENDDO

#ifdef JACOBIAN
    ! Do the same for Jacobian CH4 tracers, if any
    IF ( numJacobianTracers > 0 ) THEN
       DO N = 1, numJacobianTracers
          C(JacobianIDs(N)) = Spc(JacobianIDs(N))%Conc(I,J,L) * xnumol_CH4 / airvol_cm3
       ENDDO
    ENDIF
#endif

    ! Initialize placeholder species to 1 molec/cm3
    C(ind_DummyCH4)   = 1.0_dp
    C(ind_DummyNMVOC) = 1.0_dp
    C(ind_DummyCO)    = 0.0_dp

    ! Initialize fixed species to 1 molec/cm3
    ! These will later be set to values read via HEMCO
    C(ind_FixedCl) = 1.0_dp
    C(ind_FixedOH) = 1.0_dp
    C(ind_FixedO1D) = 1.0_dp

    ! Free pointer
    Spc => NULL()

  END SUBROUTINE oxidants_ConvertKgToMolecCm3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: oxidants_ComputeRateConstants
!
! !DESCRIPTION: Computes the rate constants used in the carbon chemical
!  mechanism.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE oxidants_ComputeRateConstants(                                    &
             I,             J,                 L,                            &
             ConcClMnd,     ConcOHMnd,         ConcO1Dmnd,                   &
             OHdiurnalFac,  PCO_fr_CH4_use,                                   &
             PCO_fr_CH4,    PCO_fr_NMVOC_use,  PCO_fr_NMVOC,                 &
             CH4_soil_k,    dtChem,            State_Chm,                    &
             State_Met                                                      )
!
! !USES:
!
    USE PhysConstants, ONLY : AVO, AIRMW
    USE State_Chm_Mod, ONLY : ChmState
    USE State_Met_Mod, ONLY : MetState
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN) :: I, J, L
    REAL(fp),       INTENT(IN) :: ConcClMnd       ! Cl conc [molec/cm3]
    REAL(fp),       INTENT(IN) :: ConcOHmnd       ! OH conc [molec/cm3]
    REAL(fp),       INTENT(IN) :: ConcO1Dmnd      ! O(1D) conc [molec/cm3]
    REAL(fp),       INTENT(IN) :: OHdiurnalFac    ! OH diurnal scale factor [1]
    LOGICAL,        INTENT(IN) :: PCO_fr_CH4_use  ! Use P(CO) fr CH4? [T/F]
    REAL(fp),       INTENT(IN) :: PCO_fr_CH4      ! P(CO) fr CH4 [molec/cm3/s]
    LOGICAL,        INTENT(IN) :: PCO_fr_NMVOC_use! Use P(CO) from NMVOC [T/F]
    REAL(fp),       INTENT(IN) :: PCO_fr_NMVOC    ! P(CO) fr NMVOC [molec/cm3/s]
    REAL(fp),       INTENT(IN) :: CH4_soil_k      ! CH4 soil-loss rate [s-1]
    REAL(fp),       INTENT(IN) :: dtChem          ! Chemistry timestep [s]
    TYPE(ChmState), INTENT(IN) :: State_Chm       ! Chemistry State object
    TYPE(MetState), INTENT(IN) :: State_Met       ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    !========================================================================
    ! oxidants_ComputeRateConstants begins here!
    !========================================================================

    ! Initialize
    K_OXIDANTS = 0.0_dp

    ! Prescribed oxidants are applied at every model level.  Preserve the
    ! supplied OH, Cl, and O(1D) concentrations without regional masking.
    C(ind_FixedOH)  = ConcOHmnd
    C(ind_FixedCl)  = ConcClMnd
    C(ind_FixedO1D) = ConcO1Dmnd

    IF ( PCO_fr_CH4_use ) THEN
       K_OXIDANTS(1) = PCO_fr_CH4 * OHdiurnalFac
    ENDIF

    K_OXIDANTS(2) = GC_OHCO()

    IF ( PCO_fr_NMVOC_use ) THEN
       K_OXIDANTS(3) = PCO_fr_NMVOC * OHdiurnalFac
    ENDIF

    ! Surface-only first-order CH4 uptake by soils
    K_OXIDANTS(4) = CH4_soil_k

  END SUBROUTINE oxidants_ComputeRateConstants
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: oxidants_ConvertMolecCm3ToKg
!
! !DESCRIPTION: Converts species concentrations (after chemistry) from kg to
!  molec/cm3 and stores into the State_Chm%Species concentration array.
!\\
! !INTERFACE:
!
  SUBROUTINE oxidants_ConvertMolecCm3ToKg( I, J, L, State_Chm, State_Met )
!
! !USES:
!
    USE Species_Mod,   ONLY : SpcConc
    USE State_Chm_Mod, ONLY : ChmState
    USE State_Met_Mod, ONLY : MetState
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)    :: I, J, L      ! Grid box indices
    TYPE(MetState), INTENT(IN)    :: State_Met    ! Meterorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm    ! Chemistry State object
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                :: N
    REAL(fp)               :: airvol_cm3, convfac_CO

    ! Pointers
    TYPE(SpcConc), POINTER :: Spc(:)

    !========================================================================
    ! oxidants__ConvertMolecCm3ToKg begins here!
    !========================================================================

    ! Point to species array
    Spc => State_Chm%Species

    ! Grid box volume [cm3]
    airvol_cm3  = State_Met%AIRVOL(I,J,L) * 1.0e+6_fp

    ! Send species from molec/cm3 to kg
    IF ( id_CH4 > 0 ) THEN
       Spc(id_CH4)%Conc(I,J,L) = C(ind_CH4) * airvol_cm3 / xnumol_CH4
    ENDIF

    IF ( id_CO > 0 ) THEN
       convfac_CO             = airvol_cm3 / xnumol_CO
       Spc(id_CO)%Conc(I,J,L) = C(ind_CO)  * convfac_CO
    ENDIF

    DO N = 1, N_TAGGED_CH4
       IF ( id_CH4_tagged(N) > 0 ) THEN
          Spc(id_CH4_tagged(N))%Conc(I,J,L) =                               &
               C(oxidants_Get_TaggedCH4KppId(N)) * airvol_cm3 / xnumol_CH4
       ENDIF
    ENDDO

#ifdef JACOBIAN
    ! Do the same for Jacobian CH4 tracers, if any
    IF ( numJacobianTracers > 0 ) THEN
       DO N = 1, numJacobianTracers
          Spc(JacobianIDs(N))%Conc(I,J,L) = C(JacobianIDs(N)) * airvol_cm3 / xnumol_CH4
       ENDDO
    ENDIF
#endif

    ! Free pointer
    Spc => NULL()

  END SUBROUTINE oxidants_ConvertMolecCm3ToKg
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: oxidants_get_co_ch4_flux
!
! !DESCRIPTION: Returns the flux of CO_CH4 in molec/cm3/s for diagnostics.
!\\
!\\
! !INTERFACE:
!
  FUNCTION oxidants_Get_COfromCH4_Flux( dtChem ) RESULT ( flux )
!
! !INPUT PARAMETERS:
!
    REAL(dp), INTENT(IN) :: dtChem       ! Chemistry timestep [s]
!
! !RETURN VALUE:
!
    REAL(dp)             :: flux         ! CO_CH4 flux [molec/cm3/s]
!EOP
!------------------------------------------------------------------------------
!BOC

    flux = C(ind_PCOfromCH4) / dtChem     ! molec/cm3 --> molec/cm3/s

  END FUNCTION oxidants_Get_COfromCH4_Flux
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: oxidants_get_cofromnmvoc_flux
!
! !DESCRIPTION: Returns the flux of CO_NMVOC in molec/cm3/s for diagnostics.
!\\
!\\
! !INTERFACE:
!
  FUNCTION oxidants_Get_COfromNMVOC_Flux( dtChem ) RESULT ( flux )
!
! !INPUT PARAMETERS:
!
    REAL(dp), INTENT(IN) :: dtChem       ! Chemistry timestep [s]
!
! !RETURN VALUE:
!
    REAL(dp)             :: flux         ! CO_NMVOC flux [molec/cm3/s]

!EOP
!------------------------------------------------------------------------------
!BOC

    flux = C(ind_PCOfromNMVOC) / dtChem   ! molec/cm3 --> molec/cm3/s

  END FUNCTION oxidants_Get_COfromNMVOC_Flux
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: oxidants_get_ch4byoh_loss
!
! !DESCRIPTION: Returns the grid-box CH4 loss rate from reaction with OH.
!\
!\
! !INTERFACE:
!
  FUNCTION oxidants_Get_CH4byOH_Loss( dtChem, airvol_m3 ) RESULT ( loss )
!
! !INPUT PARAMETERS:
!
    REAL(dp), INTENT(IN) :: dtChem       ! Chemistry timestep [s]
    REAL(dp), INTENT(IN) :: airvol_m3    ! Grid-box volume [m3]
!
! !RETURN VALUE:
!
    REAL(dp)             :: loss         ! CH4 loss [kg/s]
!EOP
!------------------------------------------------------------------------------
!BOC

    loss = C(ind_LCH4byOH) * airvol_m3 * 1.0e+6_dp /                        &
           ( dtChem * xnumol_CH4 )

  END FUNCTION oxidants_Get_CH4byOH_Loss
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: oxidants_get_ch4bycl_loss
!
! !DESCRIPTION: Returns the grid-box CH4 loss rate from reaction with Cl.
!\
!\
! !INTERFACE:
!
  FUNCTION oxidants_Get_CH4byCl_Loss( dtChem, airvol_m3 ) RESULT ( loss )
!
! !INPUT PARAMETERS:
!
    REAL(dp), INTENT(IN) :: dtChem       ! Chemistry timestep [s]
    REAL(dp), INTENT(IN) :: airvol_m3    ! Grid-box volume [m3]
!
! !RETURN VALUE:
!
    REAL(dp)             :: loss         ! CH4 loss [kg/s]
!EOP
!------------------------------------------------------------------------------
!BOC

    loss = C(ind_LCH4byCl) * airvol_m3 * 1.0e+6_dp /                        &
           ( dtChem * xnumol_CH4 )

  END FUNCTION oxidants_Get_CH4byCl_Loss
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: oxidants_get_ch4byo1d_loss
!
! !DESCRIPTION: Returns the grid-box CH4 loss rate from reaction with O(1D).
!\
!\
! !INTERFACE:
!
  FUNCTION oxidants_Get_CH4byO1D_Loss( dtChem, airvol_m3 ) RESULT ( loss )
!
! !INPUT PARAMETERS:
!
    REAL(dp), INTENT(IN) :: dtChem       ! Chemistry timestep [s]
    REAL(dp), INTENT(IN) :: airvol_m3    ! Grid-box volume [m3]
!
! !RETURN VALUE:
!
    REAL(dp)             :: loss         ! CH4 loss [kg/s]
!EOP
!------------------------------------------------------------------------------
!BOC

    loss = C(ind_LCH4byO1D) * airvol_m3 * 1.0e+6_dp /                       &
           ( dtChem * xnumol_CH4 )

  END FUNCTION oxidants_Get_CH4byO1D_Loss
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: oxidants_get_ch4bysoil_loss
!
! !DESCRIPTION: Returns the grid-box CH4 loss rate from soil uptake.
!\
!\
! !INTERFACE:
!
  FUNCTION oxidants_Get_CH4bySoil_Loss( dtChem, airvol_m3 ) RESULT ( loss )
!
! !INPUT PARAMETERS:
!
    REAL(dp), INTENT(IN) :: dtChem       ! Chemistry timestep [s]
    REAL(dp), INTENT(IN) :: airvol_m3    ! Grid-box volume [m3]
!
! !RETURN VALUE:
!
    REAL(dp)             :: loss         ! CH4 loss [kg/s]
!EOP
!------------------------------------------------------------------------------
!BOC

    loss = C(ind_LCH4bySoil) * airvol_m3 * 1.0e+6_dp /                      &
           ( dtChem * xnumol_CH4 )

  END FUNCTION oxidants_Get_CH4bySoil_Loss
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: oxidants_get_taggedch4_loss
!
! !DESCRIPTION: Returns a tagged CH4 loss rate for one chemical
! pathway.  Pathways 1--4 denote OH, Cl, O(1D), and soil uptake.
!\
!\
! !INTERFACE:
!
  FUNCTION oxidants_Get_TaggedCH4Loss( tag, pathway, dtChem, airvol_m3 )     &
       RESULT ( loss )
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(IN) :: tag          ! Compiled tag index
    INTEGER,  INTENT(IN) :: pathway      ! Chemical loss pathway
    REAL(dp), INTENT(IN) :: dtChem       ! Chemistry timestep [s]
    REAL(dp), INTENT(IN) :: airvol_m3    ! Grid-box volume [m3]
!
! !RETURN VALUE:
!
    REAL(dp)             :: loss         ! Tagged CH4 loss [kg/s]
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: KppId

    KppId = 0
    IF ( tag >= 1 .and. tag <= N_BASE_TAGGED_CH4 ) THEN
       SELECT CASE ( pathway )
          CASE ( 1 )
             KppId = ind_LCH4_tagged_byOH(tag)
          CASE ( 2 )
             KppId = ind_LCH4_tagged_byCl(tag)
          CASE ( 3 )
             KppId = ind_LCH4_tagged_byO1D(tag)
          CASE ( 4 )
             KppId = ind_LCH4_tagged_bySoil(tag)
       END SELECT
    ELSE IF ( tag > N_BASE_TAGGED_CH4 .and. tag <= N_TAGGED_CH4 ) THEN
       KppId = oxidants_Get_TaggedCH4KppId(tag) + pathway
    ENDIF

    loss = 0.0_dp
    IF ( KppId > 0 ) THEN
       loss = C(KppId) * airvol_m3 * 1.0e+6_dp /                            &
              ( dtChem * xnumol_CH4 )
    ENDIF

  END FUNCTION oxidants_Get_TaggedCH4Loss
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: oxidants_get_ch4soil_uptake_rate
!
! !DESCRIPTION: Converts a prescribed positive CH4 soil-uptake flux and the
!  instantaneous surface-grid-box CH4 mass into a first-order rate constant.
!  Keeping this parameterization in one function provides a clean interface
!  for a future process-based soil-uptake model.
!\
!\
! !INTERFACE:
!
  FUNCTION oxidants_Get_CH4SoilUptakeRate( uptakeFlux, gridArea, ch4Mass )     &
       RESULT ( rate )
!
! !INPUT PARAMETERS:
!
    REAL(fp), INTENT(IN) :: uptakeFlux ! Positive uptake flux [kg/m2/s]
    REAL(fp), INTENT(IN) :: gridArea   ! Grid-cell area [m2]
    REAL(fp), INTENT(IN) :: ch4Mass    ! Surface grid-box CH4 mass [kg]
!
! !RETURN VALUE:
!
    REAL(fp)             :: rate       ! First-order soil loss rate [s-1]
!EOP
!------------------------------------------------------------------------------
!BOC

    rate = 0.0_fp
    IF ( uptakeFlux <= 0.0_fp .or. gridArea <= 0.0_fp .or.                 &
         ch4Mass    <= 0.0_fp ) RETURN

    rate = SafeDiv( uptakeFlux * gridArea, ch4Mass, 0.0_fp )

  END FUNCTION oxidants_Get_CH4SoilUptakeRate
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: oxidants_get_cobyoh_loss
!
! !DESCRIPTION: Returns the grid-box CO loss rate from reaction with OH.
!\
!\
! !INTERFACE:
!
  FUNCTION oxidants_Get_CObyOH_Loss( dtChem, airvol_m3 ) RESULT ( loss )
!
! !INPUT PARAMETERS:
!
    REAL(dp), INTENT(IN) :: dtChem       ! Chemistry timestep [s]
    REAL(dp), INTENT(IN) :: airvol_m3    ! Grid-box volume [m3]
!
! !RETURN VALUE:
!
    REAL(dp)             :: loss         ! CO loss [kg/s]
!EOP
!------------------------------------------------------------------------------
!BOC

    loss = C(ind_LCObyOH) * airvol_m3 * 1.0e+6_dp /                         &
           ( dtChem * xnumol_CO )

  END FUNCTION oxidants_Get_CObyOH_Loss
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GC_OHCO
!
! !DESCRIPTION: Returns the rate of CO loss due to chemical reaction with OH
!\\
!\\
! !INTERFACE:
!
  FUNCTION GC_OHCO() RESULT( k )

!
! !REMARKS:
! Loss of CO due to chemical reaction with OH at any model level
!
!  DECAY RATE
!  The decay rate (KRATE) is calculated by:
!
!     OH + CO -> products (JPL 15-10)
!     k = (1 + 0.6Patm) * 1.5E-13
!
!  KRATE has units of [ molec^2 CO / cm6 / s ]^-1,
!  since this is a 2-body reaction.
!
! From JPL 2006: "The  reaction between HO and CO to yield
! H + CO2 akes place on a potential energy surface that
! contains the radical HOCO.  The yield of H and CO2 is
! diminished as the pressure rises.  The loss of reactants
! is thus the sum of two processes, an association to yield
! HOCO and the chemical activation process yielding H and
! CO2." So we now need two complicated reactions.
!
! GY( A0 = 5.9e-33, B0 = 1.,     A1 = 1.1e-12, B1 = -1.3e0,
!     A2 = 1.5e-13, B2 = 0.,     A3 = 2.1e09,  B3 = -6.1e0 )
!
!EOP
!------------------------------------------------------------------------------
!BOC
    REAL(dp) :: klo1,   klo2,   khi1,  khi2
    REAL(dp) :: xyrat1, xyrat2, blog1, blog2,   fexp1
    REAL(dp) :: fexp2,  kco1,   kco2,  TEMP300, k
    !
    klo1   = 5.9E-33_dp * K300_OVER_TEMP
    khi1   = 1.1E-12_dp * K300_OVER_TEMP**(-1.3_dp)
    xyrat1 = klo1 * NUMDEN / khi1
    blog1  = LOG10( xyrat1 )
    fexp1  = 1.0_dp / ( 1.0_dp + blog1*blog1 )
    kco1   = klo1 * NUMDEN * 0.6_dp**fexp1 / ( 1.0_dp + xyrat1 )

    klo2   = 1.5E-13_dp
    khi2   = 2.1E+09_dp * K300_OVER_TEMP**(-6.1_dp)
    xyrat2 = klo2 * NUMDEN / khi2
    blog2  = LOG10( xyrat2 )
    fexp2  = 1.0_dp / ( 1.0_dp + blog2*blog2 )
    kco2   = klo2 * 0.6_dp**fexp2 / ( 1.0_dp + xyrat2 )
    k      = kco1 + kco2

  END FUNCTION GC_OHCO
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: oxidants_InitKPPFuncs
!
! !DESCRIPTION: Stores species indices and kg/molecule in module variables
!  for fast lookup.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE oxidants_InitKPPFuncs( kgmolec_CH4, kgmolec_CO, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE State_Chm_Mod,  ONLY : Ind_
!
! !INPUT PARAMETERS:
!
    REAL(fp), INTENT(IN)  :: kgmolec_CH4   ! kg CH4 / molec CH4
    REAL(fp), INTENT(IN)  :: kgmolec_CO    ! kg CO  / molec CO
!
! !OUTPUT PARAMETERS:
!
    INTEGER,  INTENT(OUT) :: RC            ! Success or failure?
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: N

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc, spcName
    CHARACTER(LEN=4)   :: N_char

    !=================================================================
    ! oxidants_InitKPPFuncs begins here!
    !=================================================================

    ! Initialize
    RC       = GC_SUCCESS
    ErrMsg   = ''
    ThisLoc  = &
       ' -> at oxidants_InitKPPFuncs (in KPP/oxidants/oxidants_InitKPPFuncs.F90'

#ifdef JACOBIAN
    ! Get number of Jacobian tracers, if any
    numJacobianTracers = 0
    DO N = 1, 999 ! Set max to very high number < 1000
       write( N_char, '(I4.4)' ) N
       spcName = 'CH4_jac' // TRIM(N_char)
       IF ( Ind_(TRIM(spcName)) > 0 ) THEN
          numJacobianTracers = N
       ELSE
          EXIT
       ENDIF
    ENDDO

    ! Allocate Jacobian tracer mapping array and assign ids
    IF ( numJacobianTracers > 0 ) THEN
       ALLOCATE( JacobianIDs(numJacobianTracers), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating Jacobian tracer ID mapping array!'
          CALL GC_Error( errMsg, RC, thisLoc )
       ENDIF
       DO N = 1, numJacobianTracers
          write( N_char, '(I4.4)' ) N
          spcName = 'CH4_jac' // TRIM(N_char)
          JacobianIDs(N) = Ind_(TRIM(spcName))
       ENDDO
    ENDIF
#endif

    ! Define flags for species ID's
    id_CH4      = Ind_('CH4')
    id_CO       = Ind_('CO' )
    DO N = 1, N_TAGGED_CH4
       id_CH4_tagged(N) = Ind_( TRIM( TAGGED_CH4_NAMES(N) ), 'A' )
    ENDDO

    ! Set kg species / molec species
    xnumol_CH4  = kgmolec_CH4
    xnumol_CO   = kgmolec_CO

  END SUBROUTINE oxidants_InitKPPFuncs
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: oxidants_CleanupKPPFuncs
!
! !DESCRIPTION: Deallocates module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE oxidants_CleanupKPPFuncs( RC )
!
! !USES:
!
    USE ErrCode_Mod
!
! !OUTPUT PARAMETERS:
!
    INTEGER,  INTENT(OUT) :: RC            ! Success or failure?
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    !=================================================================
    ! oxidants_CleanupKPPFuncs begins here!
    !=================================================================

    ! Initialize
    RC       = GC_SUCCESS

#ifdef JACOBIAN
    IF ( ALLOCATED( JacobianIDs ) ) THEN
       DEALLOCATE( JacobianIDs, STAT=RC )
       CALL GC_CheckVar( 'oxidants_Funcs.F90:JacobianIDs', 1, RC )
    ENDIF
#endif

  END SUBROUTINE oxidants_CleanupKPPFuncs
!EOC
END MODULE oxidants_Funcs
