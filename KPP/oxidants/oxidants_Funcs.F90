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
  INTEGER, PARAMETER :: N_TAGGED_CH4 = 14
  INTEGER            :: id_CH4_sector(N_TAGGED_CH4)

  CHARACTER(LEN=7), PARAMETER :: TAGGED_CH4_NAMES(N_TAGGED_CH4) = (/       &
       'CH4_OIL', 'CH4_GAS', 'CH4_COL', 'CH4_LIV', 'CH4_LDF', 'CH4_WST', &
       'CH4_RIC', 'CH4_OTA', 'CH4_BBN', 'CH4_WTL', 'CH4_SEE', 'CH4_LAK', &
       'CH4_TER', 'CH4_RES' /)

  INTEGER, PARAMETER :: ind_CH4_sector(N_TAGGED_CH4) = (/                  &
       ind_CH4_OIL, ind_CH4_GAS, ind_CH4_COL, &
       ind_CH4_LIV, ind_CH4_LDF, ind_CH4_WST, &
       ind_CH4_RIC, ind_CH4_OTA, ind_CH4_BBN, &
       ind_CH4_WTL, ind_CH4_SEE, ind_CH4_LAK, &
       ind_CH4_TER, ind_CH4_RES /)
  INTEGER, PARAMETER :: ind_LCH4_sector_byOH(N_TAGGED_CH4) = (/             &
       ind_LCH4_OILbyOH, ind_LCH4_GASbyOH, ind_LCH4_COLbyOH, &
       ind_LCH4_LIVbyOH, ind_LCH4_LDFbyOH, ind_LCH4_WSTbyOH, &
       ind_LCH4_RICbyOH, ind_LCH4_OTAbyOH, ind_LCH4_BBNbyOH, &
       ind_LCH4_WTLbyOH, ind_LCH4_SEEbyOH, ind_LCH4_LAKbyOH, &
       ind_LCH4_TERbyOH, ind_LCH4_RESbyOH /)
  INTEGER, PARAMETER :: ind_LCH4_sector_byCl(N_TAGGED_CH4) = (/             &
       ind_LCH4_OILbyCl, ind_LCH4_GASbyCl, ind_LCH4_COLbyCl, &
       ind_LCH4_LIVbyCl, ind_LCH4_LDFbyCl, ind_LCH4_WSTbyCl, &
       ind_LCH4_RICbyCl, ind_LCH4_OTAbyCl, ind_LCH4_BBNbyCl, &
       ind_LCH4_WTLbyCl, ind_LCH4_SEEbyCl, ind_LCH4_LAKbyCl, &
       ind_LCH4_TERbyCl, ind_LCH4_RESbyCl /)
  INTEGER, PARAMETER :: ind_LCH4_sector_byO1D(N_TAGGED_CH4) = (/            &
       ind_LCH4_OILbyO1D, ind_LCH4_GASbyO1D, ind_LCH4_COLbyO1D, &
       ind_LCH4_LIVbyO1D, ind_LCH4_LDFbyO1D, ind_LCH4_WSTbyO1D, &
       ind_LCH4_RICbyO1D, ind_LCH4_OTAbyO1D, ind_LCH4_BBNbyO1D, &
       ind_LCH4_WTLbyO1D, ind_LCH4_SEEbyO1D, ind_LCH4_LAKbyO1D, &
       ind_LCH4_TERbyO1D, ind_LCH4_RESbyO1D /)
  INTEGER, PARAMETER :: ind_LCH4_sector_bySoil(N_TAGGED_CH4) = (/           &
       ind_LCH4_OILbySoil, ind_LCH4_GASbySoil, ind_LCH4_COLbySoil, &
       ind_LCH4_LIVbySoil, ind_LCH4_LDFbySoil, ind_LCH4_WSTbySoil, &
       ind_LCH4_RICbySoil, ind_LCH4_OTAbySoil, ind_LCH4_BBNbySoil, &
       ind_LCH4_WTLbySoil, ind_LCH4_SEEbySoil, ind_LCH4_LAKbySoil, &
       ind_LCH4_TERbySoil, ind_LCH4_RESbySoil /)

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
       IF ( id_CH4_sector(N) > 0 ) THEN
          C(ind_CH4_sector(N)) = Spc(id_CH4_sector(N))%Conc(I,J,L)           &
                                * xnumol_CH4 / airvol_cm3
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
       IF ( id_CH4_sector(N) > 0 ) THEN
          Spc(id_CH4_sector(N))%Conc(I,J,L) = C(ind_CH4_sector(N))           &
                                             * airvol_cm3 / xnumol_CH4
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
! !DESCRIPTION: Returns a sector-tagged CH4 loss rate for one chemical
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
    INTEGER,  INTENT(IN) :: tag          ! Sector tag index
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
    IF ( tag >= 1 .and. tag <= N_TAGGED_CH4 ) THEN
       SELECT CASE ( pathway )
          CASE ( 1 )
             KppId = ind_LCH4_sector_byOH(tag)
          CASE ( 2 )
             KppId = ind_LCH4_sector_byCl(tag)
          CASE ( 3 )
             KppId = ind_LCH4_sector_byO1D(tag)
          CASE ( 4 )
             KppId = ind_LCH4_sector_bySoil(tag)
       END SELECT
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
       id_CH4_sector(N) = Ind_( TRIM( TAGGED_CH4_NAMES(N) ), 'A' )
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
