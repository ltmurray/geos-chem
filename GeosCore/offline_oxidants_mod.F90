!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: offline_oxidants_mod.F90
!
! !DESCRIPTION: Module OFFLINE_OXIDANTS_MOD contains variables and routines
! for simulating gases whose chemistry is driven by prescribed oxidants.
!\\
!\\
! !INTERFACE:
!
MODULE Offline_Oxidants_Mod
!
! !USES:
!
  USE Error_Mod,     ONLY : Safe_Div
  USE Hco_Error_Mod, ONLY : HCO_SUCCESS, HCO_FAIL, HCO_WARNING, hp
  USE PhysConstants
  USE Precision_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Run_Offline_Oxidants
  PUBLIC :: Init_Offline_Oxidants
  PUBLIC :: Cleanup_Offline_Oxidants
!
! !REVISION HISTORY:
!  04 Apr 2022 - M.S. Long   - Initial version, based on work by B. Bukosa
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!

  ! Scalars
  LOGICAL               :: useGlobOHbmk10yr
  LOGICAL               :: useGlobOHv5
  INTEGER               :: id_CH4,     id_CH4_adv, id_CO,      id_CO_adv
  INTEGER               :: id_C12H4,   id_C13H4,   id_C12H3D,  id_C13H3D
  INTEGER               :: id_OH
  REAL(fp)              :: xnumol_CH4,   xnumol_C12H4,  xnumol_C13H4
  REAL(fp)              :: xnumol_C12H3D, xnumol_C13H3D
  REAL(fp)              :: xnumol_CO,     xnumol_OH
  REAL(fp)              :: mw_CH4,        mw_C12H4,      mw_C13H4
  REAL(fp)              :: mw_C12H3D,     mw_C13H3D

  ! Arrays
  REAL(fp), ALLOCATABLE :: sumOfCosSza(:,:)
!
! !DEFINED PARAMETERS:
!
  INTEGER,  PARAMETER   :: N_CH4_DIAGS = 16
  REAL(fp), PARAMETER   :: CM2perM2    = 1.0e+4_fp
  REAL(fp), PARAMETER   :: CM3perM3    = 1.0e+6_fp
  REAL(fp), PARAMETER   :: toMolecCm3  = ( AVO / AIRMW ) * 1.0e-3_fp
  ! Defaults overridden by extension 115 options in HEMCO_Config.rc
  REAL(fp)              :: Rs_13C_VPDB = 0.0111802_fp
  REAL(fp)              :: Rs_D_VSMOW  = 1.5576e-4_fp

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chem_oxidants_gases
!
! !DESCRIPTION: Computes the chemical loss of carbon species (sources - sinks)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Run_Offline_Oxidants( Input_Opt,  State_Met,  State_Chm,            &
                                State_Grid, State_Diag, RC                   )
!
! !USES:
!
    USE oxidants_Funcs
    USE gckpp_Global
    USE gckpp_Integrator,     ONLY : Integrate
    USE gckpp_Parameters
    USE gckpp_Precision
    USE gckpp_Rates,          ONLY : Update_Rconst
    USE ErrCode_Mod
    USE HCO_State_Mod,        ONLY : Hco_GetHcoId
    USE HCO_Utilities_GC_Mod, ONLY : HCO_GC_HcoStateOK
    USE Input_Opt_Mod,        ONLY : OptInput
    USE Species_Mod,          ONLY : SpcConc
    USE State_Grid_Mod,       ONLY : GrdState
    USE State_Chm_Mod,        ONLY : ChmState
    USE State_Diag_Mod,       ONLY : DgnState
    USE State_Met_Mod,        ONLY : MetState
    USE Time_Mod,             ONLY : Get_Ts_Chem
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
    !
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! SAVEd scalars
    LOGICAL, SAVE          :: first = .TRUE.

    ! Scalars
    LOGICAL                :: failed
    INTEGER                :: HcoID,     I
    INTEGER                :: J,         L
    INTEGER                :: NA,        N
    INTEGER                :: IERR,      KppId
    INTEGER                :: S
    REAL(fp)               :: dtChem,    facDiurnal
    REAL(fp)               :: soilLossRate
    REAL(fp)               :: atomsD, atomsH, ratio
    REAL(fp)               :: rocCH4, rocD13, rocDD
    REAL(fp)               :: tsPerDay
    REAL(f4)               :: timeAfter, timeBefore

    ! Strings
    CHARACTER(LEN=255)     :: errMsg
    CHARACTER(LEN=255)     :: thisLoc

    ! Pointers
    TYPE(SpcConc), POINTER :: Spc(:)

    ! Arrays
    INTEGER                :: ICNTRL(20)
    INTEGER                :: ISTATUS(20)
    REAL(dp)               :: RCNTRL(20)
    REAL(dp)               :: RSTATE(20)
    REAL(fp)               :: OHdiurnalFac(State_Grid%NX, State_Grid%NY)

    ! Arrays for data read in via HEMCO
    REAL(fp) :: Global_OH(     State_Grid%NX, State_Grid%NY, State_Grid%NZ)
    REAL(fp) :: Global_Cl(     State_Grid%NX, State_Grid%NY, State_Grid%NZ)
    REAL(fp) :: Global_O1D(    State_Grid%NX, State_Grid%NY, State_Grid%NZ)
    REAL(fp) :: PCO_fr_CH4  (  State_Grid%NX, State_Grid%NY, State_Grid%NZ)
    REAL(fp) :: PCO_fr_NMVOC(  State_Grid%NX, State_Grid%NY, State_Grid%NZ)
    REAL(fp) :: CH4_soil_uptake(State_Grid%NX, State_Grid%NY)

    !========================================================================
    ! Run_Offline_Oxidants begins here!
    !========================================================================

    ! Initialize
    RC       =  GC_SUCCESS
    dtChem   =  Get_Ts_Chem()
    tsPerDay =  86400.0_fp / dtChem
    Spc      => State_Chm%Species
    errMsg   = ''
    thisLoc  = &
     ' -> at Run_Offline_Oxidants (in module GeosCore/offline_oxidants_mod.F90)'

    !========================================================================
    ! First-time only safety checks
    !========================================================================
    IF ( first ) THEN

       ! Exit if the HEMCO state object has not yet been initialized
       IF ( .not. HCO_GC_HcoStateOK() ) THEN
          errMsg = 'HcoState object is not associated!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       ! Determine which OH oxidant field we are using
       ! OH is only needed when CH4 or CO are included in the carbon species
       IF ( id_CH4_adv > 0 .or. id_CO_adv > 0 ) THEN
          CALL InquireGlobalOHversion( Input_Opt, RC )
          IF ( RC /= GC_SUCCESS ) THEN
             errMsg = 'Error encountered in "InquireGlobalOHversion"!'
             CALL GC_Error( errMsg, RC, thisLoc )
             RETURN
          ENDIF
       ENDIF

       ! Read the isotope-reference ratios used by extension 115
       IF ( id_C12H4 > 0 .or. id_C13H4 > 0 .or.                         &
            id_C12H3D > 0 .or. id_C13H3D > 0 ) THEN
          CALL InquireIsotopeStandards( Input_Opt, RC )
          IF ( RC /= GC_SUCCESS ) THEN
             errMsg = 'Error encountered in "InquireIsotopeStandards"!'
             CALL GC_Error( errMsg, RC, thisLoc )
             RETURN
          ENDIF
       ENDIF

       ! Reset first-time flag
       first = .FALSE.
    ENDIF

    !========================================================================
    ! Zero diagnostic arrays so that we don't have leftover values
    ! from previous timesteps hanging around at the top of the chemgrid
    !========================================================================
    CALL oxidants_ZeroDiagArrays( State_Diag )

    !========================================================================
    ! Read chemical inputs (oxidant fields, concentrations) via HEMCO
    !========================================================================
    CALL ReadChemInputFields(                                                &
         Input_Opt     = Input_Opt,                                          &
         State_Grid    = State_Grid,                                         &
         State_Met     = State_Met,                                          &
         Global_OH     = Global_OH,                                          &
         Global_Cl     = Global_Cl,                                          &
         Global_O1D    = Global_O1D,                                         &
         PCO_fr_CH4    = PCO_fr_CH4,                                         &
         PCO_fr_NMVOC  = PCO_fr_NMVOC,                                       &
         CH4_soil_uptake = CH4_soil_uptake,                                 &
         RC            = RC                                                  )

    !========================================================================
    ! Compute OH diurnal cycle scaling factor
    ! (this scales OH by the position of the sun, and zeroes it at night)
    !
    ! OH is only needed when CH4 or CO are included in the carbon species
    !========================================================================
    IF ( id_CH4_adv > 0 .or. id_CO_adv > 0 ) THEN
       CALL Calc_Diurnal(                                                    &
            State_Grid   = State_Grid,                                       &
            State_Met    = State_Met,                                        &
            OHdiurnalFac = OHdiurnalFac                                     )

       !========================================================================
       ! HISTORY: Diagnostic archival of OH [molec/cm3]
       !========================================================================
       IF ( State_Diag%Archive_OHconcAfterChem ) THEN

          !$OMP PARALLEL DO                                                  &
          !$OMP DEFAULT( SHARED                                             )&
          !$OMP PRIVATE( I, J, L                                            )&
          !$OMP COLLAPSE( 3                                                 )
          DO L = 1, State_Grid%NZ
          DO J = 1, State_Grid%NY
          DO I = 1, State_Grid%NX

             ! Archive OH if we are in the chemistry grid [molec/cm3]
             IF ( State_Met%InChemGrid(I,J,L) ) THEN
                IF ( State_Diag%Archive_OHconcAfterChem ) THEN
                   State_Diag%OHconcAfterChem(I,J,L) = Global_OH(I,J,L)         &
                                                     * OHdiurnalFac(I,J)
                ENDIF
             ENDIF

          ENDDO
          ENDDO
          ENDDO
          !$OMP END PARALLEL DO
       ENDIF
    ENDIF

    ! Do not call KPP if neither CH4 or CO are advected species
    IF ( id_CH4_adv > 0 .or. id_CO_adv > 0 ) THEN

       !=====================================================================
       ! Main chemistry loop -- call KPP to integrate the mechanism forward
       !=====================================================================

       ! KPP analytical integrator settings
       ICNTRL     =  0
       ICNTRL(15) = -1   ! Do not call Update_SUN, Update_RCONST w/in integrator

       ! Set a flag to denote if the chemistry failed
       failed     = .FALSE.

       ! Loop over grid boxes
       !$OMP PARALLEL DO                                                     &
       !$OMP DEFAULT( SHARED                                                )&
       !$OMP PRIVATE( I, J, L, N, soilLossRate, timeBefore, timeAfter       )&
       !$OMP COLLAPSE( 3                                                    )&
       !$OMP SCHEDULE( DYNAMIC, 24                                          )
       DO L = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX

          ! Initialize PRIVATE and THREADPRIVATE loop variables
          C              = 0.0_dp                    ! Species conc. [molec/cm3]
          CFACTOR        = 1.0_dp                    ! Not used, set = 1
          k_Oxidants       = 0.0_dp                    ! Carbon mechanism rates
          NUMDEN         = State_Met%AIRNUMDEN(I,J,L)! Air density [molec/cm3]
          TEMP           = State_Met%T(I,J,L)        ! Temperature [K]
          INV_TEMP       = 1.0_dp / TEMP             ! 1/T  term for equations
          TEMP_OVER_K300 = TEMP / 300.0_dp           ! T/300 term for equations
          K300_OVER_TEMP = 300.0_dp / TEMP           ! 300/T term for equations
          SUNCOS         = State_Met%SUNCOSmid(I,J)  ! Cos(SZA) ) [1]
          timeBefore     = 0.0_fp
          timeAfter      = 0.0_fp

          !==================================================================
          ! Convert species to molec/cm3 for the KPP solver
          !==================================================================
          CALL oxidants_ConvertKgtoMolecCm3(                                   &
               I          = I,                                               &
               J          = J,                                               &
               L          = L,                                               &
               State_Met  = State_Met,                                       &
               State_Chm  = State_Chm                                       )

          ! Convert the prescribed positive MeMo uptake flux [kg/m2/s] to a
          ! surface-only first-order loss rate [s-1].  The conversion is
          ! isolated behind a function that can later host a process model.
          soilLossRate = 0.0_fp
          IF ( L == 1 .and. id_CH4 > 0 ) THEN
             soilLossRate = oxidants_Get_CH4SoilUptakeRate(                   &
                  uptakeFlux = CH4_soil_uptake(I,J),                        &
                  gridArea   = State_Grid%AREA_M2(I,J),                     &
                  ch4Mass    = Spc(id_CH4)%Conc(I,J,L)                      )
          ENDIF

          !==================================================================
          ! Update reaction rates
          !==================================================================

          ! Start measuring KPP-related routine timing for this grid box
          IF ( State_Diag%Archive_KppTime ) THEN
             CALL CPU_Time( timeBefore )
          ENDIF

          ! Compute the rate constants that will be used
          CALL oxidants_ComputeRateConstants(                                  &
               I                = I,                                         &
               J                = J,                                         &
               L                = L,                                         &
               dtChem           = dtChem,                                    &
               ConcClMnd        = Global_Cl(I,J,L),                          &
               ConcOHmnd        = Global_OH(I,J,L),                          &
               ConcO1Dmnd       = Global_O1D(I,J,L),                         &
               OHdiurnalFac     = OHdiurnalFac(I,J),                         &
               PCO_fr_CH4_use   = Input_Opt%LPCO_CH4,                        &
               PCO_fr_CH4       = PCO_fr_CH4(I,J,L),                         &
               PCO_fr_NMVOC_use = Input_Opt%LPCO_NMVOC,                      &
               PCO_fr_NMVOC     = PCO_fr_NMVOC(I,J,L),                       &
               CH4_soil_k       = soilLossRate,                             &
               State_Met        = State_Met,                                 &
               State_Chm        = State_Chm                                 )

          ! Update the array of rate constants for the KPP solver
          CALL Update_RCONST()

          !==================================================================
          ! Call the KPP integrator
          !==================================================================

          ! Integrate the mechanism forward in time
          CALL Integrate(                                                    &
               TIN      = 0.0_dp,                                            &
               TOUT     = dtChem,                                            &
               ICNTRL_U = ICNTRL,                                            &
               IERR_U   = IERR                                              )

          ! Trap potential errors
          IF ( IERR /= 1 ) failed = .TRUE.

          ! Start measuring KPP-related routine timing for this grid box
          IF ( State_Diag%Archive_KppTime ) THEN
             CALL CPU_Time( timeAfter )
          ENDIF

          !==================================================================
          ! HISTORY: Update diags for KPP statistics + reaction rates
          !==================================================================
          CALL oxidants_UpdateKppDiags(                                        &
               I          = I,                                               &
               J          = J,                                               &
               L          = L,                                               &
               ISTATUS    = ISTATUS,                                         &
               timeBefore = timeBefore,                                      &
               timeAfter  = timeAfter,                                       &
               State_Diag = State_Diag                                      )

          ! Aggregate KPP family diagnostics [molec/cm3/s].  The carbon
          ! mechanism defines PCO, LCO, and LCH4 as rate accumulators; they
          ! are diagnostics, not transported species or restart requirements.
          IF ( State_Diag%Archive_Loss ) THEN
             DO S = 1, State_Diag%Map_Loss%nSlots
                KppId = State_Diag%Map_Loss%slot2Id(S)
                State_Diag%Loss(I,J,L,S) = C(KppId) / dtChem
             ENDDO
          ENDIF

          IF ( State_Diag%Archive_Prod ) THEN
             DO S = 1, State_Diag%Map_Prod%nSlots
                KppId = State_Diag%Map_Prod%slot2Id(S)
                State_Diag%Prod(I,J,L,S) = C(KppId) / dtChem
             ENDDO
          ENDIF

          !==================================================================
          ! Convert species back to kg
          !==================================================================
          CALL oxidants_ConvertMolecCm3ToKg(                                   &
               I            = I,                                             &
               J            = J,                                             &
               L            = L,                                             &
               State_Chm    = State_Chm,                                     &
               State_Met    = State_Met                                     )

          !==================================================================
          ! HISTORY: Production and loss of CO
          !
          ! NOTE: Call functions in KPP/oxidants/oxidants_Funcs.F90 so
          ! that we avoid bringing in KPP species indices into this module.
          ! This avoids compile-time dependency errors.
          !==================================================================

          ! Production of CO from CH4
          IF ( State_Diag%Archive_ProdCOfromCH4 ) THEN
             State_Diag%ProdCOfromCH4(I,J,L) =                               &
                  oxidants_Get_COfromCH4_Flux( dtChem )
          ENDIF

          ! Production of CO from NMVOCs [molec/cm3/s]
          IF ( State_Diag%Archive_ProdCOfromNMVOC ) THEN
             State_Diag%ProdCOfromNMVOC(I,J,L) =                             &
                  oxidants_Get_COfromNMVOC_Flux( dtChem )
          ENDIF

          ! Loss of CH4 by OH [kg/s]
          IF ( State_Diag%Archive_LossCH4byOH ) THEN
             State_Diag%LossCH4byOH(I,J,L) =                                 &
                  oxidants_Get_CH4byOH_Loss( dtChem, State_Met%AIRVOL(I,J,L) )
          ENDIF

          ! Loss of CH4 by Cl [kg/s]
          IF ( State_Diag%Archive_LossCH4byCl ) THEN
             State_Diag%LossCH4byCl(I,J,L) =                                 &
                  oxidants_Get_CH4byCl_Loss( dtChem, State_Met%AIRVOL(I,J,L) )
          ENDIF

          ! Loss of CH4 by O(1D) [kg/s]
          IF ( State_Diag%Archive_LossCH4byO1D ) THEN
             State_Diag%LossCH4byO1D(I,J,L) =                                &
                  oxidants_Get_CH4byO1D_Loss( dtChem, State_Met%AIRVOL(I,J,L) )
          ENDIF

          ! Loss of CH4 by first-order soil uptake [kg/s]
          IF ( State_Diag%Archive_LossCH4bySoil .and. L == 1 ) THEN
             State_Diag%LossCH4bySoil(I,J) =                                &
                  oxidants_Get_CH4bySoil_Loss(                                &
                       dtChem, State_Met%AIRVOL(I,J,L) )
          ENDIF

          ! Loss of CO by OH [kg/s]
          IF ( State_Diag%Archive_LossCObyOH ) THEN
             State_Diag%LossCObyOH(I,J,L) =                                  &
                  oxidants_Get_CObyOH_Loss( dtChem, State_Met%AIRVOL(I,J,L) )
          ENDIF

       ENDDO
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO

       ! Diagnose stable-isotope composition from molecular abundances.
       ! Molecular weights come exclusively from species_database.yml.
       IF ( id_C12H4 > 0 .and. id_C13H4 > 0 .and.                       &
            id_C12H3D > 0 .and. id_C13H3D > 0 ) THEN
          !$OMP PARALLEL DO                                               &
          !$OMP DEFAULT( SHARED                                          )&
          !$OMP PRIVATE( I, J, L, atomsD, atomsH, ratio                  )&
          !$OMP PRIVATE( rocCH4, rocD13, rocDD                           )&
          !$OMP COLLAPSE( 3                                              )
          DO L = 1, State_Grid%NZ
          DO J = 1, State_Grid%NY
          DO I = 1, State_Grid%NX

             IF ( State_Diag%Archive_d13CCH4 ) THEN
                atomsH = Spc(id_C12H4)%Conc(I,J,L)  / mw_C12H4          &
                       + Spc(id_C12H3D)%Conc(I,J,L) / mw_C12H3D
                State_Diag%d13CCH4(I,J,L) = 0.0_fp
                IF ( atomsH > 0.0_fp ) THEN
                   ratio = ( Spc(id_C13H4)%Conc(I,J,L)  / mw_C13H4       &
                           + Spc(id_C13H3D)%Conc(I,J,L) / mw_C13H3D )    &
                           / atomsH
                   State_Diag%d13CCH4(I,J,L) =                           &
                        ( ratio / Rs_13C_VPDB - 1.0_fp ) * 1000.0_fp
                ENDIF
             ENDIF

             IF ( State_Diag%Archive_d2HCH4 ) THEN
                atomsD = Spc(id_C12H3D)%Conc(I,J,L) / mw_C12H3D         &
                       + Spc(id_C13H3D)%Conc(I,J,L) / mw_C13H3D
                atomsH = 4.0_fp * Spc(id_C12H4)%Conc(I,J,L) / mw_C12H4  &
                       + 4.0_fp * Spc(id_C13H4)%Conc(I,J,L) / mw_C13H4  &
                       + 3.0_fp * Spc(id_C12H3D)%Conc(I,J,L) / mw_C12H3D&
                       + 3.0_fp * Spc(id_C13H3D)%Conc(I,J,L) / mw_C13H3D
                State_Diag%d2HCH4(I,J,L) = 0.0_fp
                IF ( atomsH > 0.0_fp ) THEN
                   ratio = atomsD / atomsH
                   State_Diag%d2HCH4(I,J,L) =                            &
                        ( ratio / Rs_D_VSMOW - 1.0_fp ) * 1000.0_fp
                ENDIF
             ENDIF

             ! Print surface methane and stable-isotope composition at the
             ! 4 x 5 grid box containing Rochester, New York.  Retain the
             ! legacy ROC record layout for real-time comparison with the
             ! GCClassic 14.5 spinup.  Radiomethane is not yet implemented
             ! in the prescribed-oxidants mechanism, so its fields are
             ! reported explicitly as unavailable.
             IF ( I == 21 .and. J == 34 .and. L == 1 ) THEN
                rocCH4 = 0.0_fp
                rocD13 = 0.0_fp
                rocDD  = 0.0_fp
                IF ( State_Met%AD(I,J,L) > 0.0_fp ) THEN
                   rocCH4 = 1.0e9_fp * Spc(id_CH4)%Conc(I,J,L)           &
                          / State_Met%AD(I,J,L) * AIRMW / mw_CH4
                ENDIF
                atomsH = Spc(id_C12H4)%Conc(I,J,L)  / mw_C12H4          &
                       + Spc(id_C12H3D)%Conc(I,J,L) / mw_C12H3D
                IF ( atomsH > 0.0_fp ) THEN
                   ratio = ( Spc(id_C13H4)%Conc(I,J,L)  / mw_C13H4       &
                           + Spc(id_C13H3D)%Conc(I,J,L) / mw_C13H3D )    &
                           / atomsH
                   rocD13 = ( ratio / Rs_13C_VPDB - 1.0_fp ) * 1000.0_fp
                ENDIF
                atomsD = Spc(id_C12H3D)%Conc(I,J,L) / mw_C12H3D         &
                       + Spc(id_C13H3D)%Conc(I,J,L) / mw_C13H3D
                atomsH = 4.0_fp * Spc(id_C12H4)%Conc(I,J,L) / mw_C12H4  &
                       + 4.0_fp * Spc(id_C13H4)%Conc(I,J,L) / mw_C13H4  &
                       + 3.0_fp * Spc(id_C12H3D)%Conc(I,J,L) / mw_C12H3D&
                       + 3.0_fp * Spc(id_C13H3D)%Conc(I,J,L) / mw_C13H3D
                IF ( atomsH > 0.0_fp ) THEN
                   rocDD = ( atomsD / atomsH / Rs_D_VSMOW - 1.0_fp )    &
                         * 1000.0_fp
                ENDIF
100             FORMAT( 'ROC: CH4 :', F10.2, ' ppb | δ13-CH4: ', F7.1, &
                     ' ‰    |   δD-CH4: ', F7.1,                       &
                     ' ‰   |   pMC-CH4: N/A | D14CH4: N/A' )
                WRITE( 6, 100 ) rocCH4, rocD13, rocDD
             ENDIF

          ENDDO
          ENDDO
          ENDDO
          !$OMP END PARALLEL DO
       ENDIF

       !=====================================================================
       ! Exit with error if the integration failed
       ! (we cannot do this from within the parallel loop)
       !=====================================================================
       IF ( failed ) THEN
          errMsg = 'KPP integration failed!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

    ENDIF

    ! Free pointers for safety's sake
    Spc => NULL()

  END SUBROUTINE Run_Offline_Oxidants
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ReadChemInputFields
!
! !DESCRIPTION: Calls HCO_EvalFld to read the various oxidant and other
!  non-emissions fields needed for the carbon simulation
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ReadChemInputFields( Input_Opt,   State_Grid,   State_Met,      &
                                  Global_OH,   Global_Cl,    Global_O1D,     &
                                  PCO_fr_CH4,  PCO_fr_NMVOC,                &
                                  CH4_soil_uptake, RC                       )
!
! !USES:
!
    USE ErrCode_Mod
    USE HCO_Utilities_GC_Mod, ONLY : HCO_GC_EvalFld
    USE Input_Opt_Mod,        ONLY : OptInput
    USE State_Grid_Mod,       ONLY : GrdState
    USE State_Met_Mod,        ONLY : MetState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt        ! Input Options object
    TYPE(GrdState), INTENT(IN)  :: State_Grid       ! Grid State object
    TYPE(MetState), INTENT(IN)  :: State_Met        ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
    REAL(fp),       INTENT(OUT) :: Global_OH(                                &
                                    State_Grid%NX,                           &
                                    State_Grid%NY,                           &
                                    State_Grid%NZ)  ! OH conc
    REAL(fp),       INTENT(OUT) :: Global_Cl(                                &
                                    State_Grid%NX,                           &
                                    State_Grid%NY,                           &
                                    State_Grid%NZ)  ! Cl conc
    REAL(fp),       INTENT(OUT) :: Global_O1D(                               &
                                    State_Grid%NX,                           &
                                    State_Grid%NY,                           &
                                    State_Grid%NZ)  ! O(1D) conc
    REAL(fp),       INTENT(OUT) :: PCO_fr_CH4(                               &
                                    State_Grid%NX,                           &
                                    State_Grid%NY,                           &
                                    State_Grid%NZ)  ! P(CO) from CH4
    REAL(fp),       INTENT(OUT) :: PCO_fr_NMVOC(                             &
                                    State_Grid%NX,                           &
                                    State_Grid%NY,                           &
                                    State_Grid%NZ)  ! P(CO) from NMVOC
    REAL(fp),       INTENT(OUT) :: CH4_soil_uptake(                          &
                                    State_Grid%NX,                           &
                                    State_Grid%NY)  ! CH4 soil uptake [kg/m2/s]
    INTEGER,        INTENT(OUT) :: RC               ! Success/failure
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: found

    ! Strings
    CHARACTER(LEN=63)  :: dgnName
    CHARACTER(LEN=255) :: thisLoc
    CHARACTER(LEN=512) :: errMsg

    !========================================================================
    !  ReadChemInputFields begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    found   = .FALSE.
    dgnName = ''
    errMsg  = ''
    thisLoc = &
     ' -> at ReadInputChemFields (in module GeosCore/offline_oxidants_mod.F90)'

    ! Initialize fields to zero
    Global_OH     = 0.0_fp
    Global_Cl     = 0.0_fp
    Global_O1D    = 0.0_fp
    PCO_fr_CH4    = 0.0_fp
    PCO_fr_NMVOC  = 0.0_fp
    CH4_soil_uptake = 0.0_fp

    ! Fields only needed if CH4 is an advected species
    IF ( id_CH4_adv > 0 ) THEN

       !------------------------------------------------------------------------
       ! Positive MeMo soil-uptake magnitude [kg/m2/s].  This is converted
       ! to a surface-only first-order KPP loss in Run_Offline_Oxidants.
       !------------------------------------------------------------------------
       DgnName = 'CH4_SOIL_UPTAKE'
       CALL HCO_GC_EvalFld( Input_Opt, State_Grid, DgnName,                  &
                            CH4_soil_uptake, RC, found=found                )
       IF ( RC /= GC_SUCCESS .or. .not. found ) THEN
          errMsg = 'Cannot get pointer to HEMCO field ' // TRIM( DgnName )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !------------------------------------------------------------------------
       ! Cl concentration:
       ! Input via HEMCO ("SpeciesConc" collection) as [mol/mol dry]
       ! Convert to [molec/cm3] below
       !------------------------------------------------------------------------

       DgnName = 'GLOBAL_Cl'
       CALL HCO_GC_EvalFld( Input_Opt, State_Grid, DgnName,                  &
                            Global_Cl, RC,         found=found              )
       IF ( RC /= GC_SUCCESS .or. .not. found ) THEN
          errMsg = 'Cannot get pointer to HEMCO field ' // TRIM( DgnName )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       ! Convert orignal units [mol/mol dry air] to [molec/cm3]
       Global_Cl = ( Global_Cl * State_Met%AirDen ) * toMolecCm3

       !------------------------------------------------------------------------
       ! O(1D) concentration supplied by HEMCO as [mol/mol dry air]
       !------------------------------------------------------------------------
       DgnName = 'GLOBAL_O1D'
       CALL HCO_GC_EvalFld( Input_Opt, State_Grid, DgnName,                  &
                            Global_O1D, RC,         found=found             )
       IF ( RC /= GC_SUCCESS .or. .not. found ) THEN
          errMsg = 'Cannot get pointer to HEMCO field ' // TRIM( DgnName )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       Global_O1D = ( Global_O1D * State_Met%AirDen ) * toMolecCm3

    ENDIF

    ! Fields only needed if CH4 or CO are advected species
    IF ( id_CH4_adv > 0 .or. id_CO_adv > 0 ) THEN

       !------------------------------------------------------------------------
       ! OH concentration: from GEOS-Chem v5 or GEOS-Chem 10yr benchmark
       !------------------------------------------------------------------------

       IF ( useGlobOHv5 .or. useGlobOHbmk10yr ) THEN

          ! NOTE: Container name is GLOBAL_OH for both data sets!
          DgnName = 'GLOBAL_OH'
          CALL HCO_GC_EvalFld( Input_Opt, State_Grid, DgnName,               &
                               Global_OH, RC,         found=found           )
          IF ( RC /= GC_SUCCESS .or. .not. found ) THEN
             errMsg = 'Cannot get pointer to HEMCO field ' // TRIM( DgnName )
             CALL GC_Error( errMsg, RC, thisLoc )
             RETURN
          ENDIF
       ENDIF

       ! If we are using OH from recent a 10-year benchmark ("SpeciesConc")
       ! then convert OH [mol/mol dry air] to [molec/cm3].
       IF ( useGlobOHbmk10yr ) THEN
          Global_OH = ( Global_OH * State_Met%AirDen ) * toMolecCm3
       ENDIF

       ! If we are using Global_OH from GEOS-Chem v5 (e.g. for the IMI or
       ! methane simulations) then convert OH from [kg/m3] to [molec/cm3].
       IF ( useGlobOHv5 ) THEN
          Global_OH = Global_OH * xnumol_OH / CM3perM3
       ENDIF

    ENDIF

    ! Fields only needed if CO is an advected species
    IF ( id_CO_adv > 0 ) THEN

       !---------------------------------------------------------------------
       ! P(CO) from CH4
       ! Input via HEMCO ("ProdCOfromCH4 field") as [molec/cm3/s]
       !---------------------------------------------------------------------
       IF ( Input_Opt%LPCO_CH4 ) THEN
          DgnName = 'PCO_CH4'
          CALL HCO_GC_EvalFld( Input_Opt,  State_Grid, DgnName,              &
                               PCO_fr_CH4, RC,         found=found          )
          IF ( RC /= GC_SUCCESS .or. .not. found ) THEN
             errMsg = 'Cannot get pointer to HEMCO field ' // TRIM( DgnName )
             CALL GC_Error( errMsg, RC, thisLoc )
             RETURN
          ENDIF
       ENDIF

       !---------------------------------------------------------------------
       ! P(CO) from NMVOC
       ! Input via HEMCO ("ProdCOfromNMVOC" field) as [molec/cm3/s]
       !---------------------------------------------------------------------
       IF ( Input_Opt%LPCO_NMVOC ) THEN
          DgnName = 'PCO_NMVOC'
          CALL HCO_GC_EvalFld( Input_Opt,    State_Grid, DgnName,            &
                               PCO_fr_NMVOC, RC,         found=found        )
          IF ( RC /= GC_SUCCESS .or. .not. found ) THEN
             errMsg = 'Cannot get pointer to HEMCO field ' // TRIM( DgnName )
             CALL GC_Error( errMsg, RC, thisLoc )
             RETURN
          ENDIF
       ENDIF

    ENDIF

  END SUBROUTINE ReadChemInputFields
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_oxidants_gases
!
! !DESCRIPTION: Allocates and zeroes module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_Offline_Oxidants( Input_Opt,  State_Chm, State_Diag,             &
                                State_Grid, RC                                )
!
! !USES:
!
    USE oxidants_Funcs,   ONLY : oxidants_InitKPPFuncs
    USE gckpp_Global,   ONLY : SR_MW, HENRY_CR, HENRY_K0
    USE gckpp_Integrator, ONLY : ValidateOxidantsAnalyticalMechanism
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Chm_Mod,  ONLY : ChmState, Ind_
    USE State_Diag_Mod, ONLY : DgnState
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(ChmState), INTENT(IN)    :: State_Chm   ! Chemistry State object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!  This routine is called from GC_INIT_EXTRA (in GeosCore/input_mod.f)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: mechanismIsValid
    INTEGER            :: KppId,  N

    ! Strings
    CHARACTER(LEN=255) :: errMsg
    CHARACTER(LEN=255) :: thisLoc

    !========================================================================
    ! Initialize
    !========================================================================
    RC         = GC_SUCCESS
    xnumol_CH4 = 1.0_fp
    xnumol_C12H4 = 1.0_fp
    xnumol_C13H4 = 1.0_fp
    xnumol_C12H3D = 1.0_fp
    xnumol_C13H3D = 1.0_fp
    xnumol_CO  = 1.0_fp
    xnumol_OH  = 1.0_fp
    errMsg     = ''
    thisLoc    = &
     ' -> at Init_Offline_Oxidants (in module GeosCore/offline_oxidants_mod.F90)'

    ! Refuse startup if the generated KPP mechanism no longer matches the
    ! reaction and accumulator ordering assumed by the analytical solver.
    CALL ValidateOxidantsAnalyticalMechanism( mechanismIsValid, errMsg )
    IF ( .not. mechanismIsValid ) THEN
       errMsg = 'Oxidants analytical integrator validation failed: ' //       &
                TRIM( errMsg )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !========================================================================
    ! Define GEOS-Chem species indices
    !
    ! Also denote which species are advected to facilitate single-tracer runs
    !========================================================================
    id_CH4        = Ind_( 'CH4'        )
    id_CH4_adv    = Ind_( 'CH4',   'A' )
    id_C12H4      = Ind_( 'C12H4'      )
    id_C13H4      = Ind_( 'C13H4'      )
    id_C12H3D     = Ind_( 'C12H3D'     )
    id_C13H3D     = Ind_( 'C13H3D'     )
    id_CO         = Ind_( 'CO'         )
    id_CO_adv     = Ind_( 'CO',    'A' )
    id_OH         = Ind_( 'FixedOH'    )

    ! Define xnumol factors = molec/kg ratios
    ! These are useful in unit conversions
    IF ( id_CH4 > 0 ) THEN
       mw_CH4     = State_Chm%SpcData(id_CH4)%Info%MW_g
       xnumol_CH4 = AVO / ( mw_CH4 * 1.0e-3_fp )
    ENDIF
    IF ( id_C12H4 > 0 ) THEN
       mw_C12H4     = State_Chm%SpcData(id_C12H4)%Info%MW_g
       xnumol_C12H4 = AVO / ( mw_C12H4 * 1.0e-3_fp )
    ENDIF
    IF ( id_C13H4 > 0 ) THEN
       mw_C13H4     = State_Chm%SpcData(id_C13H4)%Info%MW_g
       xnumol_C13H4 = AVO / ( mw_C13H4 * 1.0e-3_fp )
    ENDIF
    IF ( id_C12H3D > 0 ) THEN
       mw_C12H3D     = State_Chm%SpcData(id_C12H3D)%Info%MW_g
       xnumol_C12H3D = AVO / ( mw_C12H3D * 1.0e-3_fp )
    ENDIF
    IF ( id_C13H3D > 0 ) THEN
       mw_C13H3D     = State_Chm%SpcData(id_C13H3D)%Info%MW_g
       xnumol_C13H3D = AVO / ( mw_C13H3D * 1.0e-3_fp )
    ENDIF
    IF ( id_CO  > 0 ) &
       xnumol_CO  = AVO / ( State_Chm%SpcData(id_CO )%Info%MW_g * 1.0e-3_fp )
    IF ( id_OH  > 0 ) &
       xnumol_OH  = AVO / ( State_Chm%SpcData(id_OH )%Info%MW_g * 1.0e-3_fp )

    !========================================================================
    ! Initialize module variables in oxidants_Funcs (in KPP)
    !========================================================================
    IF ( (id_CH4 > 0) .OR. (id_CO > 0) ) THEN
       CALL oxidants_InitKPPFuncs(                                        &
            xnumol_CH4,   xnumol_C12H4,  xnumol_C13H4,                   &
            xnumol_C12H3D, xnumol_C13H3D, xnumol_CO, RC                  )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Cannot initialize module variables in oxidants_InitCarbonChem'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
    ENDIF

    !========================================================================
    ! Initialize variables for COchemistry
    !========================================================================
    ALLOCATE( sumOfCosSza( State_Grid%NX, State_Grid%NY ), STAT=RC )
    CALL GC_CheckVar( 'offline_oxidants_mod.F90:sumOfCosSza', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    sumOfCosSza = 0.0_fp

  END SUBROUTINE Init_Offline_Oxidants
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_oxidants_gases
!
! !DESCRIPTION: Subroutine CLEANUP\_GLOBAL\_CH4 deallocates module arrays.
!  (bmy, 1/16/01)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_Offline_Oxidants( RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE oxidants_Funcs,   ONLY : oxidants_CleanupKPPFuncs
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT) :: RC          ! Success or failure?
!EOP
!------------------------------------------------------------------------------
!BOC
!
    !=================================================================
    ! Cleanup_Offline_Oxidants begins here!
    !=================================================================

    ! Initialize
    RC = GC_SUCCESS

    IF ( ALLOCATED( sumOfCosSza ) ) THEN
       DEALLOCATE( sumOfCosSza, STAT=RC )
       CALL GC_CheckVar( 'offline_oxidants_mod.F90:sumOfCosSza', 2, RC )
    ENDIF

    CALL oxidants_CleanupKPPFuncs( RC )

  END SUBROUTINE Cleanup_Offline_Oxidants
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: calc_diurnal
!
! !DESCRIPTION: Subroutine CALC\_DIRUNAL computes the sume of the cosine
!  of the solar zenith angle over a 24 hour day as well as the total
!  length of daylight to scale the offline OH concentrations.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Calc_Diurnal( State_Grid, State_Met, OHdiurnalFac )
!
! !USES:
!
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Met_Mod,  ONLY : MetState
    USE TIME_MOD,       ONLY : ITS_A_NEW_DAY
    USE TIME_MOD,       ONLY : GET_MINUTE,    GET_SECOND,      GET_HOUR
    USE TIME_MOD,       ONLY : GET_TS_CHEM,   GET_DAY_OF_YEAR, GET_LOCALTIME
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN)  :: State_Grid        ! Grid State object
    TYPE(MetState), INTENT(IN)  :: State_Met         ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
    REAL(fp),       INTENT(OUT) :: OHdiurnalFac(                             &
                                    State_Grid%NX,                           &
                                    State_Grid%NY)   ! OH diurnal scaling [1]
!
! !REVISION HISTORY:
!  12 Mar 2014 - J. Fisher - Initial version, Copied from OHNO3TIME in
!                            oxidants_mod and COSSZA in dao_mod
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! SAVEd scalars
    LOGICAL, SAVE      :: FIRST = .TRUE.

    ! Scalars
    INTEGER            :: I, J, N, NDYSTEP
    INTEGER            :: SECOND,  MINUTE, TS_SUN
    REAL*8             :: GMT_MID, TIMLOC, FACTOR
    REAL*8             :: R,       AHR,    DEC
    REAL*8             :: YMID_R,  SUNTMP_MID
    REAL*8             :: dtChem,  timestepsPerDay
!
! !DEFINED PARAMETERS:
!
    ! Coefficients for solar declination angle
    REAL*8,  PARAMETER :: A0 = 0.006918d0
    REAL*8,  PARAMETER :: A1 = 0.399912d0
    REAL*8,  PARAMETER :: A2 = 0.006758d0
    REAL*8,  PARAMETER :: A3 = 0.002697d0
    REAL*8,  PARAMETER :: B1 = 0.070257d0
    REAL*8,  PARAMETER :: B2 = 0.000907d0
    REAL*8,  PARAMETER :: B3 = 0.000148d0

    !=================================================================
    ! CALC_DIURNAL begins here!
    !=================================================================

    ! Only do at the start of a new day
    IF ( FIRST .or. ITS_A_NEW_DAY() ) THEN

       ! Zero array
       sumOfCosSza = 0d0

       ! Get time for central chemistry timestep
       TS_SUN = GET_TS_CHEM()                     ! Chemistry interval
       SECOND = GET_SECOND()                      ! Current seconds
       MINUTE = GET_MINUTE()                      ! Current minutes
       FACTOR = ( MINUTE * 60 + SECOND ) / TS_SUN ! Multiplying factor

       ! GMT at the midpoint of the chemistry time interval for first
       ! timestep of the day
       GMT_MID  = ( DBLE( GET_HOUR()        )        ) &
                + ( DBLE( TS_SUN * FACTOR ) / 3600d0 ) &
                + ( DBLE( TS_SUN / 2      ) / 3600d0 )

       ! Solar declination angle (low precision formula):
       ! Path length of earth's orbit traversed since Jan 1 [radians]
       R = ( 2d0 * PI / 365d0 ) * FLOAT( GET_DAY_OF_YEAR() - 1 )
       DEC = A0 - A1*COS(    R) + B1*SIN(    R) &
                - A2*COS(2d0*R) + B2*SIN(2d0*R) &
                - A3*COS(3d0*R) + B3*SIN(3d0*R)

       ! NDYSTEP is # of chemistry time steps
       NDYSTEP = INT( 24d0 * 3600d0 / GET_TS_CHEM() )

       ! Loop forward through NDYSTEP "fake" timesteps for this day
       DO N = 1, NDYSTEP

          ! Increment GMT (hours) to midpoint of next timestep
          IF ( N > 1 ) GMT_MID = GMT_MID + TS_SUN / 3600d0

          ! Loop over surface grid boxes
          !$OMP PARALLEL DO                                                  &
          !$OMP DEFAULT( SHARED                                             )&
          !$OMP PRIVATE( I, J, YMID_R, TIMLOC, AHR, SUNTMP_MID              )&
          !$OMP COLLAPSE( 2                                                 )
          DO J = 1, State_Grid%NY
          DO I = 1, State_Grid%NX

             ! Zero SUNTMP_MID
             SUNTMP_MID = 0d0

             ! Grid box latitude center [radians]
             YMID_R = State_Grid%YMid_R(I,J)

             ! Local time at box (I,J) [hours]
             TIMLOC = GET_LOCALTIME( I, J, 1, State_Grid, GMT=GMT_MID)

             ! Hour angle at box (I,J) [radians]
             AHR = ABS( TIMLOC - 12d0 ) * 15d00 * PI_180

             !===========================================================
             ! The cosine of the solar zenith angle (SZA) is given by:
             !
             !  cos(SZA) = sin(LAT)*sin(DEC) + cos(LAT)*cos(DEC)*cos(AHR)
             !
             ! where LAT = the latitude angle,
             !       DEC = the solar declination angle,
             !       AHR = the hour angle, all in radians.
             !
             ! If SUNCOS < 0, then the sun is below the horizon, and
             ! therefore does not contribute to any solar heating.
             !===========================================================

             ! Compute Cos(SZA)
             SUNTMP_MID = sin(YMID_R) * sin(DEC) + &
                          cos(YMID_R) * cos(DEC) * cos(AHR)

             ! sumOfCosSza is the sum of SUNTMP_MID at location (I,J)
             ! Do not include negative values of SUNTMP_MID
             sumOfCosSza(I,J) = sumOfCosSza(I,J) + MAX( SUNTMP_MID, 0d0 )

          ENDDO
          ENDDO
          !$OMP END PARALLEL DO
       ENDDO

       ! Reset first-time flag
       FIRST = .FALSE.
    ENDIF

    !========================================================================
    ! Calculate the OH diurnal scaling factor
    !========================================================================

    ! Chemistry timestep [s] and timesteps per day
    dtChem          = GET_TS_CHEM()
    timestepsPerDay = 86400.0_fp / dtChem

    ! Loop over surface grid boxes
    !$OMP PARALLEL DO                                                        &
    !$OMP DEFAULT( SHARED                                                   )&
    !$OMP PRIVATE( I, J                                                     )&
    !$OMP COLLAPSE( 2                                                       )
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Initialize loop variables
       OHdiurnalFac(I,J) = 0.0_fp

       ! Scaling factor for OH diurnal cycles - zero at night
       IF ( State_Met%SUNCOSmid(I,J) > 0.0_fp  .and.                         &
            sumOfCosSza(I,J)         > 0.0_fp ) THEN
          OHdiurnalFac(I,J) = State_Met%SUNCOSmid(I,J)                       &
                            / sumOfCosSza(I,J)                               &
                            * timestepsPerDay
       ENDIF
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE Calc_Diurnal
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: InquireGlobalOHversion
!
! !DESCRIPTION: Determines if we are using global OH from a GEOS-Chem 10-year
!  benchmark, or from GEOS-Chem v5 (needed for IMI & methane simulations).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE InquireGlobalOHversion( Input_Opt, RC )
!
! !USES:
!
    USE CharPak_Mod,          ONLY : To_UpperCase
    USE ErrCode_Mod
    USE Hco_Utilities_Gc_Mod, ONLY : HCO_GC_GetOption
    USE Input_Opt_Mod,        ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC   ! Success or failure
!
! !REMARKS:
!  This has to be called on the first chemistry timestep, because the HEMCO
!  Configuration file (HEMCO_Config.rc) will not have been read when
!  Init_Offline_Oxidants is called.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: isOHon, isOHv5on

    ! Strings
    CHARACTER(LEN=255) :: optVal
    CHARACTER(LEN=255) :: thisLoc
    CHARACTER(LEN=512) :: errMsg

    !========================================================================
    ! InquireGlobalOHversion begins here!
    !========================================================================

    ! Initialize
    RC       = GC_SUCCESS
    errMsg   = ''
    thisLoc  = &
     ' -> at InquireGlobalOHversion (in module GeosCore/offline_oxidants_mod.F90)'

    ! Test if any Global OH option is on
    ! NOTE: Update the option name with each major version!!!
    optVal   = HCO_GC_GetOption( "GLOBAL_OH_GC14", extNr=0 )
    isOHon   = ( To_UpperCase( TRIM( optVal ) ) == 'TRUE' )

    ! Test for GEOS-Chem v5 OH
    ! =T
    optVal   = HCO_GC_GetOption( "GLOBAL_OH_GCv5", extNr=0 )
    isOHv5On = ( To_UpperCase( TRIM( optVal ) ) == 'TRUE' )

    ! Set a global variable to determine which OH to use
    useGlobOHbmk10yr = isOHon
    useGlobOHv5      = isOHv5on

    IF ( Input_Opt%amIRoot ) THEN
       IF ( useGlobOHv5 .or. useGlobOHbmk10yr ) THEN
          IF ( useGlobOHbmk10yr ) WRITE( 6, 100 ) 'GLOBAL_OH_GC14'
          IF ( useGlobOHv5      ) WRITE( 6, 100 ) 'GLOBAL_OH_GCv5'
 100      FORMAT( 'Offline_Oxidants: Using global OH oxidant field option: ', a )
       ELSE
          WRITE( 6, 110 )
 110      FORMAT( 'Offline_Oxidants: Global OH is set to zero!' )
       ENDIF

    ENDIF

  END SUBROUTINE InquireGlobalOHversion
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: InquireIsotopeStandards
!
! !DESCRIPTION: Reads the reference isotope ratios used by the methane
!  isotopologue HEMCO extension so emission partitioning and diagnostics use
!  exactly the same standards.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE InquireIsotopeStandards( Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Hco_Utilities_Gc_Mod, ONLY : HCO_GC_GetOption
    USE Input_Opt_Mod,        ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: IOS
    CHARACTER(LEN=255) :: optVal
    CHARACTER(LEN=255) :: thisLoc
    CHARACTER(LEN=512) :: errMsg

    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = &
     ' -> at InquireIsotopeStandards (in module GeosCore/offline_oxidants_mod.F90)'

    optVal = HCO_GC_GetOption( 'R13_VPDB', extNr=115 )
    READ( optVal, *, IOSTAT=IOS ) Rs_13C_VPDB
    IF ( IOS /= 0 .or. Rs_13C_VPDB <= 0.0_fp ) THEN
       errMsg = 'Invalid R13_VPDB option for HEMCO extension 115!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    optVal = HCO_GC_GetOption( 'RD_VSMOW', extNr=115 )
    READ( optVal, *, IOSTAT=IOS ) Rs_D_VSMOW
    IF ( IOS /= 0 .or. Rs_D_VSMOW <= 0.0_fp ) THEN
       errMsg = 'Invalid RD_VSMOW option for HEMCO extension 115!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 100 ) Rs_13C_VPDB, Rs_D_VSMOW
 100   FORMAT( 'Offline_Oxidants: Isotope reference ratios: ',            &
               'R13_VPDB=', es13.6, ', RD_VSMOW=', es13.6 )
    ENDIF

  END SUBROUTINE InquireIsotopeStandards
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: oxidants_ZeroDiagArrays
!
! !DESCRIPTION: Zeroes diagnostic archival arrays to make sure that we
!  don't have any leftover values from the last timestep near the top of
!  the chemistry grid.  This was abstracted out of Run_Offline_Oxidants
!  in order to reduce clutter.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE oxidants_ZeroDiagArrays( State_Diag )
!
! !USES:
!
    USE State_Diag_Mod, ONLY : DgnState
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DgnState), INTENT(INOUT) :: State_Diag   ! Diagnostic State object
!
! !REVISION HISTORY:
!  06 May 2025 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC

    !========================================================================
    ! oxidants_ZeroDiagArrays begins here!
    !========================================================================

    IF ( State_Diag%Archive_OHconcAfterChem ) THEN
       State_Diag%OHconcAfterChem = 0.0_f4
    ENDIF

    IF ( State_Diag%Archive_Loss ) THEN
       State_Diag%Loss = 0.0_f4
    ENDIF

    IF ( State_Diag%Archive_Prod ) THEN
       State_Diag%Prod = 0.0_f4
    ENDIF

    IF ( State_Diag%Archive_ProdCOfromCH4 ) THEN
       State_Diag%ProdCOfromCH4 = 0.0_f4
    ENDIF

    IF ( State_Diag%Archive_ProdCOfromNMVOC ) THEN
       State_Diag%ProdCOfromNMVOC = 0.0_f4
    ENDIF

    IF ( State_Diag%Archive_LossCH4byOH ) THEN
       State_Diag%LossCH4byOH = 0.0_f4
    ENDIF

    IF ( State_Diag%Archive_LossCH4byCl ) THEN
       State_Diag%LossCH4byCl = 0.0_f4
    ENDIF

    IF ( State_Diag%Archive_LossCH4byO1D ) THEN
       State_Diag%LossCH4byO1D = 0.0_f4
    ENDIF

    IF ( State_Diag%Archive_LossCH4bySoil ) THEN
       State_Diag%LossCH4bySoil = 0.0_f4
    ENDIF

    IF ( State_Diag%Archive_d13CCH4 ) THEN
       State_Diag%d13CCH4 = 0.0_f4
    ENDIF

    IF ( State_Diag%Archive_d2HCH4 ) THEN
       State_Diag%d2HCH4 = 0.0_f4
    ENDIF

    IF ( State_Diag%Archive_LossCObyOH ) THEN
       State_Diag%LossCObyOH = 0.0_f4
    ENDIF

    IF ( State_Diag%Archive_RxnRate ) THEN
       State_Diag%RxnRate = 0.0_f4
    ENDIF

    IF ( State_Diag%Archive_RxnConst ) THEN
       State_Diag%RxnConst = 0.0_f4
    ENDIF


  END SUBROUTINE oxidants_ZeroDiagArrays
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: oxidants_UpdateKppDiags
!
! !DESCRIPTION: Updates History diagnostic arrays for KPP reaction rates,
!  rate constants, and solver statistics.  This was abstracted out of
!  routine Chem_Carbon_gases to reduce clutter.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE oxidants_UpdateKppDiags( I,         J,          L,                &
                                    ISTATUS,   timeBefore, timeAfter,        &
                                    State_Diag                              )
!
! !USES:
!
    USE gckpp_Global,     ONLY : C,      RCONST
    USE gckpp_Function,   ONLY : Fun
    USE gckpp_Parameters, ONLY : NREACT, NSPEC,  NVAR
    USE gckpp_Precision
    USE State_Diag_Mod,   ONLY : DgnState
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)    :: I, J, L      ! Grid box indices
    INTEGER,        INTENT(IN)    :: ISTATUS(20)  ! KPP input options
    REAL(f4),       INTENT(IN)    :: timeBefore   ! Time before rates + integ
    REAL(f4),       INTENT(IN)    :: timeAfter    ! Time after  rates + integ
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DgnState), INTENT(INOUT) :: State_Diag   ! Diagnostic State object
!
! !REVISION HISTORY:
!  06 May 2025 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:

    ! Scalars
    INTEGER  :: N, S

    ! Arrays
    REAL(dp) :: Aout(NREACT)
    REAL(dp) :: Vloc(NVAR)

    !========================================================================
    ! HISTORY: Archive KPP solver diagnostics
    !
    ! The exact carbon integrator takes one accepted analytical step per
    ! chemistry operator call.  Retain the standard KPP counters so solver
    ! behavior remains visible to the existing diagnostic infrastructure.
    !========================================================================
    IF ( State_Diag%Archive_KppDiags ) THEN

       ! # of integrator calls
       IF ( State_Diag%Archive_KppIntCounts ) THEN
          State_Diag%KppIntCounts(I,J,L) = ISTATUS(1)
       ENDIF

       ! # of times Jacobian was constructed
       IF ( State_Diag%Archive_KppJacCounts ) THEN
          State_Diag%KppJacCounts(I,J,L) = ISTATUS(2)
       ENDIF

       ! # of internal timesteps
       IF ( State_Diag%Archive_KppTotSteps ) THEN
          State_Diag%KppTotSteps(I,J,L) = ISTATUS(3)
       ENDIF

       ! # of accepted internal timesteps
       IF ( State_Diag%Archive_KppAccSteps ) THEN
          State_Diag%KppAccSteps(I,J,L) = ISTATUS(4)
       ENDIF

       ! # of rejected internal timesteps
       IF ( State_Diag%Archive_KppRejSteps ) THEN
          State_Diag%KppRejSteps(I,J,L) = ISTATUS(5)
       ENDIF

       ! # of LU-decompositions
       IF ( State_Diag%Archive_KppLuDecomps ) THEN
          State_Diag%KppLuDecomps(I,J,L) = ISTATUS(6)
       ENDIF

       ! # of forward and backwards substitutions
       IF ( State_Diag%Archive_KppSubsts ) THEN
          State_Diag%KppSubsts(I,J,L) = ISTATUS(7)
       ENDIF

       ! # of singular-matrix decompositions
       IF ( State_Diag%Archive_KppSmDecomps ) THEN
          State_Diag%KppSmDecomps(I,J,L) = ISTATUS(8)
       ENDIF

       ! # of singular-matrix decompositions
       IF ( State_Diag%Archive_KppTime ) THEN
          State_Diag%KppTime(I,J,L) = timeAfter - timeBefore
       ENDIF
    ENDIF

    !========================================================================
    ! HISTORY: Archive KPP reaction rates [molec cm-3 s-1]
    !
    ! See gckpp_Monitor.F90 for a list of chemical reactions
    !========================================================================
    IF ( State_Diag%Archive_RxnRate                                     .or. &
         State_Diag%Archive_SatDiagnRxnRate                           ) THEN

       ! Get equation rates (Aout)
       CALL Fun( V    = C(1:NVAR),                                           &
                 F    = C(NVAR+1:NSPEC),                                     &
                 RCT  = RCONST,                                              &
                 Vdot = Vloc,                                                &
                 Aout = Aout                                                )

       ! Archive the RxnRate diagnostic collection
       IF ( State_Diag%Archive_RxnRate ) THEN
          DO S = 1, State_Diag%Map_RxnRate%nSlots
             N = State_Diag%Map_RxnRate%slot2Id(S)
             State_Diag%RxnRate(I,J,L,S) = Aout(N)
          ENDDO
       ENDIF

       ! Archive the SatDiagnRxnRate diagnostic collection
       IF ( State_Diag%Archive_SatDiagnRxnRate ) THEN
          DO S = 1, State_Diag%Map_SatDiagnRxnRate%nSlots
             N = State_Diag%Map_SatDiagnRxnRate%slot2Id(S)
             State_Diag%SatDiagnRxnRate(I,J,L,S) = Aout(N)
          ENDDO
       ENDIF
    ENDIF

    !========================================================================
    ! HISTORY: Archive KPP reaction rate constants (RCONST).
    !
    ! The units vary.  They are already updated in Update_RCONST,
    ! and do not require a call of Fun(). (hplin, 3/28/23)
    !
    ! See gckpp_Monitor.F90 for a list of chemical reactions
    !========================================================================
    IF ( State_Diag%Archive_RxnConst ) THEN
       DO S = 1, State_Diag%Map_RxnConst%nSlots
          N = State_Diag%Map_RxnConst%slot2Id(S)
          State_Diag%RxnConst(I,J,L,S) = RCONST(N)
       ENDDO
    ENDIF

  END SUBROUTINE oxidants_UpdateKppDiags
!EOC
END MODULE Offline_Oxidants_Mod
