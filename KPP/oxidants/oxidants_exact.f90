!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Exact analytical integrator for the linear GEOS-Chem prescribed oxidants mechanism. !
! Rates are held fixed over one chemistry operator timestep.              !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
MODULE KPP_ROOT_Integrator

  USE KPP_ROOT_Precision,  ONLY : dp
  USE KPP_ROOT_Parameters
  USE KPP_ROOT_Global
  USE KPP_ROOT_Monitor,    ONLY : SPC_NAMES, FAM_NAMES, EQN_NAMES

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Integrate
  PUBLIC :: ValidateOxidantsAnalyticalMechanism

  LOGICAL, PRIVATE :: Do_Update_RCONST
  LOGICAL, PRIVATE :: Do_Update_PHOTO
  LOGICAL, PRIVATE :: Do_Update_SUN

CONTAINS

  SUBROUTINE Integrate( TIN,       TOUT,      ICNTRL_U, RCNTRL_U,          &
                        ISTATUS_U, RSTATUS_U, IERR_U )

    USE KPP_ROOT_Rates, ONLY : Update_SUN, Update_RCONST
    USE KPP_ROOT_Util,  ONLY : Integrator_Update_Options

    REAL(dp), INTENT(IN)            :: TIN
    REAL(dp), INTENT(IN)            :: TOUT
    INTEGER,  INTENT(IN),  OPTIONAL :: ICNTRL_U(20)
    REAL(dp), INTENT(IN),  OPTIONAL :: RCNTRL_U(20)
    INTEGER,  INTENT(OUT), OPTIONAL :: ISTATUS_U(20)
    REAL(dp), INTENT(OUT), OPTIONAL :: RSTATUS_U(20)
    INTEGER,  INTENT(OUT), OPTIONAL :: IERR_U

    INTEGER  :: ICNTRL(20), ISTATUS(20), IERR
    REAL(dp) :: RSTATUS(20), dt

    ICNTRL  = 0
    ISTATUS = 0
    RSTATUS = 0.0_dp
    IERR    = 1

    ! Match the standard KPP integrator convention.  GEOS-Chem passes -1
    ! because it updates rates immediately before calling Integrate.
    ICNTRL(15) = 5
    IF ( PRESENT( ICNTRL_U ) ) THEN
       WHERE ( ICNTRL_U /= 0 ) ICNTRL = ICNTRL_U
    ENDIF

    CALL Integrator_Update_Options( ICNTRL(15),       Do_Update_RCONST,     &
                                    Do_Update_PHOTO,  Do_Update_SUN )

    VAR => C(1:NVAR)
    FIX => C(NVAR+1:NSPEC)

    IF ( Do_Update_SUN    ) CALL Update_SUN()
    IF ( Do_Update_RCONST ) CALL Update_RCONST()

    dt = TOUT - TIN
    IF ( dt < 0.0_dp ) THEN
       IERR = -1
    ELSE IF ( dt > 0.0_dp ) THEN
       CALL OxidantsAnalyticalStep( dt )
       ISTATUS(1) = 1 ! One analytical function evaluation
       ISTATUS(3) = 1 ! One total step
       ISTATUS(4) = 1 ! One accepted step
    ENDIF

    VAR => NULL()
    FIX => NULL()

    IF ( PRESENT( ISTATUS_U ) ) ISTATUS_U = ISTATUS
    IF ( PRESENT( RSTATUS_U ) ) RSTATUS_U = RSTATUS
    IF ( PRESENT( IERR_U    ) ) IERR_U    = IERR

  END SUBROUTINE Integrate

  SUBROUTINE OxidantsAnalyticalStep( dt )

    REAL(dp), INTENT(IN) :: dt

    INTEGER, PARAMETER :: nTaggedCH4 = 14
    INTEGER, PARAMETER :: indTaggedCH4(nTaggedCH4) = (/                    &
         ind_CH4_OIL, ind_CH4_GAS, ind_CH4_COL, &
         ind_CH4_LIV, ind_CH4_LDF, ind_CH4_WST, &
         ind_CH4_RIC, ind_CH4_OTA, ind_CH4_BBN, &
         ind_CH4_WTL, ind_CH4_SEE, ind_CH4_LAK, &
         ind_CH4_TER, ind_CH4_RES /)
    INTEGER, PARAMETER :: indTaggedLossOH(nTaggedCH4) = (/                  &
         ind_LCH4_OILbyOH, ind_LCH4_GASbyOH, ind_LCH4_COLbyOH, &
         ind_LCH4_LIVbyOH, ind_LCH4_LDFbyOH, ind_LCH4_WSTbyOH, &
         ind_LCH4_RICbyOH, ind_LCH4_OTAbyOH, ind_LCH4_BBNbyOH, &
         ind_LCH4_WTLbyOH, ind_LCH4_SEEbyOH, ind_LCH4_LAKbyOH, &
         ind_LCH4_TERbyOH, ind_LCH4_RESbyOH /)
    INTEGER, PARAMETER :: indTaggedLossCl(nTaggedCH4) = (/                  &
         ind_LCH4_OILbyCl, ind_LCH4_GASbyCl, ind_LCH4_COLbyCl, &
         ind_LCH4_LIVbyCl, ind_LCH4_LDFbyCl, ind_LCH4_WSTbyCl, &
         ind_LCH4_RICbyCl, ind_LCH4_OTAbyCl, ind_LCH4_BBNbyCl, &
         ind_LCH4_WTLbyCl, ind_LCH4_SEEbyCl, ind_LCH4_LAKbyCl, &
         ind_LCH4_TERbyCl, ind_LCH4_RESbyCl /)
    INTEGER, PARAMETER :: indTaggedLossO1D(nTaggedCH4) = (/                 &
         ind_LCH4_OILbyO1D, ind_LCH4_GASbyO1D, ind_LCH4_COLbyO1D, &
         ind_LCH4_LIVbyO1D, ind_LCH4_LDFbyO1D, ind_LCH4_WSTbyO1D, &
         ind_LCH4_RICbyO1D, ind_LCH4_OTAbyO1D, ind_LCH4_BBNbyO1D, &
         ind_LCH4_WTLbyO1D, ind_LCH4_SEEbyO1D, ind_LCH4_LAKbyO1D, &
         ind_LCH4_TERbyO1D, ind_LCH4_RESbyO1D /)
    INTEGER, PARAMETER :: indTaggedLossSoil(nTaggedCH4) = (/                &
         ind_LCH4_OILbySoil, ind_LCH4_GASbySoil, ind_LCH4_COLbySoil, &
         ind_LCH4_LIVbySoil, ind_LCH4_LDFbySoil, ind_LCH4_WSTbySoil, &
         ind_LCH4_RICbySoil, ind_LCH4_OTAbySoil, ind_LCH4_BBNbySoil, &
         ind_LCH4_WTLbySoil, ind_LCH4_SEEbySoil, ind_LCH4_LAKbySoil, &
         ind_LCH4_TERbySoil, ind_LCH4_RESbySoil /)

    INTEGER  :: N
    REAL(dp) :: ch4Initial, ch4Loss, ch4Remaining
    REAL(dp) :: coInitial, coLoss, coRemaining
    REAL(dp) :: kCH4OH, kCH4Cl, kCH4O1D, kCH4Soil, kCH4
    REAL(dp) :: kCO, pCOCH4, pCONMVOC, pCO, producedCO
    REAL(dp) :: lossFraction, producedLossFraction

    !---------------------------------------------------------------------
    ! CH4: exact decay and exact allocation to the four loss channels.
    !---------------------------------------------------------------------
    ch4Initial = MAX( C(ind_CH4), 0.0_dp )
    kCH4OH     = MAX( RCONST(1) * C(ind_FixedOH),  0.0_dp )
    kCH4Cl     = MAX( RCONST(2) * C(ind_FixedCl),  0.0_dp )
    kCH4O1D    = MAX( RCONST(3) * C(ind_FixedO1D), 0.0_dp )
    kCH4Soil   = MAX( RCONST(4),                    0.0_dp )
    kCH4       = kCH4OH + kCH4Cl + kCH4O1D + kCH4Soil

    IF ( kCH4 > 0.0_dp .and. ch4Initial > 0.0_dp ) THEN
       lossFraction = OneMinusExpNeg( kCH4 * dt )
       ch4Loss      = ch4Initial * lossFraction
       ch4Remaining = ch4Initial - ch4Loss

       C(ind_LCH4byOH)   = C(ind_LCH4byOH)                         &
                          + ch4Loss * ( kCH4OH   / kCH4 )
       C(ind_LCH4byCl)   = C(ind_LCH4byCl)                         &
                          + ch4Loss * ( kCH4Cl   / kCH4 )
       C(ind_LCH4byO1D)  = C(ind_LCH4byO1D)                        &
                          + ch4Loss * ( kCH4O1D  / kCH4 )
       C(ind_LCH4bySoil) = C(ind_LCH4bySoil)                       &
                          + ch4Loss * ( kCH4Soil / kCH4 )
       C(ind_LCH4)       = C(ind_LCH4) + ch4Loss
       C(ind_CH4)        = MAX( ch4Remaining, 0.0_dp )
    ELSE
       C(ind_CH4) = ch4Initial
    ENDIF

    !---------------------------------------------------------------------
    ! Sector-tagged CH4: exact decay with the same four rate constants.
    ! Tagged carbon terminates in fixed DummyCO until tagged CO is coupled.
    !---------------------------------------------------------------------
    DO N = 1, nTaggedCH4
       ch4Initial = MAX( C(indTaggedCH4(N)), 0.0_dp )
       IF ( kCH4 > 0.0_dp .and. ch4Initial > 0.0_dp ) THEN
          lossFraction = OneMinusExpNeg( kCH4 * dt )
          ch4Loss      = ch4Initial * lossFraction
          ch4Remaining = ch4Initial - ch4Loss

          C(indTaggedLossOH(N)) = C(indTaggedLossOH(N))                 &
                                + ch4Loss * ( kCH4OH / kCH4 )
          C(indTaggedLossCl(N)) = C(indTaggedLossCl(N))                 &
                                + ch4Loss * ( kCH4Cl / kCH4 )
          C(indTaggedLossO1D(N)) = C(indTaggedLossO1D(N))               &
                                 + ch4Loss * ( kCH4O1D / kCH4 )
          C(indTaggedLossSoil(N)) = C(indTaggedLossSoil(N))             &
                                  + ch4Loss * ( kCH4Soil / kCH4 )
          C(indTaggedCH4(N)) = MAX( ch4Remaining, 0.0_dp )
       ELSE
          C(indTaggedCH4(N)) = ch4Initial
       ENDIF
    ENDDO


    !---------------------------------------------------------------------
    ! CO: exact solution of dCO/dt = P - k*CO.  The loss of CO produced
    ! during the step is included so that production - loss equals dCO.
    !---------------------------------------------------------------------
    coInitial = MAX( C(ind_CO), 0.0_dp )
    pCOCH4    = MAX( RCONST(5) * C(ind_DummyCH4),   0.0_dp )
    kCO       = MAX( RCONST(6) * C(ind_FixedOH),    0.0_dp )
    pCONMVOC  = MAX( RCONST(7) * C(ind_DummyNMVOC), 0.0_dp )
    pCO       = pCOCH4 + pCONMVOC
    producedCO = pCO * dt

    IF ( kCO > 0.0_dp ) THEN
       lossFraction = OneMinusExpNeg( kCO * dt )
       producedLossFraction = ProductionLossIntegral( kCO * dt )
       coLoss = coInitial * lossFraction                                  &
                + producedCO * producedLossFraction
    ELSE
       coLoss = 0.0_dp
    ENDIF

    coRemaining = coInitial + producedCO - coLoss
    C(ind_CO)            = MAX( coRemaining, 0.0_dp )
    C(ind_PCOfromCH4)    = C(ind_PCOfromCH4) + pCOCH4 * dt
    C(ind_PCOfromNMVOC)  = C(ind_PCOfromNMVOC) + pCONMVOC * dt
    C(ind_PCO)           = C(ind_PCO) + producedCO
    C(ind_LCObyOH)       = C(ind_LCObyOH) + coLoss
    C(ind_LCO)           = C(ind_LCO) + coLoss

  END SUBROUTINE OxidantsAnalyticalStep

  PURE FUNCTION OneMinusExpNeg( x ) RESULT( value )

    REAL(dp), INTENT(IN) :: x
    REAL(dp)             :: value

    IF ( x < 1.0e-4_dp ) THEN
       value = x * ( 1.0_dp - x * ( 0.5_dp - x *                         &
               ( 1.0_dp / 6.0_dp - x / 24.0_dp ) ) )
    ELSE
       value = 1.0_dp - EXP( -x )
    ENDIF

  END FUNCTION OneMinusExpNeg

  PURE FUNCTION ProductionLossIntegral( x ) RESULT( value )

    ! Fraction of uniformly produced material that is lost during dt:
    ! 1 - (1-exp(-x))/x.  Use a series to avoid cancellation for small x.
    REAL(dp), INTENT(IN) :: x
    REAL(dp)             :: value

    IF ( x < 1.0e-4_dp ) THEN
       value = x * ( 0.5_dp - x * ( 1.0_dp / 6.0_dp - x *               &
               ( 1.0_dp / 24.0_dp - x / 120.0_dp ) ) )
    ELSE
       value = 1.0_dp - OneMinusExpNeg( x ) / x
    ENDIF

  END FUNCTION ProductionLossIntegral

  SUBROUTINE ValidateOxidantsAnalyticalMechanism( isValid, message )

    LOGICAL,          INTENT(OUT) :: isValid
    CHARACTER(LEN=*), INTENT(OUT) :: message

    isValid = .FALSE.
    message = ''

    IF ( NSPEC /= 88 .or. NVAR /= 82 .or. NVARACT /= 16 .or.             &
         NREACT /= 63 .or. NFAM /= 3 ) THEN
       message = 'Unexpected oxidants KPP dimensions or reaction count'
       RETURN
    ENDIF

    IF ( TRIM( SPC_NAMES(ind_CH4)          ) /= 'CH4'          .or.       &
         TRIM( SPC_NAMES(ind_CO)           ) /= 'CO'           .or.       &
         TRIM( SPC_NAMES(ind_PCOfromCH4)   ) /= 'PCOfromCH4'   .or.       &
         TRIM( SPC_NAMES(ind_PCOfromNMVOC) ) /= 'PCOfromNMVOC' .or.       &
         TRIM( SPC_NAMES(ind_LCH4byOH)     ) /= 'LCH4byOH'     .or.       &
         TRIM( SPC_NAMES(ind_LCH4byCl)     ) /= 'LCH4byCl'     .or.       &
         TRIM( SPC_NAMES(ind_LCH4byO1D)    ) /= 'LCH4byO1D'    .or.       &
         TRIM( SPC_NAMES(ind_LCH4bySoil)   ) /= 'LCH4bySoil'   .or.       &
         TRIM( SPC_NAMES(ind_LCObyOH)      ) /= 'LCObyOH'      .or.       &
         TRIM( SPC_NAMES(ind_PCO)          ) /= 'PCO'          .or.       &
         TRIM( SPC_NAMES(ind_LCO)          ) /= 'LCO'          .or.       &
         TRIM( SPC_NAMES(ind_LCH4)         ) /= 'LCH4'         .or.       &
         TRIM( SPC_NAMES(ind_FixedOH)      ) /= 'FixedOH'      .or.       &
         TRIM( SPC_NAMES(ind_FixedCl)      ) /= 'FixedCl'      .or.       &
         TRIM( SPC_NAMES(ind_FixedO1D)     ) /= 'FixedO1D'     .or.       &
         TRIM( SPC_NAMES(ind_DummyCH4)     ) /= 'DummyCH4'     .or.       &
         TRIM( SPC_NAMES(ind_DummyNMVOC)   ) /= 'DummyNMVOC'   .or.       &
         TRIM( SPC_NAMES(ind_DummyCO)      ) /= 'DummyCO' ) THEN
       message = 'Unexpected oxidants KPP species ordering'
       RETURN
    ENDIF

    IF ( TRIM( FAM_NAMES(1) ) /= 'PCO'  .or.                              &
         TRIM( FAM_NAMES(2) ) /= 'LCO'  .or.                              &
         TRIM( FAM_NAMES(3) ) /= 'LCH4' ) THEN
       message = 'Unexpected oxidants KPP family ordering'
       RETURN
    ENDIF

    IF ( INDEX( EQN_NAMES(1), 'CH4 + FixedOH'  ) == 0 .or.                 &
         INDEX( EQN_NAMES(1), 'LCH4byOH'       ) == 0 .or.                 &
         INDEX( EQN_NAMES(2), 'CH4 + FixedCl'  ) == 0 .or.                 &
         INDEX( EQN_NAMES(2), 'LCH4byCl'       ) == 0 .or.                 &
         INDEX( EQN_NAMES(3), 'CH4 + FixedO1D' ) == 0 .or.                 &
         INDEX( EQN_NAMES(3), 'LCH4byO1D'      ) == 0 .or.                 &
         INDEX( EQN_NAMES(4), 'CH4 --> LCH4bySoil' ) == 0 .or.             &
         INDEX( EQN_NAMES(5), 'DummyCH4 --> CO' ) == 0 .or.                &
         INDEX( EQN_NAMES(6), 'CO + FixedOH'   ) == 0 .or.                 &
         INDEX( EQN_NAMES(6), 'LCObyOH'        ) == 0 .or.                 &
         INDEX( EQN_NAMES(7),  'DummyNMVOC --> CO' ) == 0 .or.            &
         INDEX( EQN_NAMES(8),  'CH4_OIL + FixedOH' ) == 0 .or.            &
         INDEX( EQN_NAMES(8),  'LCH4_OILbyOH' ) == 0 .or.                 &
         INDEX( EQN_NAMES(63), 'CH4_RES --> LCH4_RESbySoil' ) == 0 .or.   &
         INDEX( EQN_NAMES(63), 'LCH4_RESbySoil' ) == 0 ) THEN
       message = 'Oxidants equations changed; analytical solver requires review'
       RETURN
    ENDIF

    isValid = .TRUE.
    message = 'Oxidants analytical mechanism signature validated'

  END SUBROUTINE ValidateOxidantsAnalyticalMechanism

END MODULE KPP_ROOT_Integrator
