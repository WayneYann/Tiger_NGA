!-----------------------------------------------------------------------
! ======= example file  =======
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     THIS SUBROUTINE COMPUTES RATES OF PRODUCTION CDOT
!     IN [KMOLE/(M^3S)]. THE PARAMETERS W ( REACTION RATE ),
!     K ( RATE COEFFICIENT ) AND M ( THIRD BODY CONCENTRATIONS ) ARE
!     JUST WORK SPACE FOR THIS FUNCTION.
!     C CONTAINS THE CONCENTRATIONS OF NON STEADY STATE SPECIES IN
!     [KMOLE/M^3] AND IS WORKSPACE FOR THE STEADY STATE CONCENTRATIONS,
!     WHICH ARE COMPUTED IN THIS FUNCTION.
!     TEMP IS THE TEMPERATURE IN [K] AND
!     PRESSURE IS THE PRESSURE IN [PA].
!     CALLED FUNCTIONS ARE 'GETLINDRATECOEFF', 'COMPSTEADYSTATES',
!     'CTCHZERO'
!-----------------------------------------------------------------------
SUBROUTINE PRODRATES( CDOT, W, K, C, M, TEMP, PRESSURE )
  IMPLICIT NONE
  DOUBLE PRECISION CDOT(0), W(0), K(0), C(0), M(0), TEMP, PRESSURE
  
END SUBROUTINE PRODRATES

!-----------------------------------------------------------------------
!	FILLS 'MM' WITH SPECIES MOLAR MASS IN KG/KMOLE
!-----------------------------------------------------------------------
SUBROUTINE GETMOLARMASS( MM )
  IMPLICIT NONE
  DOUBLE PRECISION MM(0)

END SUBROUTINE GETMOLARMASS

!-----------------------------------------------------------------------
!	FILLS 'NAMES' WITH SPECIES IDENTIFIER/KG
!-----------------------------------------------------------------------
SUBROUTINE GETSPECIESNAMES( NAMES )
  IMPLICIT NONE
  CHARACTER *20 NAMES(0)

END SUBROUTINE GETSPECIESNAMES

!-----------------------------------------------------------------------
!	FILLS 'NSPECIES' WITH NUMBER OF SPECIES 
!-----------------------------------------------------------------------
SUBROUTINE GETNSPECIES( NSPECIES )
  IMPLICIT NONE
  INTEGER NSPECIES
  NSPECIES=0

END SUBROUTINE GETNSPECIES

!-----------------------------------------------------------------------
!	FILLS 'NREACTIONS' WITH NUMBER OF REACTIONS 
!-----------------------------------------------------------------------
SUBROUTINE GETNREACTIONS( NREACTIONS )
  IMPLICIT NONE
  INTEGER NREACTIONS
  NREACTIONS=0
END SUBROUTINE GETNREACTIONS

!-----------------------------------------------------------------------
!     RETURNS THE NUMBER OF THE NON STEADY STATE SPECIES
!-----------------------------------------------------------------------
SUBROUTINE GETNSPECS(NSPECIES_NONS)
  IMPLICIT NONE      
  INTEGER NSPECIES_NONS
  NSPECIES_NONS=0
END SUBROUTINE GETNSPECS

!-----------------------------------------------------------------------
!     THIS SUBROUTINE COMPUTES THE STEADY STATE CONCENTRATIONS FROM
!     THE CONCENTRATIONS OF COMPUTED SPECIES AND RATE COEFFICIENTS.
!     CONCENTRATIONS OF COMPUTED SPECIES MAY NOT BE ALTERED.
!-----------------------------------------------------------------------
SUBROUTINE COMPSTEADYSTATES( K, C, M, TEMP, PRESSURE )
  IMPLICIT NONE
  DOUBLE PRECISION K(0), C(0), M(0)
  DOUBLE PRECISION TEMP, PRESSURE

END SUBROUTINE COMPSTEADYSTATES

!-----------------------------------------------------------------------
!     THIS FUNCTION COMPUTES ENTHALPY 'H' AND HEAT CAPACITY 'CP' AS
!     FUNCTION OF TEMPERATURE T FOR ALL NON STEADY STATE SPECIES
!     IN UNITS [J/KG] and [J/KG K], RESPECTIVELY.
!     THE PARAMETER H AND CP SHOULD PROVIDE WORKSPACE OF LENGTH 43
!-----------------------------------------------------------------------
SUBROUTINE COMPTHERMODATA( H, CP, T )
  IMPLICIT NONE
  DOUBLE PRECISION H(0), CP(0), T

END SUBROUTINE COMPTHERMODATA

!-----------------------------------------------------------------------
!	    FILLS 'KOVEREPS' WITH KOVEREPS
!-----------------------------------------------------------------------
SUBROUTINE GETKOVEREPS( KOVEREPS )
  IMPLICIT NONE
  DOUBLE PRECISION KOVEREPS(0)
END SUBROUTINE GETKOVEREPS

!-----------------------------------------------------------------------
!	FILLS 'MUCOEFF' WITH MUECOEFF IN KG/(M*S)
!-----------------------------------------------------------------------
SUBROUTINE GETMUCOEFF( MUCOEFF )
  IMPLICIT NONE
  DOUBLE PRECISION MUCOEFF(0)

END SUBROUTINE GETMUCOEFF
