#!/bin/bash

#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !IROUTINE: singleOxidantsSpecies.sh
#
# !DESCRIPTION: Updates oxidants simulation configuration files so that
#  a simulation with a single species can be performed.
#
# !CALLING SEQUENCE:
#  ./singleOxidantsSpecies.sh <species-to-retain> <path-to-rundir>
#
# !REVISION HISTORY:
#  14 Sep 2023 - R. Yantosca - Initial version
#  See the subsequent Git history with the gitk browser!
#EOP
#------------------------------------------------------------------------------
#BOC

# List of all advected species in the oxidants simulation
ALL_SPECIES=(CO CH4)


function isItemInList() {

    #=========================================================================
    # Test if an item is in a list.
    #
    # Arguments:
    # ${1}: The item
    # ${2}: The list
    #
    # Returns (via $?)
    # 0 if item is in the list
    # 1 if item is not in the list
    #
    # See stackoverflow.com/questions/8063228/check-if-a-variable-exists-in-a-list-in-bash
    #=========================================================================
    echo "${2}" | tr ' '  '\n' | grep -F -x -q "${1}"
}


function keyValueUpdate() {

    #=========================================================================
    # Runs an sed command to update a "key: value" pair in a file
    #
    # Arguments:
    # ${1} : Key
    # ${2} : Value
    # ${3} : Replacement value
    # ${4) : File in which the text is found
    #=========================================================================
    cmd="s/${1}: ${2}/${1}: ${3}/"
    sed -i -e "$cmd" "${4}"
}


function speciesToExclude() {

    #========================================================================
    # Returns the list of species to exclude
    #
    # Arguments:
    # ${1} : Species to retain
    #
    # Returns:
    # List of species to exclude
    #========================================================================

    # All species
    list=("${ALL_SPECIES[@]}")

    # Keep all species except for the 1st argument
    for i in "${!list[@]}"; do
        if [[ "x${list[i]}" == "x${1}" ]]; then
            unset 'list[i]'
        fi
    done

    # Print the result (that's how we return strings)
    result=""
    for spc in "${list[@]}"; do
        result+="${spc} "
    done
    echo $result
}


function updateGeosChemConfig() {

    #========================================================================
    # Removes advected species from geoschem_config.yml
    #
    # Arguments:
    # ${1} : List of species to exclude
    # ${2} : Path to the run directory
    #========================================================================

    # File to modify
    file="${2}/geoschem_config.yml"

    # Remove advected species (use \< and \> for an exact match).
    for spc in ${1}; do
    cmd="/\s*\-\s*\<${spc}\>/d"
        sed -i -e "${cmd}" "${file}"
    done
}


function updateHemcoConfig() {

    #========================================================================
    # Removes advected species from geoschem_config.yml
    #
    # Arguments:
    # ${1} : List of species to exclude
    # ${2} : Path to the run directory
    #========================================================================

    # File to be modified
    file="${2}/HEMCO_Config.rc"

    # True/False values in the HEMCO extension switches section
    true="      true"
    false="      false"

    # If CO is in the exclude list, turn off CO options
    isItemInList "CO" "${1}"
    if [[ $? == 0 ]]; then
        key="--> USE_CO_DATA            "
        keyValueUpdate "${key}" "${true}" "${false}" "${file}"
    fi

    # If CH4 is in the exclude list, turn off CH4 options
    isItemInList "CH4" "${1}"
    if [[ $? == 0 ]]; then
        key="--> USE_CH4_DATA           "
        keyValueUpdate "${key}" "${true}" "${false}" "${file}"
    fi

}


function updateHemcoDiagn() {

    #========================================================================
    # Comments out lines for unused species in HEMCO_Diagn.rc
    #
    # Arguments:
    # ${1} : List of species to exclude
    # ${2} : Path to the run directory
    #========================================================================

    # File to modify
    file="${2}/HEMCO_Diagn.rc"

    # Remove entries for excluded species
    exclude=("${1}")
    for spc in ${exclude[@]}; do
        sed -i "/Emis${spc}_/d" "${file}"
    done
    sed -i "/#####/d" "${file}"
}


function updateHistory() {

    #========================================================================
    # Removes entries in HISTORY.rc for unused species
    #
    # Arguments:
    # ${1} : List of species to exclude
    # ${2} : Path to the run directory
    # ${3} : List of species to include
    #========================================================================

    # File to be modified
    file="${2}/HISTORY.rc"

    # For GCHP: remove entries for species to be excluded
    exclude=("${1}")
    for spc in ${exclude[@]}; do
        sed -i "/\_${spc} /d"   "${file}" # trailing space required
        sed -i "/Emis${spc}_/d" "${file}"
    done

    # Restore Collection.fields line
    if [[ ! ${3} =~ "CH4" ]]; then
	oldline="                            'Emis${3}_Total"
	newline="Emissions.fields:           'Emis${3}_Total"
	sed -i "s|$oldline|$newline|g" "${file}"

	oldline="                            'SpeciesConcVV_${3}"
	newline="SpeciesConc.fields:         'SpeciesConcVV_${3}"
	sed -i "s|$oldline|$newline|g" "${file}"

	oldline="                            'BudgetEmisDryDepFull_${3}"
	newline="Budget.fields:              'BudgetEmisDryDepFull_${3}"
	sed -i "s|$oldline|$newline|g" "${file}"

	oldline="                              'CloudConvFlux_${3}"
	newline="CloudConvFlux.fields:         'CloudConvFlux_${3}"
	sed -i "s|$oldline|$newline|g" "${file}"
    fi

    # Retain ProdLoss diagnostics only when the diagnostic target species
    # (the X in ProdXfromY or LossXbyY) is advected.  These diagnostics are
    # ordinary HISTORY fields, so filtering the generated HISTORY.rc is more
    # robust than relying on KPP family membership.  Do not replace these
    # explicit names with Prod_?PRD? or Loss_?LOS?: those wildcards expand all
    # compiled KPP families, even when the target species is not advected.
    firstProdLoss=""
    case "${3}" in
        CH4)
            sed -i "/Prod_CO/d"          "${file}"
            sed -i "/Loss_CO/d"          "${file}"
            sed -i "/ProdCOfromCH4/d"   "${file}"
            sed -i "/ProdCOfromNMVOC/d" "${file}"
            sed -i "/LossCObyOH/d"       "${file}"
            firstProdLoss="Loss_CH4"
            ;;
        CO)
            sed -i "/Loss_CH4/d"       "${file}"
            sed -i "/LossCH4byOH/d"    "${file}"
            sed -i "/LossCH4byCl/d"    "${file}"
            sed -i "/LossCH4byO1D/d"   "${file}"
            sed -i "/LossCH4bySoil/d"  "${file}"
            firstProdLoss="Prod_CO"
            ;;
    esac

    # If deletion removed the original ProdLoss.fields line, promote the
    # first retained diagnostic so that both GCClassic and GCHP syntax remain
    # valid.  Preserve everything after the diagnostic name (including the
    # GCHP component label).
    if [[ -n "${firstProdLoss}" ]]; then
        sed -i -E \
            "s/^[[:space:]]*'${firstProdLoss}/  ProdLoss.fields:          '${firstProdLoss}/" \
            "${file}"
    fi
}


function updateExtData() {

    #========================================================================
    # Removes entries in ExtData.rc for unused species (GCHP only)
    #
    # Arguments:
    # ${1} : List of species to exclude
    # ${2} : Path to the run directory
    #========================================================================

    # File to be modified
    file="${2}/ExtData.rc"

    # Skip if there is no ExtData.rc file (e.g. for GCClassic)
    [[ ! -f "${file}" ]] && return 0

    # If CH4 is in the exclude list, remove CH4 entries.
    isItemInList "CH4" "${1}"
    if [[ $? == 0 ]]; then
        sed -i "/^GHGI_/d"              "${file}"
        sed -i "/^MEX_/d"               "${file}"
        sed -i "/^CAN_/d"               "${file}"
        sed -i "/^GFEI_/d"              "${file}"
        sed -i "/^EDGAR8_CH4_/d"        "${file}"
        sed -i "/CMIP6_CH4_/d"          "${file}"
        sed -i "/CMIP6_BB_CH4/d"        "${file}"
        sed -i "/^UPDATED_GFED4_CH4/d"  "${file}"
        sed -i "/^JPLW_CH4/d"           "${file}"
        sed -i "/^CH4_SEEPS/d"          "${file}"
        sed -i "/^CH4_RES_DAM/d"        "${file}"
        sed -i "/^CH4_RES_SFC/d"        "${file}"
        sed -i "/^CH4_TERMITES/d"       "${file}"
        sed -i "/^CH4_SOIL_UPTAKE/d"    "${file}"
        sed -i "/^\#CH4_/d"             "${file}"
        sed -i "/RCP3PD_CH4/d"          "${file}"
        sed -i "/RCP45_CH4/d"           "${file}"
        sed -i "/RCP60_CH4/d"           "${file}"
        sed -i "/RCP85_CH4/d"           "${file}"
	sed -i "/^EMIS_SF/d"            "${file}"
        sed -i "/^OH_SF/d"              "${file}"
        sed -i "/^MANURE_SF/d"          "${file}"
        sed -i "/^RICE_SF/d"            "${file}"
        sed -i "/^EDGAR_SEASONAL_SF/d"  "${file}"
        sed -i "/^CONUS_/d"             "${file}"
    fi

    # If CO is in the exclude list, remove CO entries
    isItemInList "CO" "${1}"
    if [[ $? == 0 ]]; then
        sed -i "/AEIC19_DAILY_CO /d"    "${file}"  # trailing space required
        sed -i "/AEIC19_MONMEAN_CO /d"  "${file}"  # trailing space required
        sed -i "/APEI_CO/d"             "${file}"
        sed -i "/^CEDS_CO_/d"           "${file}"
        sed -i "/CMIP6_CO_/d"           "${file}"
        sed -i "/\#DICE_/d"             "${file}"
        sed -i "/EDGAR_CO_/d"           "${file}"
        sed -i "/EPA16_CO_/d"           "${file}"
        sed -i "/HTAP_CO_/d"            "${file}"
        sed -i "/RCP3PD_CO /d"          "${file}"
        sed -i "/RCP45_CO /d"           "${file}"
        sed -i "/RCP60_CO /d"           "${file}"
        sed -i "/RCP85_CO /d"           "${file}"
        sed -i "/NEI99_DOW_CO /d"       "${file}"
	sed -i "/LIQFUEL_/d"            "${file}"
    fi

}


function singleOxidantsSpecies() {

    #========================================================================
    # Main function
    #
    # Arguments:
    # ${1} : Species that you wish to retain
    # ${2} : Path to the run directory
    #========================================================================

    # Error check arguments
    if [[ "x${1}" == "x" ]]; then
        echo "Need to pass the species to retain!"
        exit 1
    fi

    # Path to the run directory
    if [[ "x${2}" == "x" ]]; then
        rundir="."
    else
        rundir="${2}"
    fi

    # Get species to include/exclude
    exclude=$(speciesToExclude "${1}")
    include=${1}

    # Update configuration files
    updateGeosChemConfig "${exclude}" "${rundir}"
    updateHemcoConfig    "${exclude}" "${rundir}"
    updateHemcoDiagn     "${exclude}" "${rundir}"
    updateHistory        "${exclude}" "${rundir}" "${include}"
    updateExtData        "${exclude}" "${rundir}"
}


function methaneIsotopologues() {

    #========================================================================
    # Configure a CH4-only prescribed-oxidants run with stable methane
    # isotopologues.  Extension 115 partitions each hierarchy-resolved CH4
    # sector using its gridded d13C and dD source signatures.
    #
    # Arguments:
    # ${1} : Path to the run directory
    #========================================================================

    local rundir="${1:-.}"
    local geos="${rundir}/geoschem_config.yml"
    local hemco="${rundir}/HEMCO_Config.rc"
    local hdiag="${rundir}/HEMCO_Diagn.rc"
    local hist="${rundir}/HISTORY.rc"
    local tmp

    # Transport total CH4 plus all four stable isotopologues.
    sed -i '/^[[:space:]]*- CH4[[:space:]]*$/a\      - C12H4\n      - C13H4\n      - C12H3D\n      - C13H3D' "${geos}"

    # Activate the isotope partitioning extension.  All ordinary inventory
    # records continue to target CH4 only; extension 115 creates the
    # isotopologue-specific mass fluxes online.
    sed -i -E 's|^(115[[:space:]]+CH4_Isotopologues[[:space:]]+:[[:space:]]+)off|\1on |' "${hemco}"

    # The standard carbon restart contains total CH4 but not the four isotope
    # tracers.  Request CH4 explicitly instead of expanding ?ADV? to every
    # transported species.  GEOS-Chem will then initialize each absent isotope
    # from its BackgroundVV value in species_database.yml.
    sed -i -E \
        's|^(\*[[:space:]]+SPC_)[[:space:]]+(.*[[:space:]])SpeciesRst_\?ADV\?|\1CH4           \2SpeciesRst_CH4|' \
        "${hemco}"

    # Mirror every CH4 diagnostic for all four isotopologues and select
    # extension 115, which is where their emissions are added.  Isotopologue
    # biomass burning is explicitly stored in extension category 9.
    tmp=$(mktemp)
    awk '
        {
            print
            if ( $1 ~ /^EmisCH4_/ && $2 == "CH4" ) {
                name = $1
                n = split( "C12H4 C13H4 C12H3D C13H3D", spc, " " )
                for ( i = 1; i <= n; i++ ) {
                    $1 = name
                    sub( /^EmisCH4_/, "Emis" spc[i] "_", $1 )
                    $2 = spc[i]
                    $3 = 115
                    if ( name == "EmisCH4_Total" ) $3 = -1
                    if ( name == "EmisCH4_BiomassBurn" ) $4 = 9
                    longName = $8
                    sub( /^CH4_/, spc[i] "_", longName )
                    printf "%-27s %-7s %5d %4d %4d %3d  %-8s  %s\n",     \
                           $1, $2, $3, $4, $5, $6, $7, longName
                }
            }
        }
    ' "${hdiag}" > "${tmp}" || return 1
    mv "${tmp}" "${hdiag}"

    # Archive isotope concentrations, aggregate chemical losses, isotope
    # ratios, and the emission checks by default.
    sed -i "/^[[:space:]]*'SpeciesConc',/a\             'Isotopologues'," "${hist}"

    sed -i "/'Loss_CH4/a\                             'Loss_C12H4                    ',\n                             'Loss_C13H4                    ',\n                             'Loss_C12H3D                   ',\n                             'Loss_C13H3D                   '," "${hist}"
}
