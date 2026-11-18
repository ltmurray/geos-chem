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

# CH4 tracers used to separate emissions by source sector
CH4_SECTOR_SPECIES=(CH4_OIL CH4_GAS CH4_COL CH4_LIV CH4_LDF CH4_WST CH4_RIC \
                    CH4_OTA CH4_BBN CH4_WTL CH4_SEE CH4_LAK CH4_TER CH4_RES)


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


function setupTaggedCH4Sectors() {

    #========================================================================
    # Adds sector-tagged CH4 tracers to a CH4-only oxidants simulation.
    # Total CH4 remains an independent tracer for closure checks.  Source
    # sector tracers use the category assignments from the historical tagged
    # CH4 simulation.  Soil uptake is a first-order loss, not a source tracer.
    #
    # Arguments:
    # ${1} : Path to the run directory
    #========================================================================

    rundir="${1}"
    config="${rundir}/geoschem_config.yml"
    hemco="${rundir}/HEMCO_Config.rc"
    diagn="${rundir}/HEMCO_Diagn.rc"
    history="${rundir}/HISTORY.rc"

    # Add the source-sector tracers after total CH4.
    for ((n=${#CH4_SECTOR_SPECIES[@]}-1; n>=0; n--)); do
        spc="${CH4_SECTOR_SPECIES[n]}"
        sed -i "/      - CH4$/a\\      - ${spc}" "${config}"
    done

    # Add allocation factors for combined Canadian oil-and-gas emissions.
    # HEMCO does not support complements of scalar scale factors, so these
    # values must be kept complementary.  They can later be replaced by
    # spatially varying oil- and gas-fraction fields.
    sed -i "/1 NEGATIVE -1.0/a\\
# Canadian oil-and-gas fractions (CAN_OIL_FRACTION + CAN_GAS_FRACTION = 1)\\
2901 CAN_OIL_FRACTION 0.5 - - - xy 1 1\\
2902 CAN_GAS_FRACTION 0.5 - - - xy 1 1" "${hemco}"

    # Duplicate each total-CH4 base-emission line into its historical source
    # sector.  Split combined-sector fields between tracers with explicitly
    # named HEMCO scale factors.  Category 9 biomass burning is handled
    # separately below because the default inventory is supplied by a HEMCO
    # extension.
    awk '
        function tag(cat) {
            if (cat == "1")  return "CH4_OIL"
            if (cat == "2")  return "CH4_GAS"
            if (cat == "3")  return "CH4_COL"
            if (cat == "4")  return "CH4_LIV"
            if (cat == "5")  return "CH4_LDF"
            if (cat == "6")  return "CH4_WST"
            if (cat == "7")  return "CH4_RIC"
            if (cat == "8")  return "CH4_OTA"
            if (cat == "9")  return "CH4_BBN"
            if (cat == "10") return "CH4_WTL"
            if (cat == "11") return "CH4_SEE"
            if (cat == "12") return "CH4_LAK"
            if (cat == "13") return "CH4_TER"
            if (cat == "15") return "CH4_RES"
            return ""
        }
        function allocation(name, cat) {
            if (name ~ /^CAN_OIL_GAS_/ && cat == "1") return "2901"
            if (name ~ /^CAN_OIL_GAS_/ && cat == "2") return "2902"
            return ""
        }
        {
            print
            if ($1 == "112" && $2 == "GFAS_CH4_2D") {
                $2 = "GFAS_CH4_BBN_2D"
                print
                next
            }
            if ($1 != "0" || $2 == "CH4_Emis_Prior") next
            target = ""
            for (n = 3; n <= NF; n++) {
                if ($n == "CH4") {
                    target = n
                }
            }
            if (target == "") next
            baseName = $2
            baseScale = $(target+1)
            nCat = split($(NF-1), cats, "/")
            for (i = 1; i <= nCat; i++) {
                spc = tag(cats[i])
                if (spc == "") continue
                $2 = baseName "_TAG_" cats[i]
                $target = spc
                $(NF-1) = cats[i]
                if (nCat > 1) {
                    scale = allocation(baseName, cats[i])
                    if (scale == "") {
                        printf "No tagged-sector allocation factor for %s category %s\n", \
                               baseName, cats[i] > "/dev/stderr"
                        exit 2
                    }
                    $(target+1) = baseScale "/" scale
                } else {
                    $(target+1) = baseScale
                }
                print
            }
        }
    ' "${hemco}" > "${hemco}.tmp"
    mv "${hemco}.tmp" "${hemco}"

    # Add the biomass-burning tracer to the GFED/GFAS extension species list.
    sed -i 's/CH4\/CO/CH4\/CH4_BBN\/CO/g' "${hemco}"

    # Archive the total emission flux into each sector tracer as an independent
    # check on the category-based total-CH4 emission diagnostics.
    for spc in "${CH4_SECTOR_SPECIES[@]}"; do
        printf "Emis%s_Total  %s  -1  -1  -1  2  kg/m2/s  %s_emissions_from_all_sources\n" \
               "${spc}" "${spc}" "${spc}" >> "${diagn}"
    done

    # Add per-sector concentration-loss diagnostics.  Soil loss is 2-D;
    # oxidant loss diagnostics are 3-D.
    for spc in "${CH4_SECTOR_SPECIES[@]}"; do
        sed -i "/'LossCH4bySoil/a\\                             'Loss${spc}bySoil               '," "${history}"
        sed -i "/'LossCH4byO1D/a\\                             'Loss${spc}byO1D                '," "${history}"
        sed -i "/'LossCH4byCl/a\\                             'Loss${spc}byCl                 '," "${history}"
        sed -i "/'LossCH4byOH/a\\                             'Loss${spc}byOH                 '," "${history}"
    done
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
    if [[ "x${1}" == "xCH4_SECTORS" ]]; then
        include="CH4"
    else
        include=${1}
    fi
    exclude=$(speciesToExclude "${include}")

    # Update configuration files
    updateGeosChemConfig "${exclude}" "${rundir}"
    updateHemcoConfig    "${exclude}" "${rundir}"
    updateHemcoDiagn     "${exclude}" "${rundir}"
    updateHistory        "${exclude}" "${rundir}" "${include}"
    updateExtData        "${exclude}" "${rundir}"

    if [[ "x${1}" == "xCH4_SECTORS" ]]; then
        setupTaggedCH4Sectors "${rundir}"
    fi
}
