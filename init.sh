#!/bin/bash
pathadd() {
    if [[ ":$PATH:" != *":$1:"* ]]; then
        PATH="${PATH:+"$PATH:"}$1"
        export PATH
    fi
}
# get the directory of the script
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]:-$0}")" &>/dev/null && pwd 2>/dev/null)"
if ! command -v lsb_release &> /dev/null
then
    source /etc/os-release
    distro=$NAME
    os_version=$VERSION_ID
else
    distro=$(lsb_release -i | cut -f2)
    os_version=$(lsb_release -r | cut -f2)
fi
distro=${distro//[[:space:]]/}
distro="${distro//Linux/}"
distro="${distro//linux/}"
echo "Setting up CROWN for $distro Version $os_version"
# check if the distro is centos
if [[ "$distro" == "CentOS" ]]; then
    # if the first number of os_version is a 7, we are on centOS 7
    if [[ ${os_version:0:1} == "7" ]]; then # if uname -a | grep -E 'el7' -q
        # source /cvmfs/sft-nightlies.cern.ch/lcg/views/dev3/latest/x86_64-centos7-gcc11-opt/setup.sh
        # source /cvmfs/sft-nightlies.cern.ch/lcg/views/dev3/latest/x86_64-centos7-clang12-opt/setup.sh
        # source /cvmfs/sft-nightlies.cern.ch/lcg/views/dev3/latest/x86_64-centos7-gcc11-dbg/setup.sh
        echo "CentOS 7 is EOL, running on LCG 105, support will be dropped soon"
        source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-centos7-gcc11-opt/setup.sh
    else
        echo "Unsupported CentOS version, exiting..."
        return 0
    fi
elif [[ "$distro" == "RedHatEnterprise" || "$distro" == "Alma" || "$distro" == "Rocky" ]]; then
    if [[ ${os_version:0:1} == "8" ]]; then # elif uname -a | grep -E 'el8' -q
        echo "Unsupported CentOS version, exiting..."
        return 0
    elif [[ ${os_version:0:1} == "9" ]]; then # elif uname -a | grep -E 'el8' -q
        source /cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-dbg/setup.sh
    else
        echo "Unsupported CentOS version, exiting..."
        return 0
    fi
elif [[ "$distro" == "Ubuntu" ]]; then
    if [[ ${os_version:0:2} == "20" ]]; then
        source /cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-ubuntu2004-gcc9-opt/setup.sh
    elif [[ ${os_version:0:2} == "22" ]]; then
        source /cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-ubuntu2204-gcc11-opt/setup.sh
    else
        echo "Unsupported Ubuntu version, exiting..."
        return 0
    fi
else
    echo "You are not running on CentOS or Ubuntu, exiting..."
    return 0
fi
# add ~/.local/bin to path if it is not already there
pathadd "${HOME}/.local/bin/"
# set the cmake generator to Ninja
# export CMAKE_GENERATOR="Ninja"
export CMAKE_GENERATOR="Unix Makefiles"
# set the compiler optimization for cling to O2, this
# will result in about 20% faster JIT for the snapshot generation
export EXTRA_CLING_ARGS='-O2'

# clone a given analysis if an argument is given
if [ -z "$1" ]; then
    echo "No configuration clone"
else
    if [[ "$1" == "tau" && ! -d "${SCRIPT_DIR}/analysis_configurations/tau" ]]; then
        echo "Cloning analysis tau into ${SCRIPT_DIR}/analysis_configurations/tau"
        git clone git@github.com:KIT-CMS/TauAnalysis-CROWN.git "${SCRIPT_DIR}/analysis_configurations/tau"
    elif [[ "$1" == "run3analysis" && ! -d "${SCRIPT_DIR}/analysis_configurations/run3analysis" ]]; then
        echo "Cloning analysis run3analysis into ${SCRIPT_DIR}/analysis_configurations/run3analysis"
        git clone https://github.com/KIT-CMS/Run3Analysis-CROWN "${SCRIPT_DIR}/analysis_configurations/run3analysis"
    elif [[ "$1" == "whtautau" && ! -d "${SCRIPT_DIR}/analysis_configurations/whtautau" ]]; then
        echo "Cloning analysis whtautau into ${SCRIPT_DIR}/analysis_configurations/whtautau"
        git clone git@github.com:KIT-CMS/WHTauTauAnalysis-CROWN.git "${SCRIPT_DIR}/analysis_configurations/whtautau"
    elif [[ "$1" == "boosted_h_tautau" && ! -d "${SCRIPT_DIR}/analysis_configurations/boosted_h_tautau" ]]; then
        echo "Cloning analysis boosted_h_tautau into ${SCRIPT_DIR}/analysis_configurations/boosted_h_tautau"
        git clone git@github.com:KIT-CMS/BoostedHiggsTauTauAnalysis-CROWN.git "${SCRIPT_DIR}/analysis_configurations/boosted_h_tautau"
    elif [[ "$1" == "s" && ! -d "${SCRIPT_DIR}/analysis_configurations/s" ]]; then
        echo "Cloning analysis s-channel into ${SCRIPT_DIR}/analysis_configurations/s"
        git clone git@github.com:nfaltermann/CROWNs.git "${SCRIPT_DIR}/analysis_configurations/s"
    fi
fi
