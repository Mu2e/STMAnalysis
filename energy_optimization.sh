rm mwdResolution.log
rootFile=""
rm dig.owner.STMMWD.version.sequencer.art
for fclfile in `ls mwd_M*.fcl` ; do
    echo $fclfile
    mu2e -c $fclfile -s ./dig.owner.STMZP.version.sequencer.art
    rootFile=`echo $fclfile | cut -d . -f 1`
    echo $rootFile
    rootFile=${rootFile}".root"
    echo $rootFile
    mu2e -c Production/JobConfig/stm/PlotSTMDigisSpectrum.fcl -s dig.owner.STMMWD.version.sequencer.art
    mu2e -c Production/JobConfig/stm/data/makeSTMHits.fcl -s dig.owner.STMMWD.version.sequencer.art
    mu2e -c Production/JobConfig/stm/plotSTMSpectrum.fcl -s mcs.owner.STM.version.sequencer.art --TFileName $rootFile
    root -l -b -q MWDResolution.C\(\"$rootFile\"\)
    done
