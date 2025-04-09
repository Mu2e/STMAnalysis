mu2e -c Offline/Sources/test/FromSTMTestBeamData.fcl -s /pnfs/mu2e/tape/phy-raw/raw/mu2e/STM_ELBE_2022/HPGeRaw/dat/$1/*/*.dat -n 1000
mu2e -c Production/JobConfig/stm/data/zeroSuppression.fcl -s dig.owner.description.version.sequencer.art
mu2e -c Production/JobConfig/stm/data/mwd.fcl -s dig.owner.STMZP.version.sequencer.art
mu2e -c Production/JobConfig/stm/PlotSTMDigisSpectrum.fcl -s dig.owner.STMMWD.version.sequencer.art
root -l Offline/STMReco/src/digiGain.C\(\"H\"\)
