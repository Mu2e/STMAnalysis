mu2e -c Offline/Sources/test/FromSTMTestBeamData.fcl -s /pnfs/mu2e/tape/phy-raw/raw/mu2e/STM_ELBE_2022/LaBrRaw/dat/$1/*/*.dat -n 3
mu2e -c Offline/STMReco/fcl/pedestalSubtraction.fcl -s dig.owner.description.version.sequencer.art
mu2e -c Production/JobConfig/stm/data/zeroSuppressionLaBr.fcl -s dig.owner.PedSubWaveforms.version.sequencer.art --nskip 2
mu2e -c Production/JobConfig/stm/data/mwdLaBr.fcl -s dig.owner.STMZP.version.sequencer.art
