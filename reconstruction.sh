#. digisLaBr.sh $1
python3 STMCalibration/calibrationSTM.py
mu2e -c Production/JobConfig/stm/data/makeSTMHits.fcl -s dig.owner.STMMWD.version.sequencer.art
mu2e -c Production/JobConfig/stm/plotSTMSpectrum.fcl -s mcs.owner.STM.version.sequencer.art
root -l MWDResolution.C\(\"stmSpectrum.root\"\)
