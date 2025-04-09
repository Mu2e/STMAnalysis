rm mwdSize.log
rootFile=""
rm dig.owner.STMMWD.version.sequencer.art
for fclfile in `ls mwd_M*.fcl` ; do
    echo $fclfile
    mu2e -c $fclfile -s ./dig.owner.STMZP.version.sequencer.art
    done
