scriptName="makePUReWeightJSON.py"

python "${scriptName}" --nominal=MyDataPileupHistogram_16preVFP.root --up=MyDataPileupHistogram_16preVFP_up.root --down=MyDataPileupHistogram_16preVFP_down.root --makePlot -o 2016preVFP_UL/puWeights.json --gzip --mcprofile=2016UL_25ns --format=correctionlib --name=Collisions16_UltraLegacy_goldenJSON

python "${scriptName}" --nominal=MyDataPileupHistogram_16postVFP.root --up=MyDataPileupHistogram_16postVFP_up.root --down=MyDataPileupHistogram_16postVFP_down.root --makePlot -o 2016postVFP_UL/puWeights.json --gzip --mcprofile=2016UL_25ns --format=correctionlib --name=Collisions16_UltraLegacy_goldenJSON

python "${scriptName}" --nominal=MyDataPileupHistogram_17.root --up=MyDataPileupHistogram_17_up.root --down=MyDataPileupHistogram_17_down.root --makePlot -o 2017_UL/puWeights.json --gzip --mcprofile=2017UL_25ns --format=correctionlib --name=Collisions17_UltraLegacy_goldenJSON

python "${scriptName}" --nominal=MyDataPileupHistogram_18.root --up=MyDataPileupHistogram_18_up.root --down=MyDataPileupHistogram_18_down.root --makePlot -o 2018_UL/puWeights.json --gzip --mcprofile=2018UL_25ns --format=correctionlib --name=Collisions18_UltraLegacy_goldenJSON

#for file in */puWeights.json.gz
#do
#    mv $file ../../POG/LUM/$file
#done

