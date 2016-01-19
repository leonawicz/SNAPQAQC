#! /bin/bash

model=$1 # e.g., cccma_cgcm3_1.sresa1b, model directory
period=$2 # historical or not (i.e., projected), indicate if years spexcified below are historical CRU 3.2 data

if [ "$period" == "historical" ]
then
 modelOut=CRU32 # renaming destination model directory 
 yr1=1950
 yr2=2007
else
 modelOut=$model # leaving destination directory same as source
 yr1=2008
 yr2=2100
fi

inDir=/atlas_scratch/apbennett/IEM/FinalCalib/$model/Maps # source Maps directory
outDir=/atlas_scratch/mfleonawicz/alfresco/IEM/outputs/FinalCalib/$modelOut/Maps # destination Maps directory
mkdir -p $outDir

for year in `seq $yr1 $yr2`;
do
 out=$outDir/$year
 mkdir $out
 cp $inDir/Age*$year.tif $out/
 cp $inDir/FireScar*$year.tif $out/
 cp $inDir/Veg*$year.tif $out/
done
