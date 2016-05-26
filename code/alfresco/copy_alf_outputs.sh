#! /bin/bash

model=$1 # e.g., cccma_cgcm3_1.sresa1b, model directory
period=$2 # historical or not (i.e., projected), indicate if years specified below are historical CRU 3.2 data

if [ "$period" == "historical" ]
then
 modelOut=CRU32 # renaming destination model directory 
 yr1=1900
 yr2=1949
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
 in=$inDir/*$year
 mkdir $out
 cp $in/Age*$year.tif $out/
 cp $in/FireScar*$year.tif $out/
 cp $in/Veg*$year.tif $out/
done
