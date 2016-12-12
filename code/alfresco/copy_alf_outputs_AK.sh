#! /bin/bash

rcp=$1 #e.g. rcp60
model=$2 # e.g., cccma_cgcm3_1.sresa1b, model directory
period=$3 # historical or not (i.e., projected), indicate if years specified below are historical CRU 3.2 data
gbm=$4 # 3m or 5m

if [ "$period" == "historical" ]
then
 modelOut=CRU32 # renaming destination model directory 
 yr1=1900
 yr2=2013
else
 modelOut=$model # leaving destination directory same as source
 yr1=2014
 yr2=2099
fi

inDir=/big_scratch/shiny/Runs_Statewide/paul.duffy_at_neptuneinc.org/*$rcp*$model/Maps # source Maps directory
outDir=/atlas_scratch/mfleonawicz/alfresco/CMIP5_Statewide/outputs/$gbm/$rcp.$modelOut/Maps # destination Maps directory
mkdir -p $outDir

for year in `seq $yr1 $yr2`;
do
 out=$outDir/$year
 #in=$inDir/*$year # are inputs already organized into year directories?
 in=$inDir # if not.
 mkdir $out
 cp $in/Age*$year.tif $out/
 cp $in/FireScar*$year.tif $out/
 cp $in/Veg*$year.tif $out/
done
