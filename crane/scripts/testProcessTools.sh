#/bin/bash

# A script that aims to test all the main functionality of the processing tools.
#

# variables
DPATH=/home/workspace/ccalmr46/pturner/db29/
ST=2012-4-30
ET=2012-5-18
ET_SHORT=2012-5-1
runTag=run130
STAFILE="" #"-t real_stations.sta"
i=1; l=`printf "%03d" $i` # for log file numbering

# *** skillExtract
if [ -e obs ]; then
  rm -rf obs
fi
if [ -e $runTag ]; then
  rm -rf $runTag
fi
#  - all mode
skillExtract -A -r $runTag -s $ST -e $ET -d $DPATH$runTag/outputs/ &> log"$l"_skillExtract1.txt
i=$[i+1];l=`printf "%03d" $i`
rm -rf obs

#  - obs mode
skillExtract -r $runTag -s $ST -e $ET -d $DPATH$runTag/outputs/ &> log"$l"_skillExtract2.txt
i=$[i+1];l=`printf "%03d" $i`

#  - online obs mode (do not store obs on disk)
skillExtract -N -r $runTag -s $ST -e $ET -d $DPATH$runTag/outputs/ &> log"$l"_skillExtract3.txt
i=$[i+1];l=`printf "%03d" $i`

# *** generateSaltIntrusion
generateSaltIntrusion -r $runTag -d $DPATH$runTag/outputs/ -s $ST -e $ET &> log"$l"_generateSaltIntrusion.txt
i=$[i+1];l=`printf "%03d" $i`

# *** skillPlot

#  - without observations
if [ -e obs_old ]; then
  rm -rf obs_old
fi
mv obs obs_old
skillPlot -s $ST -e $ET -o $runTag/images $runTag &> log"$l"_skillPlot1.txt
i=$[i+1];l=`printf "%03d" $i`

#  - online obs mode  (do not read obs from disk, fetch from netcdfcache)
skillPlot -N -s $ST -e $ET -o $runTag/images $runTag &> log"$l"_skillPlot2.txt
i=$[i+1];l=`printf "%03d" $i`

#  - normal mode
if [ -e obs ]; then
  rm -rf obs
fi
mv obs_old obs
skillPlot -s $ST -e $ET -o $runTag/images $runTag &> log"$l"_skillPlot3.txt
i=$[i+1];l=`printf "%03d" $i`

# *** extractStation

# *** extractTransect
extractTransect -r $runTag -d $DPATH$runTag/outputs/ -v salt -t nc_auv_along.bp -n ncAUValong -s $ST -e $ET &> log"$l"_extractTransect1.txt
i=$[i+1];l=`printf "%03d" $i`

extractTransect -r $runTag -d $DPATH$runTag/outputs/ -v salt -t nc_auv_east.bp -n ncAUVeast -s $ST -e $ET &> log"$l"_extractTransect2.txt
i=$[i+1];l=`printf "%03d" $i`

extractTransect -r $runTag -d $DPATH$runTag/outputs/ -v salt -t nc_auv_west.bp -n ncAUVwest -s $ST -e $ET &> log"$l"_extractTransect3.txt
i=$[i+1];l=`printf "%03d" $i`

# *** extractSlab
extractSlab -r $runTag -d $DPATH$runTag/outputs/ -v salt -n slab -k 1 -s $ST -e $ET_SHORT &> log"$l"_extractSlab1.txt
i=$[i+1];l=`printf "%03d" $i`

# *** fetch AUV obs data
fetchAUVData -v salt -i 46 &> log"$l"_fetchAUVData1.txt
i=$[i+1];l=`printf "%03d" $i`

# *** extractTrack
TMP="AUV/data/track/AUV46"
extractTrack -r $runTag -d $DPATH$runTag/outputs/ -v salt -n AUV46 -i `ls $TMP*` &> log"$l"_extractTract1.txt
i=$[i+1];l=`printf "%03d" $i`

# *** makeSlabPlots
makeSlabPlots -r $runTag -s $ST -e $ET_SHORT -o $runTag/images/slab/  -b [326000,352000,284000,302000] $runTag/data/slab/slab_salt_s1_*.nc &> log"$l"_makeSlabPlots.txt
i=$[i+1];l=`printf "%03d" $i`

# *** makeTimeSeriesPlots

# *** makeTrackPlots
TMP="$runTag/data/track/AUV46"
makeTrackPlot -o $runTag/images/track/ -r $runTag `ls $TMP*` &> log"$l"_makeTrackPlot.txt
i=$[i+1];l=`printf "%03d" $i`

# *** makeTransects
TMP="$runTag/data/transect/ncAUValong"
makeTransects -s $ST -e $ET_SHORT -o $runTag/images/transect/ `ls $TMP*` &> log"$l"_makeTransects.txt
i=$[i+1];l=`printf "%03d" $i`

# *** makeAUVTransect
makeAUVTransect -r $runTag -v salt -o $runTag/images/auv -s $ST -e $ET_SHORT  &> log"$l"_makeAUVTransect.txt
i=$[i+1];l=`printf "%03d" $i`

echo '* TESTS FINISHED'
