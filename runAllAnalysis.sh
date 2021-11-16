#!/bin/bash
dirname=$1
#Nom du dossier dans lequel il y aura les dossiers à analyser
newdirname=${dirname%?????????????????????????}
#nom du dossier mais simplifié
outputDirName=histoResults_${newdirname}
baseDir=$PWD


mkdir ${newdirname}
cd ${newdirname}
alien_cp /alice/cern.ch/user/a/aliperf/alibi_nightlies/O2DPG_pp_minbias_testbeam.sh/${dirname}/mcarchive.tar.gz mcarchive.tar.gz

tar -xvf mcarchive.tar.gz
rm mcarchive.tar.gz

mkdir ${outputDirName}

outputDir=${PWD}/${outputDirName}


echo "I am working in the directory " $PWD

filenb=0

for i in $(seq 25)
do
  if [[ -d "tf$i" ]]
  then
    echo "entering directory tf$i"
    cd tf$i


    cp $baseDir/StudyMFTTracks.C .
    cp $baseDir/StudyMFTPurity.C .
    cp $baseDir/MFTdictionary.bin .

    i2=$(printf "%03d" $i)

    filenameout=${outputDir}/outputfile_studyTracks_${i2}.root
    if [[ -f "mfttracks.root" && -f "mftclusters.root" ]]
    then
    #root -l -b -q 'StudyMFTTracks.C('\"${filenameout}\"','\"sgn_${i}_Kine.root\"')' && root -l -b -q 'StudyMFTPurity.C('\"${filenameout}\"','\"sgn_${i}_Kine.root\"')' &
      root -l -b -q 'StudyMFTTracks.C('\"${filenameout}\"','\"sgn_${i}_Kine.root\"')' &
      filenb=$((filenb+1))
    else
      echo "The directory tf$i does not contain mfttracks.root or mftclusters.root"
      continue
    fi
    cd ..
  fi

done

sleep 6m

if [ $filenb -eq 0 ]; then
  echo "No files mfttracks.root or mftclusters.root, exiting"
  exit 1
fi

cd ${outputDir}
n=$(ls |wc -l)

end=$((SECONDS+600))

until [ $n == $filenb ]
  do
    n=$(ls |wc -l)
    sleep 2s
    if [ $SECONDS -gt $end ]; then
      echo "The script has been stuck for 10m, let's break it"
      break
    fi
  done

filenamemerged=outputfile_studyTracks_merged.root
filestomerge=outputfile_studyTracks_???.root
hadd -f $filenamemerged ${filestomerge}

cp $baseDir/EvalEffAndPurity.C .

root -l -b -q 'EvalEffAndPurity.C('\"${filenamemerged}\"','\"${newdirname}_EvalEffAndPurity\"')'
