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
#cd mcarchive

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
      continue
    fi
    cd ..
  fi

done

sleep 6m

cd ${outputDir}
n=$(ls |wc -l)

end=$((SECONDS+600))

until [ $n == $filenb ]
  do
    n=$(ls |wc -l)
    sleep 2s
    if [ $SECONDS -gt $end ]; then
      echo "Le script est bloqué depuis 10m, je m'arrête"
      break
    fi
  done

filenamemerged=outputfile_studyTracks_merged.root
filestomerge=outputfile_studyTracks_???.root
hadd -f $filenamemerged ${filestomerge}

#line 46: cd: 03-11-2021-18:00/histoResults_03-11-2021-18:00: No such file or directory
