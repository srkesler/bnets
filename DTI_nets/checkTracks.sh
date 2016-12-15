#!/bin/bash
# Shelli Kesler 8/26/16
# Batch save screen capture of tractography results from TrackVis for quality check
# Requires GraphicsMagick (http://www.graphicsmagick.org)

export PATH=$PATH:/Applications/TrackVis.app/Contents/MacOS

SUBJS=$(cat /fullpath/Paths.txt) #txt file containing full path to each subject folder e.g. /Desktop/Study/Subj1

mkdir /fullpath/Montage #path where screenshots will be stored for combining

for subjdir in ${SUBJS}
do 
	cd ${subjdir}
	titlestr=${subjdir: -6} #6 = no. characters from end of path used to label screencap
	track_vis dti.trk -ny 65 -title ${titlestr} -na -sc ${titlestr}trk.png 1 14 14 3 -disable_log #65 = coronal frame; 14 = image out of 360 rotation 
	cp -a ${titlestr}trk*.png ../Montage
done

cd ../Montage
gm convert *.png checkDTItrks.pdf 
