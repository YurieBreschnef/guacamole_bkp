#!/bin/bash
echo "bash: creating video from .png files.."
ffmpeg -r 20 -s 1200x600 -i %d.png  -crf 25 combo.mp4 
echo "bash: creating video from .png files done."
