#!/bin/bash
echo "bash: creating video from .png files.."
ffmpeg -r 10 -s 600x600 -i %d.png  -crf 25 buo.mp4 
echo "bash: creating video from .png files done."
