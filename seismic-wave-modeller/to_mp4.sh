#!/bin/bash

ITERATION=$1
ffmpeg -framerate 30 -pattern_type glob -i 'output/${ITERATION}/vmag_*.png' -c:v libx264 -pix_fmt yuv420p output/${ITERATION}.mp4
# gif
ffmpeg -i output/${ITERATION}.mp4 -c:v libx264 -pix_fmt yuv420p output/${ITERATION}.gif