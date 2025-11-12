#!/bin/bash
# Usage: ./to_mp4.sh <iteration>
ITERATION=$1
ffmpeg -framerate 30 -pattern_type glob -i "output/${ITERATION}/vmag_*.png" -c:v libx264 -pix_fmt yuv420p output/${ITERATION}.mp4
# gif
ffmpeg -i output/${ITERATION}.mp4 -vf "fps=10,scale=720:-1:flags=lanczos" -c:v gif output/${ITERATION}.gif