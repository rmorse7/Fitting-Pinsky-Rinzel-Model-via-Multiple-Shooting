#!/bin/bash
# declare STRING variable
STRING="Setting up environment variables..."
#print variable on a screen
echo $STRING

export LD_LIBRARY_PATH=/home/gabbiani/NAG/cll6i261dl/lib:/home/gabbiani/NAG/cll6i261dl/rtl/intel64
export NAG_KUSARI_FILE=/home/gabbiani/NAG/cll6i261dl/license.lic


echo "done"
