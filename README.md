# feedback-plate-analyzer
Python scripts to analyze the sound of finger-piezo feedback plates, making a pitch map of the plate

Requirements:
scipy, numpy, matplotlib, json, sounddevice, colorcet, mayavi (for 3D models)

To use:
1. Record a series of sound samples from the plate. The samples will be laid out in a grid in the visualization, so best to sample sounds in a regular grid of positions on the plate, for example a 9x9 grid. If you produce a single sustained sound on each sample position, the software should be able to segment the recording automatically.
To test the package you can use the supplied sound file (210629_1334_30x30_plate_9x9_soft_finger_5ptape_reverb_handheld.wav)

2. Run seg_peaks_export.py on the sound file, like this
python seg_peaks_export.py -sf 210629_1334_30x30_plate_9x9_soft_finger_5ptape_reverb_handheld
(omit the .wav extension as it is added internally when reading the sound file)
The python script has additional command line parameters that allow you to set
minimum and maximum frequency of the analysis, decibel threshold for the segmentation, minimum segment length for the segmentation, rms analysis framelength (might also help adjusting the segmentation), minimum spectral peak distance for the analysis, number of columns in the layout grid, and other parameters. Run python seg_peaks_export.py -h for a full list of parameters.

3. Run visualize_circles.py to create the 2D map with one circle per audio sample taken on the plate. Like this:
python visualize_circles.py -sf 210629_1334_30x30_plate_9x9_soft_finger_5ptape_reverb_handheld
The python script has additional command line parameters that allow you to set display parameters, also some that will override settings done in the precious step (for example spectral peak distance). Run visualize_circles.py -h for a full list of parameters.

4. Optionally run visualize_3d.py to create a 3D map, like this
python visualize_3d.py -sf 210629_1334_30x30_plate_9x9_soft_finger_5ptape_reverb_handheld
The python script has additional command line parameters that allow you to set display parameters. Run visualize_3d.py -h for a full list of parameters.
