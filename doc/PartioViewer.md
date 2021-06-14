# PartioViewer

The simulators can export the particle simulation data using the partio file format. The PartioViewer can read such a file and render the particle data using OpenGL. This tool is able to handle multiphase data and rigid body data. It can create image sequences and movies (using ffmpeg).

To visualize a sequence of partio files or a single file, call (the index in the file name is used for the sequence): 
```
PartioViewer fluid_data_1.bgeo
```

This tool is also able to read a complete output directory:
```
PartioViewer output/DamBreakModel
```
In this case the tool searches for the partio files of multiple phases in the subdirectory "partio" and for rigid body data in "rigid_bodies".

Note: To generate videos you must tell PartioViewer where it can find the ffmpeg executable.

#### Command line options:

* -h, --help: Print help
* --renderSequence: Render a sequence from startFrame to endFrame as jpeg.
* --renderVideo: Render a sequence from startFrame to endFrame as video.This function requires ffmpeg which must be in the PATH or the ffmpegPath parameter must be set.
* --noOverwrite: Do not overwrite existing frames when using --renderSequence option. Existing frames are not loaded at all which accelerates the image sequence generation.
* -o, --outdir arg: Output directory for images
* --rbData arg: Rigid body data to visualize (bin file)
* --ffmpegPath arg: Path of the ffmpeg excutable.
* --width arg: Width of the image in pixels. (default: 1024)
* --height arg: Height of the image in pixels. (default: 768)
* --fps arg: Frame rate of video. (default: 25)
* -r, --radius arg: Particle radius (default: 0.025)
* -s, --startFrame arg: Start frame (only used if value is >= 0) (default: -1)
* -e, --endFrame arg: End frame (only used if value is >= 0) (default: -1)
* --colorField arg: Name of field that is used for the color. (default: velocity)
* --colorMapType arg: Color map (0=None, 1=Jet, 2=Plasma) (default: 1)
* --renderMinValue arg: Min value of field. (default: 0.0)
* --renderMaxValue arg: Max value of field. (default: 10.0)
* --camPos arg: Camera position (e.g. --camPos "0 1 5") (default: 0 3 10)
* --camLookat arg: Camera lookat (e.g. --camLookat "0 0 0") (default: 0 0 0)

#### Hotkeys

* Space: pause/contiunue simulation
* r: reset simulation
* w: wireframe rendering of meshes
* i: print all field information of the selected particles to the console
* s: save current frame as jpg image
* v: generate video 
* j: generate image sequence
* +: step to next frame
* -: step to previous frame
* ESC: exit