Producing a video with Gephi and ffmpeg
=======================================
0. in Gephi, install the Scripting Plugin <https://github.com/gephi/gephi/wiki/Scripting-Plugin>
1. copy "first.py" and "second.sh" to a fresh directory, cd to that directory, start Gephi from there
2. in Gephi: 
    * open graph without (!) merging edges
    * under "Preview", select proper layout style
    * under "Overview", start ForceAtlas2 visualization and keep it running
    * open "Console" and type ``execfile("./first.py")``
    * run a command like ``makevideo(ts_min=0.0, ts_max=0.05, frames=50, secs0=5, secs=1, output_prefix="./frame_", output_format=".png")`` and wait for it to finish, using the proper parameters (start time, end time, number of frames to be generated, seconds to wait for layout to converge before the first frame, seconds to wait before each further frame)
3. in a shell, run ``bash second.sh`` to produce the movie ``out.mp4`` from the images ``frame_XXX.png``
