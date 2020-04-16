from time import sleep
def makevideo(ts_min, ts_max, frames, secs0, secs, output_prefix, output_format):
    length = ts_max - ts_min
    interval = length / (frames-1) # time interval between each frame
    t = ts_min # set the cursor to the beginning
    setVisible(g.filter(S <= t).filter(E >= t).filter(degree >= 1))
    sleep(secs0)
    exportGraph("%s%s%s" % (output_prefix, 0, output_format))
    for i in range(1, frames): # let's start to count from image 0 to image FRAMES
        t += interval 
        setVisible(g.filter(S <= t).filter(E >= t).filter(degree >= 1))
        sleep(secs)
        exportGraph("%s%s%s" % (output_prefix, i, output_format))

