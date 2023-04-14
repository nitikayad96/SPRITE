import sprite_utils as su

if __name__ == "__main__":

    cal_load_iter = int(sys.argv[1])

    # wait 4 seconds, roughly the amount of time required for the correlator code to initialize
    # then begin hot load slewing
    time.sleep(4)
    for i in range(cal_load_iter):

        # switch hot load into place for 15 seconds
        su.switch_cal_load()
        # wait 5 minutes for next cycle
        if i + 1 < cal_load_iter:
            time.sleep(5 * 60 - 15)
        else:
            break