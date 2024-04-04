from numpy import arange, hstack

def inclusive_range(start, end, step):
    r = arange(start, end, step)
    if r[-1] + step == end:
        r = hstack((r, [end]))
    return r
