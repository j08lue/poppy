import os

def parse_coords(s):
    result = map(float,(s.replace('m','-')).split(','))
    if len(result) == 1:
        return result[0]
    else:
        return result

def parse_pattern(pattern):
    if os.path.isdir(pattern):
        pattern = os.path.join(pattern,'*.pop.h.????-??.nc')
    return pattern
