from music21 import *

def getIntervals(xmlFile):
    '''
    get intervals
    :param xmlFile:
    :return:
    '''
    xmlScore = converter.parse(xmlFile)

    intervals = []

    for part in xmlScore:
        try:
            if part.flat.notes[0].hasLyrics():
                p = part.flat.notes
                for pp in p.elements:
                    if pp.duration.isGrace:
                        p.remove(pp)

                for ii in range(len(p)-1):
                    intervals.append(interval.notesToChromatic(p[ii],p[ii+1]).semitones)
        except:
            pass

    return intervals

# sAlt.show() # show first 5 measures