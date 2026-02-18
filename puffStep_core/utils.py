import sys, datetime

def newmsg(x):
    sys.stderr.write(str(datetime.datetime.now()) +": .. " + x + "...\n")

def bdgmsg(x, collapsed=False):
    if not collapsed:
        return sys.stderr.write(str(datetime.datetime.now()) +": ..printing " + x + " as expanded single-step bedGraph...\n")
    else:
        return sys.stderr.write(str(datetime.datetime.now()) +": ..printing " + x + " as collapsed var-step bedGraph...\n")
