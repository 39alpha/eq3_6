import re
import sys

lineregex = re.compile(r'^\d+(,\d+)?c\d+(,\d+)?$')
wasregex = re.compile(r'^< ')
nowregex = re.compile(r'^> ')
wsregex = re.compile(r'\s+')
sepregex = re.compile(r'^---$')
emptyregex = re.compile(r'^$')


def islinenum(line):
    return lineregex.match(line)


def iswas(line):
    return wasregex.match(line)


def isnow(line):
    return nowregex.match(line)


def isseparator(line):
    return sepregex.match(line)


def isempty(line):
    return emptyregex.match(line)


def fail(content):
    sys.stdout.write(content)
    sys.exit(1)


def positiveish(n):
    return n == ' ' or n == '+'


def different_sign(w, n):
    return w != n and (w == '-' or n == '-')


def bothzero(was, now):
    w, n = abs(float(was)), abs(float(now))
    return w == n and w == 0. and n == 0.


def checkfield(was, now):
    if was == now:
        return True
    elif len(was) == 0 or len(now) == 0:
        return False
    elif different_sign(was[0], now[0]) and bothzero(was, now):
        return True
    else:
        return False


def checkblock(inblock, outblock):
    if len(inblock) != len(outblock):
        return False

    for (was, now) in zip(inblock, outblock):
        for (wasfield, nowfield) in zip(wsregex.split(was), wsregex.split(now)):
            if not checkfield(wasfield, nowfield):
                return False

    return True


def main():
    content = sys.stdin.read()

    wasblock = []
    nowblock = []
    for line in content.split('\n'):
        if islinenum(line):
            if not checkblock(wasblock, nowblock):
                fail(content)
            wasblock = []
            nowblock = []
        elif iswas(line):
            wasblock.append(line[1:])
        elif isnow(line):
            nowblock.append(line[1:])
        elif isseparator(line) or isempty(line):
            continue
        else:
            print(line)
            fail('error: ill-formed diff output\n')


if __name__ == '__main__':
    main()
