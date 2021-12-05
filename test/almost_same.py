import re
import sys
from math import fabs

lineregex = re.compile(r'^\d+(,\d+)?c\d+(,\d+)?$')
wasregex = re.compile(r'^< ')
nowregex = re.compile(r'^> ')
fieldregex = re.compile(r'[\s\|]+')
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
    if isinstance(content, list):
        content = content.join('\n')
    sys.stdout.write(content)
    sys.exit(1)


def positiveish(n):
    return n == ' ' or n == '+'


def isclose(was, now):
    try:
        w, n = float(was), float(now)
        return fabs(w - n) <= 1e-3 * max(fabs(w), fabs(n))
    except ValueError:
        return False


def checkfield(was, now):
    if was == now:
        return True
    elif len(was) == 0 or len(now) == 0:
        return False
    else:
        return isclose(was, now)


def checkblock(inblock, outblock, errors):
    if len(inblock) != len(outblock):
        return False

    for (was, now) in zip(inblock, outblock):
        wasfields, nowfields = fieldregex.split(was), fieldregex.split(now)
        for (wasfield, nowfield) in zip(wasfields, nowfields):
            if not checkfield(wasfield, nowfield):
                errors.append('<{}\n---\n>{}\n'.format(was, now))
                return False

    return True


def main():
    content = sys.stdin.read()

    wasblock = []
    nowblock = []
    errors = []
    for line in content.split('\n'):
        if islinenum(line):
            if not checkblock(wasblock, nowblock, errors):
                fail(errors)
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
