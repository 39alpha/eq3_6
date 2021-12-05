import re
import sys
from math import fabs

changeregex = re.compile(r'^\d+(,\d+)?[cd]\d+(,\d+)?$')
wasregex = re.compile(r'^< ')
nowregex = re.compile(r'^> ')
fieldregex = re.compile(r'[\s\|,]+')
sepregex = re.compile(r'^---$')
emptyregex = re.compile(r'^$')


def ischange(line):
    return changeregex.match(line)


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
        content = '\n'.join(content)
    sys.stdout.write(content)
    sys.exit(1)


def positiveish(n):
    return n == ' ' or n == '+'


def isclose(was, now, rtol=1e-3, atol=0.0):
    try:
        w, n = float(was), float(now)
        diff = fabs(w - n)
        cutoff = max(1e-3 * max(fabs(w), fabs(n)), atol)
        message = 'abs({} - {}) = {} > {}'.format(w, n, diff, cutoff)
        return diff <= cutoff, message
    except ValueError:
        return False, '{} != {}'.format(was, now)


def checkfield(was, now):
    if was == now:
        return True, ''
    elif len(was) == 0 or len(now) == 0:
        return False, '{} != {}'.format(was, now)
    else:
        return isclose(was, now)


def checkblock(inblock, outblock, errors):
    if len(inblock) != len(outblock):
        return False

    for (was, now) in zip(inblock, outblock):
        wasfields, nowfields = fieldregex.split(was), fieldregex.split(now)
        for (wasfield, nowfield) in zip(wasfields, nowfields):
            valid, message = checkfield(wasfield, nowfield)
            if not valid:
                message = '< {}\n> {}Error: {}\n\n'.format(was, now, message)
                errors.append(message)
                return False

    return True


def main():
    wasblock = []
    nowblock = []
    errors = []

    for line in sys.stdin.readlines():
        if ischange(line):
            checkblock(wasblock, nowblock, errors)
            wasblock = []
            nowblock = []
        elif iswas(line):
            wasblock.append(line[1:].strip())
        elif isnow(line):
            nowblock.append(line[1:].strip())
        elif isseparator(line) or isempty(line):
            continue
        else:
            print(line)
            fail('error: ill-formed diff output\n')

    checkblock(wasblock, nowblock, errors)

    if len(errors) != 0:
        fail(errors)


if __name__ == '__main__':
    main()
