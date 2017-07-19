def n2s(n, b):
    digits = []
    while n > 0:
        digits.insert(0, str(n % b))
        n = n / b
    return ''.join(digits)


def minus(x, y, b):
    x_int = int(x, b)
    y_int = int(y, b)
    return n2s(x_int - y_int, b).zfill(len(x))


def answer(n, b):
    ids = []
    idx = 0
    while n not in ids:
        ids.append(n)
        x = ''.join(sorted(n, reverse=True))
        y = x[::-1]
        n = minus(x, y, b)
        idx += 1

    # one more time
    prev = n
    x = ''.join(sorted(n, reverse=True))
    y = x[::-1]
    n = minus(x, y, b)

    if prev == n:
        return 1
    else:
        return idx - ids.index(prev)


# test
# print answer('210022', 3)
# print answer('1211', 10)
print answer('1011', 10)